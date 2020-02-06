#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <fvar.hpp>
  ofstream mcmc_par("MCMC_chains.txt");
  dvariable posfun(const dvariable&x,const double eps)
  {
    if (x>=eps) 
    {
      return x;
    }
    else
    {
      dvariable y=1.0-x/eps;
      return eps*(1./(1.+y+y*y));
    }
  }
  
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <wgharp.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  A.allocate("A");
  K_act.allocate("K_act");
  L_K.allocate("L_K");
  U_K.allocate("U_K");
  M_act.allocate("M_act");
  L_M.allocate("L_M");
  U_M.allocate("U_M");
  M0_act.allocate("M0_act");
  L_M0.allocate("L_M0");
  U_M0.allocate("U_M0");
 ad_comm::change_datafile_name("wgharp.quo");
  n_pred.allocate("n_pred");
  quotas.allocate(0,1,"quotas");
 ad_comm::change_datafile_name("wgharp.cat");
  n.allocate("n");
  Ctmp.allocate(1,n,-1,1,"Ctmp");
  years.allocate(1,n);
  C.allocate(0,n+n_pred,-1,1);
 ad_comm::change_datafile_name("wgharp.pma");
  Ptmp.allocate(1,n,1,A,"Ptmp");
  P.allocate(1,n+n_pred,1,A);
 ad_comm::change_datafile_name("wgharp.f");
  ftmp.allocate(1,n,"ftmp");
  f.allocate(1,n+n_pred);
 ad_comm::change_datafile_name("wgharp.est");
  n_est.allocate("n_est");
  production_file.allocate(1,n_est,1,3,"production_file");
  pup_prod.allocate(1,n_est,1,3);
 ad_comm::change_datafile_name("wgharp.pri");
  priors.allocate(1,3,1,2,"priors");
  mub.allocate(1,3);
  sdb.allocate(1,3);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  K.allocate(L_K,U_K,K_act,"K");
  M.allocate(L_M,U_M,M_act,"M");
  M0.allocate(L_M0,U_M0,M0_act,"M0");
  b.allocate(1,3,"b");
  #ifndef NO_AD_INITIALIZE
    b.initialize();
  #endif
  Ntot.allocate(1,n+1,"Ntot");
  #ifndef NO_AD_INITIALIZE
    Ntot.initialize();
  #endif
  Nout.allocate(1,n+n_pred,1,3,"Nout");
  #ifndef NO_AD_INITIALIZE
    Nout.initialize();
  #endif
  N.allocate(0,n+n_pred,1,A,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  N0.allocate(0,n+n_pred,"N0");
  N1.allocate(0,n+n_pred,"N1");
  N0_2003.allocate("N0_2003");
  D.allocate("D");
  Dnew.allocate("Dnew");
  tmp.allocate("tmp");
  #ifndef NO_AD_INITIALIZE
  tmp.initialize();
  #endif
  em.allocate("em");
  #ifndef NO_AD_INITIALIZE
  em.initialize();
  #endif
  em0.allocate("em0");
  #ifndef NO_AD_INITIALIZE
  em0.initialize();
  #endif
  Ft.allocate("Ft");
  #ifndef NO_AD_INITIALIZE
  Ft.initialize();
  #endif
  Ntotmax.allocate("Ntotmax");
  #ifndef NO_AD_INITIALIZE
  Ntotmax.initialize();
  #endif
  l.allocate("l");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  int i,j;
  // Ordner fangster
  C(0) = 0.0;
  for(i=1;i<=(int) n;i++)
  {
        C(i) = Ctmp(i);             // Historical catches
	f(i) = ftmp(i);
	P(i) = Ptmp(i);
 }
 
  for(i=n+1;i<=(int) n+n_pred;i++)
  {
    C(i,-1) = C(i-1,-1)+1;            // Year
    C(i,0) = quotas(0);               // Quota 0-group
    C(i,1) = quotas(1);               // Quota 1+
    f(i) = f(n);
    P(i) = Ptmp(n);
  }
  
  for(i=1;i<=(int) n;i++)
    years(i) = C(i,-1);
  // Transforms pup production estimates
  for(i=1;i<=(int) n_est;i++)
  {
    pup_prod(i,1) = production_file(i,1) - years(1) + 1;	// Year
    pup_prod(i,2) = production_file(i,2);			// Mean
    pup_prod(i,3) = pup_prod(i,2)*production_file(i,3);		// SD
  }
  // Extract priors
  mub = column(priors,1);
  sdb = column(priors,2);
  cout << "Priors:\n" << "mean = " << mub << endl << "sd   = " << sdb << endl;
 
  
}

void model_parameters::userfunction(void)
{
  l =0.0;
  int i,j,k,maxind;
  // Mortalities
  em = exp(-M);
  em0 = exp(-M0);
  N0 = 0;
  // Concatination of parameters into b
  b(1) = K;
  b(2) = M;
  b(3) = M0;
  // Initiate of N in year 0 (1945)
  for(j=1;j<=(int) A;j++)
    N(0,j) = exp(-j*M);			// Adults
  N(0,A) /= 1-em;			// Correct A+ group
  N(0) = K*N(0)/sum(N(0));		// Normalize vector to K
  N1(0) = K;			
  N0(0) = (1-em)/em0*K;			// To balance natural mortality of 1+ group
  // Calculate population trajectory
  for(i=1;i<=(int) n + n_pred;i++)
  {
    N(i,1) = (N0(i-1)-C(i-1,0))*em0;			// 0-group from last year
    for(j=2;j<=(int) A;j++)
    {
      N(i,j) = N(i-1,j-1)*(1-C(i-1,1)/N1(i-1))*em;	// Pro-rata distribution of catches
    }  
    N(i,A) += N(i-1,A)*(1-C(i-1,1)/N1(i-1))*em;      	// A+ group correction
    // Ensures that N > 0
    N1(i) = posfun(sum(N(i))-C(i,1),1) + C(i,1);	 
    N(i) = N1(i)*N(i)/sum(N(i));			// Scales up the full age structure
    // Reqruitement equation
	//Ft = 1 - (1-f(i))*pow(N1(i-1)/1500000,2.5);
    for(j=1;j<=(int) A;j++)
      N0(i) += 0.5*f(i)*P(i,j)*N(i,j);
    // Calculate the D(1+) statistic
    if(C(i,-1)==C(n+1,-1))     
    {
      D = 1/(N1(i));
      N0_2003 = N0(i);
    }
    if(C(i,-1)==C(n+16,-1)){  
      D *= (N1(i));
      Dnew = (N1(i)+N0(i));
    }    
    // Ensures that pup production is larger than pup catch
      N0(i) = posfun(N0(i)-C(i,0),1) + C(i,0);	 
  }
  Ntotmax = 0;
  for(i=1;i<=(int) (n+1);i++)
   {
    Ntot(i) = N0(i) + N1(i);
	if(Ntot(i) > Ntotmax){
		Ntotmax = Ntot(i);
		maxind = i;
	}
   }
   Dnew /= Ntot(maxind);    
  // Likelihood contribution
  l = 0.0;
  // Pup production estimates
  for(i=1;i<=(int) n_est;i++)
  {
    tmp = (pup_prod(i,2) - N0(pup_prod(i,1)))/pup_prod(i,3);
    l += -.5*tmp*tmp;
  }
  // Likelihood penalty to avoid negative N
  for(i=1;i<=(int) n;i++)
    l -= .001*(fabs(N1(i)-10000)-(N1(i)-10000))*(fabs(N1(i)-10000)-(N1(i)-10000));
  // Add contribution from prior distibutions
  l += -.5*norm2(elem_div(b-mub,sdb)) - sum(log(sdb));
  l = -l;
  cout << setprecision(4);
  for(int i=1;i<=(int) n+n_pred;i++)
  {
    Nout(i,1) = C(i,-1);
    Nout(i,2) = N0(i);  
    Nout(i,3) = N1(i);  
  }
  if(mceval_phase())   // check if the program is in an mceval phase
    write_mcchains();  // and if so write the chains using the function write_mcchains
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << setprecision(4);
   report << Nout << endl;
  //report << N << endl;
}

void model_parameters::write_mcchains(void)
{
  mcmc_par << N0 << N1 << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  arrmblsize = 1000000L;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
