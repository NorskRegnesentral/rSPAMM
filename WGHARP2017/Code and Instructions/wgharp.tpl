// File: wgharp.tpl    WGHARP assessment model used at Archangelsk meeting + updated with 
//			prior distributions on model parameters
// Date: 29.11.04
// Author: Hans Julius Skaug (skaug@imr.no)
// Modified: Tor Arne Øigård in 2011 (toro@imr.no)

DATA_SECTION

  // Fixed parameters in the model
  init_int	A			// Max-age in model
  
  // Parameter bounds and active switch
  init_int    K_act                     // Active switch for K 
  init_number L_K                       // Lower bound on K
  init_number U_K                       // Upper bound on K
  init_int    M_act			// Active switch for M
  init_number L_M                      	// Lower bound on M
  init_number U_M                      	// Upper bound on M
  init_int    M0_act			// Active switch for M0
  init_number L_M0                   	// Lower bound on M0
  init_number U_M0                      // Upper bound on M0


  // Prediction under given quota regime
  !!USER_CODE ad_comm::change_datafile_name("wgharp.quo");
  init_int n_pred                       // Number of different quota rules
  init_vector quotas(0,1)       	// Quotas per years (pups,adults)

  // Pup and adult catch
  !!USER_CODE ad_comm::change_datafile_name("wgharp.cat");
  init_int n                    // Number of years with catch data
  init_matrix Ctmp(1,n,-1,1)    // Matrix with columns (year,pups,adults) 
  vector years(1,n)		// Years with data
  matrix C(0,n+n_pred,-1,1)     // Matrix with columns (year,pups,adults)
  
  // P-matrix - maturity curve
  !!USER_CODE ad_comm::change_datafile_name("wgharp.pma");
  init_matrix Ptmp(1,n,1,A)		// Proportion of females mature at age
  matrix P(1,n+n_pred,1,A)
  
  // Time varying fecundity rates
  !!USER_CODE ad_comm::change_datafile_name("wgharp.f");
  init_vector ftmp(1,n)       	// Quotas per years (pups,adults)
  vector f(1,n+n_pred)
  
  // Pup production estimates
  !!USER_CODE ad_comm::change_datafile_name("wgharp.est");
  init_int n_est                    		// Number of years with estimates
  init_matrix production_file(1,n_est,1,3) 	// Matrix with 3 cols (year,est,cv) 
  matrix pup_prod(1,n_est,1,3)			// Matrix with 3 cols (index,est,sd) 

  // Prior information
  !!USER_CODE ad_comm::change_datafile_name("wgharp.pri");
  init_matrix priors(1,3,1,2) 	// Matrix of priors (prior_mean,prior_sd)
  vector mub(1,3)		// Prior mean of b
  vector sdb(1,3)		// Prior sd of b

 
PARAMETER_SECTION
    
  init_bounded_number K(L_K,U_K,K_act)   	// Starting population
  init_bounded_number M(L_M,U_M,M_act)		// Natural mortality for adults
  init_bounded_number M0(L_M0,U_M0,M0_act)	// Natural mortality for pups
  
  vector b(1,3)					// Concatination of parameters
  vector Ntot(1,n+1)
  matrix Nout(1,n+n_pred,1,3)
  matrix N(0,n+n_pred,1,A)                      // Population matrix
  sdreport_vector N0(0,n+n_pred)                         // 0-group
  sdreport_vector N1(0,n+n_pred)                // 1+ populasjon
  sdreport_number N0_2003                	// 1+ populasjon today
  sdreport_number D				// N_(y+10)/N_y eqv.
  sdreport_number Dnew				// N_pred/Nmax


  // General parameters
  number tmp
  number em
  number em0
  number Ft
  number Ntotmax
  

  objective_function_value l

PRELIMINARY_CALCS_SECTION

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
 
  
GLOBALS_SECTION
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

  

PROCEDURE_SECTION
  
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


REPORT_SECTION
  report << setprecision(4);

   report << Nout << endl;
  //report << N << endl;

TOP_OF_MAIN_SECTION
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  arrmblsize = 1000000L;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);


FUNCTION write_mcchains
  mcmc_par << N0 << N1 << endl;
