DATA_SECTION
	//create a command line option
	//-sim seed to simulate data
	int sim;
	int rseed;
	LOCAL_CALCS
		sim=0; //default 0 read in data
		rseed=0; //default rseed value
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1) //option match is a true false
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
		
	END_CALCS
	// read in data and control file from base pm.dat file
	init_adstring datafile;
	init_adstring ctlfile;
	//open data file
	!! ad_comm::change_datafile_name(datafile);
	init_int syr;
	init_int nyr; 
	init_int nage;
	
	//growth pars note not used in the production model but fun to read in
	init_number linf;
	init_number vonbk;
	init_number a;
	init_number b;
	
	//time series data
	init_ivector iyr(syr,nyr);  //integer type vector
	init_vector yt(syr,nyr); //observed relative abundace
	init_vector ct(syr,nyr); // observed catch
	init_matrix A(syr,nyr,1,nage); //observed age composition but not used here
	//end of file.
	init_int eof;
	//!!cout<<yt<<endl;
	!!cout<<""<<endl;

	//Warn user if eof != 999 and quit.
	LOC_CALCS
		if(eof != 999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			exit(1);
		}
	END_CALCS
	//Open control file
	!! ad_comm::change_datafile_name(ctlfile);
	init_int npar;
	init_vector ipar(1,npar);
	init_int eofc;
	LOC_CALCS
		if(eofc != 999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			exit(1);
		}
	END_CALCS
	int iter;
	//!!cout<<yt<<endl;
	//!!ad_exit(1);	
	
PARAMETER_SECTION
	init_number log_k;  //log of carrying capacity.
	init_bounded_number r(0.001,1.000);
	init_number log_q;	//catchability coefficient.
	init_number log_precsig(3);//precision process error 3 phase estimation
	init_number log_prectau(3);//precision obs error
	init_bounded_vector xt(syr,nyr,-15.,15.,2);//Begin this in phase 2

	LOCAL_CALCS
		log_k = log(ipar(1));
		r = ipar(2);
		log_q = log(ipar(3));
		log_precsig=log(1./ipar(4));
		log_prectau=log(1./ipar(5));
	END_CALCS
	
	objective_function_value f;
	
	number k;
	number q;
	number fpen;
	number sig;
	number tau;

	//proportion of biomass relative to k
	vector pt(syr,nyr+1);
	//vector bt(syr,nyr+1);
	vector it(syr,nyr);	//predicted cpue
	vector nu(syr,nyr);	//cpue residuals

	sdreport_number depletion;
	
PRELIMINARY_CALCS_SECTION
	
	if(sim)
	{
		run_simulation_model();
	}
	
PROCEDURE_SECTION
	//*** Main function calls  ***//
	// 1) population_dynamics
	population_dynamics();
	// 2) observation_model
	observation_model();	
	// 3) calc_objective_function
	calc_objective_function();
	if(mceval_phase()) mcmc_output();
	//****************************//

FUNCTION population_dynamics
	{
	//transorm parameters
	k = mfexp(log_k);
	//initialize population at unfished
	pt(syr) = 1.0;
	//bt(syr) = k;
	fpen = 0.0;
	int i;
	for(i=syr; i<=nyr; i++) //for (i in syr:nyr) R code
	{
		pt(i+1)=posfun((pt(i)+r*pt(i)*(1.-pt(i))-ct(i)/k)*mfexp(xt(i)),0.001,fpen);
		//bt(i+1)=posfun((bt(i)+r*bt(i)*(1.-bt(i)/k)-ct(i))*mfexp(xt(i)),0.001,fpen);
	}
	depletion = pt(nyr);
	//cout<<pt<<endl;
	//exit(1);
	}
	
FUNCTION observation_model
	{
		q = mfexp(log_q);
		it = q*k*pt(iyr); //NB iyr is an ivector
		nu = log(yt) - log(it);
	}

FUNCTION calc_objective_function
	{
		dvar_vector likevec(1,8);
		//likelihoods
		//obs residuals
		tau = sqrt(1./mfexp(log_prectau));
		likevec(1) = dnorm(nu,tau);
		// process residuals
		sig=sqrt(1./mfexp(log_precsig));
		likevec(2)=dnorm(xt,sig);
		//penlty for biomass going to 0. THis is necessary because of the stock reduction approach
		likevec(3)=1.e5*fpen;

		//priors
		//lognormal priors for r,k
		likevec(4)=dlnorm(k,8,0.25);
		likevec(5)=dlnorm(r,-1.38,0.51);
		//uniform prior on q
		likevec(6)=-log(q);
		likevec(7)=dgamma(1.0/(tau*tau),1.71,0.0086);
		likevec(8)=dgamma(1.0/(sig*sig),3.79,0.0102);;
		
		// total objective function
		f=sum(likevec);
		//add penalties
		//if(fpen>0)cout<<"Fpen = "<<fpen<<endl;
	}

FUNCTION mcmc_output

	// code to do forward projections
	// assume fishing at fmsy for 20 yrs
	double tac; // creates a number called tac that is not differentiable
	int i; //creates integer i for looping
	dvector newpt(nyr+1,nyr+22); //creates a vector for future depletion levels not differentiable
	if(iter==0)
	{
		ofstream ofs("projection.mcmc");
	}
	iter++;
	newpt(nyr+1)=value(pt(nyr+1)); // initialize with the last depletion value note the type transformation
	for(i=nyr+1; i<=nyr+21; i++) // project forward
	{
		tac=value(r/2.*k)*newpt(i); // set tac for a yr as fmsy* current biomass
		newpt(i+1)=newpt(i)+value(r)*newpt(i)*(1.-newpt(i))-tac/value(k);
	}
	ofstream ofs("projection.mcmc",ios::app);
	ofs<<newpt<<endl;

	/*if(iter==0)
	{
		ofstream ofs("lpars.mcmc");
		ofs<<"r\t k\t q\t precsig\t prectau\t"<<endl;
	}
	iter++;
	ofs("lpars.mcmc",ios::app);
	ofs<<r<<"\t"<<k<<"\t"<<q<<"\t"<<mfexp(log_precsig)<<"\t"<<mfexp(log_prectau)<<endl;
	*/

FUNCTION run_simulation_model
	
	random_number_generator rng(rseed);
	dvector tmp_nu(syr,nyr);
	dvector tmp_xt(syr,nyr);
	tmp_nu.fill_randn(rng);
	tmp_xt.fill_randn(rng);

	tau = sqrt(1./mfexp(log_prectau));
	sig = sqrt(1./mfexp(log_precsig));

	xt=tmp_xt*sig;
	
	population_dynamics();
	observation_model();
	nu=tmp_nu*tau;
	yt=value(elem_prod(it,exp(nu)));
	cout<<yt<<endl;
	cout<<"sim"<<endl;


REPORT_SECTION
	// capitalied report associates the names with the values in the report
	REPORT(fpen);
	REPORT(pt);
	//extracting values from dvar_vectors
	//placing into a dvector.
	dvector bt = value(pt(iyr)*k);
	REPORT(bt);
	REPORT(iyr);
	REPORT(nu);
	REPORT(sig);
	REPORT(tau);
	REPORT(r);
	REPORT(k);
	REPORT(q);
	REPORT(ct);
	REPORT(yt);
	REPORT(it);
	double fmsy=value(r/2.);
	double msy=value(r*k/4.);
	double tac=value(fmsy*pt(nyr+1)*k);
	REPORT(fmsy);
	REPORT(msy);
	REPORT(tac);
	
	
TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) \
	report << #object "\n" \
	<< object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>
	//#include <stats.cxx>

	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

