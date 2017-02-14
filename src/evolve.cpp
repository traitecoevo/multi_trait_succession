#include <cstdio>
#include <string>
#include <cstring>

#include "arg_parser.h"
#include "AdaptiveDynamics_nD.h"
#include "AdaptiveDynamics_2D.h"
#include "EBT_Metapopulation.h"

#include "MatrixCPP.h"
#include "Collect.h"

void reset_array23(double Array[2][3],double X00, double X01, double X02, double X10, double X11, double X12);
void reset_array22(double Array[2][2], double X00, double X01, double X10, double X11);

namespace
{
	const char * invocation_name = 0;
	const char * const Program_name    = "E";
}

void show_error( const char * msg, const int errcode = 0, const bool help = false ) throw() {
	if( msg && msg[0] != 0 )
	{
		std::fprintf( stderr, "%s: %s", Program_name, msg );
		if( errcode > 0 ) std::fprintf( stderr, ": %s", std::strerror( errcode ) );
		std::fprintf( stderr, "\n" );
	}
	if( help && invocation_name && invocation_name[0] != 0 )
		std::fprintf( stderr, "Try `%s --help' for more information.\n", invocation_name );
}

void internal_error( const char * msg ) throw() {
	char buf[80];
	std::snprintf( buf, sizeof( buf ), "internal error: %s.\n", msg );
	show_error( buf );

	exit(3);
}

using namespace std;

//******************************************************************************************************
int main( int argc, char * argv[] ) throw() {
// Argument parser
	invocation_name = argv[0];
	static const Arg_parser::Option options[] =
	{
			{ 'd', "timedist", Arg_parser::yes   },
			{ 'p', "productivity", Arg_parser::yes   },
			{ 'x', "parameter", Arg_parser::yes   },
			{ 'a', "multiplyBy", Arg_parser::yes   },
			{ 't', "trait", Arg_parser::yes   },
			{ 'v', "value", Arg_parser::yes   },
			{ 'f', "path_to_community", Arg_parser::yes   },
			{ 'r', "resolution", Arg_parser::yes   },
			{ 's', "step", Arg_parser::yes   },
			{ 'L', "print Landscape", Arg_parser::yes   },
			{ 'V', "print Viability", Arg_parser::yes   },
			{ 'P', "print profile", Arg_parser::yes   },

	};
	Arg_parser parser( argc, argv, options );
	if( parser.error().size() )	{ show_error( parser.error().c_str(), 0, true ); return 1; } // bad option

	int i, numRes, maxRes =EBT_MAX_N0_SPECIES;

// Setup fitness solver (ebt)
	sim_params p;
	EBT_Metapopulation ebt;
	ebt.setup(&p);
// Setup PopulationDynamics
	PopulationDynamics popDyn;
	popDyn.setup(&ebt);
// Setup CANONICAL EQN
	AdaptiveDynamics_nD ADN;
	ADN.setup(&popDyn);
// Setup 2D analysis
	AdaptiveDynamics_2D AD2;
	AD2.setup(&popDyn);
// Allocate memory for other Stuff
	double** Res =  M2Dd_alloc(maxRes, TRAIT_DIM);
	double* Mut =  Vd_alloc(TRAIT_DIM);
	double* N =  Vd_alloc(maxRes);
	string name;
	vector<int> traitList(2,0); traitList[1]=1;  // Used when printing trait values to screen

// Display options
	ebt.print=1;
	AD2.setDisplayWorkingDetails(1);
	popDyn.setDisplayWorkingDetails(1);
	ADN.setDisplayWorkingDetails(1);

// Set criteria for population dynamics
	popDyn.maxPopulationIterations =1;

// Set default parameters
	double disturbanceInterval = 60;
	double siteProductivity=1;
	int parameter=0;
	double multiplyBy=1;
	int whichTrait =2;
	double mutationVariance = 1e-3;
	double mutationRate=0.1;
	double ImmigrationRate = 1.0;
	int requiredSteps = 1005;
	int printFrequency = 20;
	double traitRange[2][3]; // trait vector
	reset_array23(traitRange, -2.5, 1.5, 0.005, log10(0.1), log10(100), 0.005);
	double value =1;
	string path_to_community = "";
	int resolution = 0;
	bool Viability=0;
	bool PopulationProfile =0;
	bool Landscape = 0;
	int step=0;

// PROCESS PASSED ARGUMENTS - these override others in file
	for( int argind = 0; argind < parser.arguments(); ++argind )
	{
		const int code = parser.code( argind );
		if( !code ) break;					// no more options
		switch( code )
		{
			case 'p': siteProductivity= atof(parser.argument( argind ).c_str()); break;
			case 'd': disturbanceInterval= atof(parser.argument( argind ).c_str()); break;
			case 'x': parameter= atoi(parser.argument( argind ).c_str()); break;
			case 'a': multiplyBy= atof(parser.argument( argind ).c_str()); break;
			case 't': whichTrait= atoi(parser.argument( argind ).c_str()); break;
			case 'v': value= atof(parser.argument( argind ).c_str()); break;
			case 'f': path_to_community= parser.argument( argind ).c_str(); break;
			case 'r': resolution= atoi(parser.argument( argind ).c_str()); break;
			case 's': step= atoi(parser.argument( argind ).c_str()); break;
			case 'V': Viability= atoi(parser.argument( argind ).c_str()); break;
			case 'P': PopulationProfile= atoi(parser.argument( argind ).c_str()); break;
			case 'L': Landscape= atoi(parser.argument( argind ).c_str()); break;
			default : internal_error( "uncaught option" );
		}
	} // end process options

// Resolution options
	ebt.set_max_no_cohorts(4000);
	p.ENV_DIM = 30;
	ebt.set_eps(5, 1e-4);	// Error control Used in adaptive sampling of environment
	ebt.set_eps(6, 1e-6);   // Minimum step size in resident ODE
	ebt.set_eps_cohort_splt(50, 1000);
	ebt.set_eps(1, 1e-3);   // Error control on resident ode
	ebt.set_max_cohort_fitness(0.05);

	if(resolution > 0){ // High resolution
		ebt.set_max_cohort_fitness(0.01);
		ebt.set_eps(1, 5e-3);
	}

	// SELECT TRAIT TO EVOLVE & SET OPTIONS
	vector<double> covarMatrixDiagonal(4,0);
	if(whichTrait ==0){
		covarMatrixDiagonal[0]=mutationVariance;
	} else if(whichTrait ==1){
		covarMatrixDiagonal[1]=mutationVariance;
	} else if(whichTrait ==2){
		covarMatrixDiagonal[0]=covarMatrixDiagonal[1]=mutationVariance;
		reset_array23(traitRange, -2.5,1.5, 0.05, log10(0.1), log10(100), 0.05);
		cout << "changing array size" <<endl;
		requiredSteps = 7505;
		printFrequency = 200;
	} else {
		cout << "bad trait selection " << whichTrait <<endl;
		exit(1);
	}

if(path_to_community ==""){ // New simulation, set pars and determine directory from pars
// GO TO OUTPUT DIRECTORY
	go_to_dir("output");
	go_to_dir("data");
	go_to_dir("[" + stringify(disturbanceInterval) + "," + stringify(siteProductivity) + "]");

	if(whichTrait ==0){
		go_to_dir("lcc");
		go_to_dir("[" + stringify(value) + "]");
	} else if(whichTrait ==1){
		go_to_dir("hsp");
		go_to_dir("[" + stringify(value) + "]");
	} else if(whichTrait ==2){
		go_to_dir("2trait");
	} else {
		cout << "bad trait selection " << whichTrait <<endl;
		exit(1);
	}


// ADJUST OTHER PARAMS IF REQUESTED
	if(multiplyBy ==1.0){  // Use baseline parameter values
		go_to_dir("base");
	}
	else{  // change `parameter` by amount `multiplyBy`

		// Create a vector of pointers to different parameters
		int Npars =18;
		vector<double*> Pars(Npars, NULL);
		vector<string> names(Npars, "");
		int cc=0;
		Pars[cc]= &p.seed_mass;	names[cc] = "seed";
		Pars[++cc]= &p.wood_dens;	names[cc] = "wood_dens";// Leaf area per sapwood area
		Pars[++cc]= &p.theta;		names[cc] = "theta";			// Leaf area per sapwood area
		Pars[++cc]= &p.a1;			names[cc] = "a1";					// Height - leaf area scaling
		Pars[++cc]= &p.a3;			names[cc] = "a3";					// Root leaf scaling
		Pars[++cc]= &p.a4;			names[cc] = "a4";					// LMA - LL scaling
		Pars[++cc]= &p.c_Rs;		names[cc] = "c_Rs";				// Sapwood respiration
		Pars[++cc]= &p.c_Rr;		names[cc] = "c_Rr";				// Root respiration
		Pars[++cc]= &p.c_Rl;		names[cc] = "c_Rl";				// Leaf respiration
		Pars[++cc]= &p.c_ext;		names[cc] = "c_ext";    	// Light extinction coefficient
		Pars[++cc]= &p.c_acc;		names[cc] = "c_acc";			// Accessory cost of reproduction - multiplication factor
		Pars[++cc]= &p.Pi_0;		names[cc] = "pi_0";				// Survival during dispersal
		Pars[++cc]= &p.c_d0;		names[cc] = "c_d0";     	// Baseline structural mortality rate
		Pars[++cc]= &p.c_d2;		names[cc] = "c_d2";     	// Baseline for growth mortality rate
		Pars[++cc]= &p.c_d3;		names[cc] = "c_d3";				// Coefficient for dry mass production in mortality function
		Pars[++cc]= &p.Eta;		names[cc] = "eta";      		// Canopy shape parameters
		Pars[++cc]= &p.c_r1;		names[cc] = "c_r1";      	// Reproductive allocation - max
		Pars[++cc]= &p.c_r2;		names[cc] = "c_r2";      	// Reproductive allocation - steepness

		if( (parameter<0) | (parameter>= Npars)){
			cout << "bad parameter selection " << parameter<<endl;
			exit(1);
		}

		// Change parameter value
		cout << "Adjusting " << names[parameter] << ": " << *Pars[parameter] << "\t";
		*Pars[parameter]=(*Pars[parameter])*multiplyBy;
			go_to_dir("["+names[parameter] + "," + stringify(multiplyBy) + "]");
		cout<<*Pars[parameter] << "\n";
	}

// Set disturbance interval and productivity
	p.log_mean_disturbance_interval=log10(disturbanceInterval);
	p.c_p1=siteProductivity*p.c_p1;

// print parameters file
	ofstream ofp("params.m");
	p.print_sim_params(ofp);
	ofp.close();

// Set baseline trait values
	for(i = 0; i<maxRes; i++){
		Res[i][0] = 0;  Res[i][1] = 0;
		Res[i][2]= log10(608.0); Res[i][3] = log10(3.8E-05);

		// Choose specific value of other trait when doing 1D evolution
		if(whichTrait==0)
			Res[i][1] = value; // Set height to specified value
		if(whichTrait==1)
			Res[i][0] = value; // Set LCC to specified value
	}

	if(Viability){
		double traitRangeWide[2][3]; // trait vector
		reset_array23(traitRangeWide, -3.0, 2, 0.005, log10(0.5), log10(100), 0.005);
		AD2.followTraitViabilityContour2D(Res[0], traitRangeWide, "Viability");
	}

	popDyn.S(Res[0], Res[0],1);

	// Print landscape for time zero with empty environment
	if(whichTrait < 2)
		ADN.outputFitnessLandscape1D(Res[0], 0, Res[0], whichTrait, traitRange[whichTrait], "T" + stringify(whichTrait) + "-" + stringify(0), 0);
	else
		ADN.outputFitnessLandscape2D(Res[0], 0, traitRange, "T-" + stringify(0), 0);

	// Evolve traits with stochastic model
	numRes=1;
	ADN.runStochasticModelNonEquil(Res[0], numRes, covarMatrixDiagonal, traitRange, mutationRate, ImmigrationRate, 0 , requiredSteps, printFrequency, "Stoch.txt");

} else{ // Use existing simulation

	go_to_dir(path_to_community);

	// Load parameters file
	string paramsFile = "params.m";
	if(FileExist(paramsFile.c_str()))
		p.load(paramsFile);
	else{
		cout << "Cannot find parameters file" << endl;
		exit(1);
	}

	// Set baseline trait values
	for(i = 0; i<maxRes; i++){
		Res[i][0] = -1;  Res[i][1] = 1;
		Res[i][2]= log10(608.0); Res[i][3] = log10(3.8E-05);
	}

	popDyn.S(Res[0], Res[0],1);

	double Time = numRes = 0;

	if(Viability){
		double traitRangeWide[2][3]; // trait vector
		reset_array23(traitRangeWide, -3.0, 2, 0.005, log10(0.5), log10(100), 0.005);
		AD2.followTraitViabilityContour2D(Res[0], traitRangeWide, "Viability");
	}

	if(step > 0){  // load existing simulation at step
		Time = step;
		ADN.inputTraitsFromStochFile(Res[0], numRes, Time, "Stoch.txt", Time);
		cout<<numRes<<endl;
		for(i = 0; i<numRes; i++)
			for(int j = 0; j<4; j++)
				cout<<Res[i][j] << "\t";
	}
	else
		numRes=0;

	if(PopulationProfile)
		popDyn.printPopulationProfile(Res[0], numRes, "T-" + stringify(Time), 1);

	if(Landscape){
		if(whichTrait < 2)
			ADN.outputFitnessLandscape1D(Res[0], numRes, Res[0], whichTrait, traitRange[whichTrait], "T" + stringify(whichTrait) + "-" + stringify(Time), 0);
		else
			ADN.outputFitnessLandscape2D(Res[0], numRes, traitRange, "T-" + stringify(Time), 0);
	}
}

	delete []N; delete []Mut; M2DdFree(Res);
	cout << "END" << endl;
	return 0;
}


void reset_array22(double Array[2][2], double X00, double X01, double X10, double X11) {
	Array[0][0] = X00;	Array[0][1] = X01;
	Array[1][0] = X10;	Array[1][1] = X11;
}

void reset_array23(double Array[2][3], double X00, double X01, double X02, double X10, double X11, double X12) {
	Array[0][0] = X00;	Array[0][1] = X01;	Array[0][2] = X02;
	Array[1][0] = X10;	Array[1][1] = X11;  Array[1][2] = X12;
}


string return_dir_name(int dir) {if(dir<10)
	return("0"+ stringify(dir));
else
	return(stringify(dir));
}
