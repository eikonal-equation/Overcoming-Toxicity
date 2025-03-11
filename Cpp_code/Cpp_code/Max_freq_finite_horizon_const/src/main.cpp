/*==============================================================================
 * File: main.cpp
 *
 * Author: MingYi Wang
 *
 * Description: This file initializes all the global variables and executes the
 * the corresponding example from the command line.
 *
 *============================================================================*/

 //-----------------------Project specific header files---------------------------
#include "CSolver.h"

//----------------------------Libraries------------------------------------------
#include <chrono>

using namespace std;

int main()
{
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 16;
	const double gMax_population = 1;
	const double gMax_frequency = 1;
	const double gEpsilon = 0.2;
	const double gRks = 0.85;
	const double gGamma = 1.0;
	//const double gLamb = 0.5;
	//const double gTf = 1.0 / gLamb;
	const double gTf = 3;
	const double gRho = 0.65;
	const double gTol = 1e-6;
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	cout << "Running Max Frequency finite horizon Solver for constitutive killers: " << endl;

	//array<double,5> lamb_list{0.8,0.85,0.9,0.95,1.0};
	//array<double,5> lamb_list{0.7,0.825,0.875,0.925,0.975};
	//array<double,10> lamb_list{0.7,0.8,0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0};
	// array<double,5> rho_list{0.5,0.55,0.6,0.65,0.7};
	//array<double,8> rho_list{0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59};
	//array<double,8> rho_list{0.525,0.535,0.545,0.555,0.565,0.575,0.585,0.595};

	// array<double,19> rho_list{0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.60,0.65,0.70};

	// array<double,1> lamb_list{0.75};
	// //array<double,1> rho_list{0.65};
	// for (int ii = 0; ii < 1; ii++){
	// 	for (int jj = 0; jj < 19; jj++){
	// 		double gLamb = lamb_list[ii];
	// 		double gRho = rho_list[jj];
	// 		double gTf = 1.0/gLamb;
	// 		//Build the main class of our PDE solver
	// 		PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gTf, gRho);

	// 		//Calling the main PDE solver
	// 		//myExample.Main_Solver();
	// 		myExample.Inf_limit_Solver();

	// 		// Writing Grid Parameters to file
	// 		std::string file_label = std::to_string(gRefinement_factor);
	// 		myExample.writeDomainToFile("output/Max_freq_finite_horizon_const");
	// 	}
	// }
	
	
	
	
	//Build the main class of our PDE solver
	PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gTf, gRho);

	//Calling the main PDE solver
	//myExample.Main_Solver();
	myExample.Inf_limit_Solver();


	//Writing Grid Parameters to file
	std::string file_label = std::to_string(gRefinement_factor);
	myExample.writeDomainToFile("Max_freq_finite_horizon_const");
	
	// std::string file_label = std::to_string(gRefinement_factor);
	// myExample.writeDomainToFile("Max_freq_finite_horizon_single");


	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
