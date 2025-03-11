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

	//T = [0.5, 2, 4 ,8, 16]; rho = [0.8062   0.4225   0.1785   0.03186   0.001015]
	const double gTf = 16;
	const double gRho = 0.001015;
	const double gTol = 1e-6;
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	cout << "Running Max Frequency finite horizon Solver for constitutive killers: " << endl;


	//Build the main class of our PDE solver
	PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gTf, gRho);

	//Calling the main PDE solver
	//myExample.Main_Solver(); //Finite-horizon problem solver only 
	//myExample.Inf_limit_Solver(); //Infinite limit computation with boost-up
	myExample.Inf_limit_Solver_OutputIter(); //Infinite limit computation with one dilution after another

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
