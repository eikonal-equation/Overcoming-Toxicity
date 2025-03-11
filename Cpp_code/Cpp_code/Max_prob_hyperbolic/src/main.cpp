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
    setbuf(stdout, NULL);
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 16;
	const double gMax_population = 1;
	const double gMax_frequency = 1;
	const double gEpsilon = 0.2;
	const double gRks = 0.85;
	const double gGamma = 1.0;
	//const double gLamb = 1.0;
	//const double gRho = 0.65;
	const double gDefeat_thres = 0.005;
	const double gVictory_thres = 0.95;
	const double gTol = 1e-6;
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	cout << "Running Max Prob hyperbolic Solver: " << endl;
	// //Build the main class of our PDE solver
	// PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gDefeat_thres, gVictory_thres, gEpsilon, gRks, gGamma, gLamb, gRho);

	// //Calling the main PDE solver
	// myExample.Main_VPI_Solver();


	// // Writing Grid Parameters to file
	// std::string file_label = std::to_string(gRefinement_factor);
	// myExample.writeDomainToFile("Max_prob_hyperbolic");

	array<double,6> lamb_list{0.75,0.8,0.85,0.9,0.95,1.0};
	array<double,5> rho_list{0.5,0.55,0.6,0.65,0.7};
	

	for (int ii = 0; ii < 6; ii++){
		for (int jj = 0; jj < 5; jj++){
			if (ii ==0  && jj == 0){
				auto start1 = chrono::high_resolution_clock::now();
				double gLamb = lamb_list[ii];
				double gRho = rho_list[jj];
				PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gDefeat_thres, gVictory_thres, gEpsilon, gRks, gGamma, gLamb, gRho);
				// double val = myExample.Solver_with_output();
				// ult_map[ii][jj] = val;
				myExample.Main_VPI_Solver();

				// Writing Grid Parameters to file
				std::string file_label = std::to_string(gRefinement_factor);
				myExample.writeDomainToFile("Max_prob_hyperbolic");
				auto end1 = chrono::high_resolution_clock::now();
				
				std::chrono::duration<double> elapsed_seconds = end1 - start1;
				cout << "THIS runtime: " << elapsed_seconds.count() << "s\n";
				}
				
		}
	}
	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
