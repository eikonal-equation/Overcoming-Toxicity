/*=============================================================================NaN
 * Copyright (C) 2025 MingYi Wang
 *
 * This program is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
 *============================================================================*/

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
#include "WriteToFile.h"
//----------------------------Libraries------------------------------------------
#include <chrono>

using namespace std;

int main()
{
    setbuf(stdout, NULL);
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 16; //Global refinement factor with respect to the 3D (x,y,t) grid (x&y for spatial; t for time)
	const double gMax_population = 1; //Maximum nomralized population size (max y value for the spatial grid)
	const double gMax_frequency = 1; //Maximum nomralized frequency of the killer (max x value for the spatial grid)
	const double gEpsilon = 0.2; //The cost of constitutive toxin production
	const double gRks = 0.85; //Intrinsic cell division rates ratio between the killer and the sensitive
	const double gGamma = 1.0; //Rescaled killing rate
	const double gTol = 1e-6; //Global Tolerance for the convergence of the iterative solver
	const double gDefeat_thres = 0.01; //Global defeat threshold (the value of the value function below which the system is considered to be in the defeat zone)
	const double gVictory_thres = 1.0 - gDefeat_thres; //Global victory threshold (the value of the value function above which the system is considered to be in the victory zone)
	
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	const double gLamb = 1.0; // Global arrival rate of the killer cells (lamb = 1.0 means that the arrival rate is equal to the inverse growth rate of the sensitive cells)
	const double gRho = 0.65; // Global survival rate
	cout << "Running the strategically-optimal killer solver (hyperbolic BC; random dilutions): " << endl;
	//Build the main class of our PDE solver
	PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gLamb, gRho);

	//Calling the main PDE solver
	myExample.Main_VPI_Solver();

	// Writing Grid Parameters to file
	myExample.writeDomainToFile("Sto_strategic_opt_hyperbolic_lamb_" + to_string(gLamb) + "_rho_" + to_string(gRho));
	cout << "Successfully completed!" << endl;

	// For Fig.8 in the main text:
	// Uncomment the following lines to run the samples with a increasing sequence of lambda values
	//
	// array<double,6> lamb_list{0.75,0.8,0.85,0.9,0.95,1.0};
	// array<double,5> rho_list{0.5,0.55,0.6,0.65,0.7};
	// for (int ii = 0; ii < 6; ii++){
	// 	for (int jj = 0; jj < 5; jj++){
	// 		if (ii > 0 || jj > 0){
	// 			auto start1 = chrono::high_resolution_clock::now();
	// 			double gLamb = lamb_list[ii];
	// 			double gRho = rho_list[jj];
	// 			cout << "Running the strategically-optimal killer solver (hyperbolic BC; random dilutions): " << endl;
	// 			//Build the main class of our PDE solver
	// 			PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gDefeat_thres, gVictory_thres, gEpsilon, gRks, gGamma, gLamb, gRho);
	// 			//Calling the main PDE solver
	// 			myExample.Main_VPI_Solver();

	// 			// Writing Grid Parameters to file
	// 			myExample.writeDomainToFile("Sto_strategic_opt_hyperbolic_lamb_" + to_string(gLamb) + "_rho_" + to_string(gRho));
	// 			cout << "Successfully completed!" << endl;

	// 			auto end1 = chrono::high_resolution_clock::now();
				
	// 			std::chrono::duration<double> elapsed_seconds = end1 - start1;
	// 			cout << "One set of lamb & rho runtime: " << elapsed_seconds.count() << "s\n";
	// 		}
	// 	}
	// }
	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
