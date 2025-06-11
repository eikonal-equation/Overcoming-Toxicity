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
	
	const string gPolicy_name = "constitutive";
	//const string gPolicy_name = "tactical"; // Uncomment this line to use the tactical policy

	array<double,10> lamb_list{0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0};
	array<double,17> rho_list{0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.60};
	
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();
	for (int i = 0; i < 10; i++){
		double gLamb = lamb_list[i];
		const string str1 = to_string(int(gLamb * 1000 + 1e-12));
		const string gFile_name = "../Tactically_optimal/Sto_tactic_valfn_lamb" + str1 + ".dat";
		for (int j = 0; j < 1; j++){
			double gRho = rho_list[j];
			cout << "Running Probabilistic Performance Metric Solver: " << endl;
			//Build the main class of our PDE solver
			PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol,
				gDefeat_thres, gVictory_thres, gEpsilon, gRks, gGamma, gLamb, gRho, gPolicy_name, gFile_name);

			//Calling the main PDE solver
			myExample.Main_PE_Solver();

			// Writing Grid Parameters to file
			std::string file_label = std::to_string(gRefinement_factor);
			myExample.writeDomainToFile("Perf_hyper");

			cout << "Successfully completed!" << endl;
		}
	}

	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
