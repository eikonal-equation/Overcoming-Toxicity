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
	//=============================================Initilization=====================================================
	const int gRefinement_factor = 16; //Global refinement factor with respect to the 3D (x,y,t) grid (x&y for spatial; t for time)
	const double gMax_population = 1; //Maximum nomralized population size (max y value for the spatial grid)
	const double gMax_frequency = 1; //Maximum nomralized frequency of the killer (max x value for the spatial grid)
	const double gEpsilon = 0.2; //The cost of constitutive toxin production
	const double gRks = 0.85; //Intrinsic cell division rates ratio between the killer and the sensitive
	const double gGamma = 1.0; //Rescaled killing rate
	const double gTf = 1; //Rescaled inter-dilution time (E.g., T = 1 the inter-dilution time equal to the inverse growth rate of sensitive cells)
	const double gRho = 0.65; //Rescaled survival rate
	const double gTol = 1e-6; //Tolerance for the convergence of the iterative solver
	
	//=============================================Start of Our PDE Solver=====================================================
	auto start = chrono::high_resolution_clock::now();

	//Build the main class of our PDE solver
	PDE_Solver myExample(gRefinement_factor, gMax_population, gMax_frequency, gTol, gEpsilon, gRks, gGamma, gTf, gRho);

	// Uncomment the following two lines to run the main linear PDE solver
	// cout << "Running the constitutive killer solver (regular dilutions): " << endl;
	// myExample.Main_Solver();

	// Uncomment the following two lines to run the PDE solver + the Infinite limit solver
	//cout << "Running the constitutive killer solver followed by the infinite limit solver (regular dilutions): " << endl;
	//myExample.Inf_limit_Solver();

	// Running the main linear PDE solver with the infinite limit solver 
	// (that also outputs the map of exclusion) for the limiting behavior
	cout << "Running the constitutive killer solver followed by the infinite limit solver that also outputs the map of competitive exclusion (regular dilutions): " << endl;
	myExample.Inf_limit_Solver_OutputIter();

	// Writing Grid Parameters to file
	myExample.writeDomainToFile("Regular_dilutions_constitutive");

	auto end = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Total runtime: " << elapsed_seconds.count() << "s\n";
	cout << "Successfully completed!" << endl;
	//=============================================END of Our PDE Solver=====================================================
	return 0;
}
