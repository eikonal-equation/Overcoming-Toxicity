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
  * File: CSolver.cpp
  *
  * Author: MingYi Wang
  * 
  * Description: This file contains the acutal implementation of the class "CPDE_Solver"
  * which constructs the PDE solver for the strategically-optimal control problem.
  * It carrys out a Gauss–Seidel Value-Policy-Iteration (GS-VPI) algorithm to 
  * solve for the value function and recover the corresponding strategically-optimal policy 
  * based on a first-order semi-Lagrangian discretization.
  * Here, we specifically handle the case of a vertical boundary condition.
  *
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "CSolver.h"
#include "WriteToFile.h"

// This function computes the drift part of the population growth dynamics 
// based on the current state (aX, aY) (frequency of the killer: "f", total population: "N") and control (aCtrl; toxin-production rate "a").
// It returns an array of size 2, where the first element is the drift in the x-direction (frequency of the killer)
// and the second element is the drift in the y-direction (total population).
array<double, 2> PDE_Solver::f_drift(const double aX, const double aY, const int aCtrl)
{
	array<double, 2> drift;
	drift[0] = aX * (1 - aX) * ((1 - aY) * (fRks * (1 - fEpsilon * double(aCtrl)) - 1) + double(aCtrl) * fGamma * aX * aY);
	drift[1] = aY * (1 - aY) * (1 + aX * (fRks * (1 - fEpsilon * double(aCtrl)) - 1)) - double(aCtrl) * fGamma * pow(aY, 2) * aX * (1 - aX);
	return drift;
}


// This function initializes the value function and policy matrices.
void PDE_Solver::InitializeMat() {
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 0; j < fNy + 1; j++) {// y loop;
			if (j == 0) {//zero population
				fMyValfn[i][j] = 0;
			}
			else if (i < fDefeat_indx) {
				fMyValfn[i][j] = 0; //Defeat zone: zero killers eventually
			}
			else if (i > fVic_indx) {
				fMyValfn[i][j] = 1; //Victory zone: all killers eventually
			}
			else {
				fMyValfn[i][j] = 1e-6; // small value for the rest of the grid points
			}
			fMyPolicy[i][j] = 0; // initialize all policies to 0 (no toxin production)
		}
	}
}

// This function pre-computes all the coefficients and foot of characteristics of the population growth dynamics
// for each grid point (i, j) and control "a" (0 or 1).
// It fills the 4D fullstate_loc_mat with the future state locations for each (i,j,a) tuple, for each state pair (x,y) in the grid.
// and fullstate_ind_mat with the corresponding indices.
// 4D: i -> x ; j -> y ; a (0,1) -> ctrl ; jj (0,1,2) -> (x,y,rho*y)
void PDE_Solver::precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat)
{
	for (int i = fDefeat_indx; i < fVic_indx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
		for (int j = 1; j < fNy + 1; j++) {// y loop 
			double yj = fy0 + double(j) * fDy; // y coordinate at the j-th grid point
			array<double, 3> future_state;
			array<double, 2> drift;
			for (int a = 0; a < 2; a++) { // Loop over the two controls (0 and 1)
				// Start with the current state (xi, yj)
				future_state[0] = xi;
				future_state[1] = yj;
				// Compute the drift based on the current state and control
				drift = f_drift(xi, yj, a);

				for (int ii = 0; ii < 2; ii++) {
					future_state[ii] += fTau * drift[ii]; // update the future state using a 1st-order approximation
				}
				future_state[2] = fRho * future_state[1]; // y_rho = rho*y

				// compute the indices for the future state
				int x_indx = find_index(future_state[0], fDx, fx0, fxmax, fNx);
				int y_indx = find_index(future_state[1], fDy, fy0, fymax, fNy);
				int yrho_indx = find_index(future_state[2], fDy, fy0, fymax, fNy);

				array<int, 3> future_indx{ x_indx, y_indx, yrho_indx };
				// Fill the fullstate_loc_mat and fullstate_ind_mat with the future state and indices
				for (int jj = 0; jj < 3; jj++) {
					fullstate_loc_mat[i][j][a][jj] = future_state[jj];
					fullstate_ind_mat[i][j][a][jj] = future_indx[jj];
				}
			}
		}
	}
}

// This function carries out the Gauss–Seidel (GS) Value iterations of computing the value function
// for a single point in the grid (i, j) based on the pre-computed future states and indices.
// aI, aJ (inputs) : i-th and j-th grid points in the x and y directions, respectively.
// afullstate_loc_mat: 4D array containing the future state locations for each (i,j,a) tuple for each state pair
// afullstate_ind_mat: 4D array containing the corresponding indices for each (i,j,a) tuple for each state pair
// aGrid: CGrid_2D object that provides the bilinear interpolation functionality.
// Output: an array of size 2 containing the value function at the current grid point (i, j) for each control (0 and 1).
array<double,4> PDE_Solver::ValueIte_single_point(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	array<double, 4> returnVal;
	int counter = 0; // Initialize the counter for the return array
	for (int a = 0; a < 2; a++) { // Loop over the two controls (0 and 1)
		// Get the indices and locations of the future state for the current control
		int x_indx = afullstate_ind_mat[aI][aJ][a][0];
		int y_indx = afullstate_ind_mat[aI][aJ][a][1];
		int yrho_indx = afullstate_ind_mat[aI][aJ][a][2];
		double x_loc = afullstate_loc_mat[aI][aJ][a][0];
		double y_loc = afullstate_loc_mat[aI][aJ][a][1];
		double yrho_loc = afullstate_loc_mat[aI][aJ][a][2];

		// Generate the 4-point stencil for bilinear interpolation
		array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn, x_indx, y_indx);
		// Compute the value function at the current grid point (i, j) for control a using bilinear interpolation
		double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

		// Generate the 4-point stencil for bilinear interpolation for the diluted one (rho*y)
		array<double, 4> aVal_sq_rho = aGrid.Stencil_for_Bilinear_Interp(fMyValfn, x_indx, yrho_indx);
		// Compute the value function at the current grid point (i, j) for control a using bilinear interpolation
		// for the diluted one (rho*y)
		double aVal_rho = aGrid.Bilinear_Interp(aVal_sq_rho, x_indx, yrho_indx, x_loc, yrho_loc);

		returnVal[counter] = aVal;
		returnVal[counter + 1] = aVal_rho;
		counter = counter + 2;
	}
	return returnVal;
}

// This function carrys out the Gauss–Seidel (GS) Value iterations of computing the value function
// for an entire iteration (i.e., looping over all girdpoints).
// afullstate_loc_mat: 4D array containing the future state locations for each (i,j,a) tuple for each state pair
// afullstate_ind_mat: 4D array containing the corresponding indices for each (i,j,a) tuple for each state pair
// aGrid: CGrid_2D object that provides the bilinear interpolation functionality.
// Output: the maximum error of the value function update across all grid points in this iteration.
double PDE_Solver::ValueIte_single_step(const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	double new_err = 0; // Initialize the maximum error for this iteration
	array<double, 4> InterpVals;
	for (int i = fDefeat_indx; i < fVic_indx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			array<double, 2> w_possible;
			int counter = 0; // Initialize the counter for the return array
			
			// Call ValueIte_single_point() to compute the value function at the current grid point (i, j)
			InterpVals
				= ValueIte_single_point(i, j, afullstate_loc_mat, afullstate_ind_mat, aGrid);

			// Compute the expected value function based on the probabilities of arrival and not arrival
			for (int ii = 0; ii < 2; ii++) {
				w_possible[ii] = fProb_not_arrival * InterpVals[counter] + fProb_arrival * InterpVals[counter+1];
				counter = counter + 2;
			}

			// Determine the maximum value and the corresponding policy (0 or 1)
			double max_val = 0;
			bool mypolicy = 0;
			if (w_possible[0] >= w_possible[1]) {
				max_val = w_possible[0];
				mypolicy = 0;
			}
			else {
				max_val = w_possible[1];
				mypolicy = 1;
			}

			// GS-VI: update the value function and policy if the maximum value is greater than the current value
			if (max_val > fMyValfn[i][j]) {

				double this_change = max_val - fMyValfn[i][j];

				fMyValfn[i][j] = max_val;

				new_err = max(new_err, this_change);

				fMyPolicy[i][j] = mypolicy;
			}
		}
	}
	return new_err;
}

// This is the main tactically-optimal control problem solver 
// with Gauss–Seidel Value-Policy-Iteration (GS-VPI) algorithm
void PDE_Solver::Main_VPI_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][3]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);
	CPolicyEvaluation myPEclass(fNx, fNy, fx0, fy0, fDx, fDy, fDefeat_thres, fVictory_thres, fDefeat_indx, fVic_indx, fProb_not_arrival, fProb_arrival);

	InitializeMat(); // Initialize the value function and policy matrices

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	double old_err = -1; // Initialize the old error to -1
	double new_err = -1; // Initialize the new error to -1
	int PE_counter = 0; // Counter for the number of policy evaluation steps
	
	for (long iter = 0; iter < fMaxIter; ++iter)
	{
		// Perform a single step of the Gauss–Seidel Value iteration and ouput the error
		new_err = ValueIte_single_step(myfull_loc_mat, myfull_ind_mat, myGrid);
		if (iter == 0) { old_err = new_err; } // Set the old error to the new error on the first iteration
		if (new_err < 0.1 && new_err < fDelta * old_err) {// If the error is small enough, perform policy evaluation
			old_err = new_err;

			// Uncomment the following line to use the direct solver for policy evaluation
			//myPEclass.PolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn);

			// Perform approximate policy evaluation
			myPEclass.ApproximatePolicyEvaluation(fMyPolicy, myfull_loc_mat, myfull_ind_mat, fMyValfn,1e-6);
			PE_counter++; // Increment the policy evaluation counter
		}
		cout << "Iteration #: " << iter + 1 << " with err = " << new_err << "\n";
		cout << "# of PE: " << PE_counter << "\n";
		if (new_err < fTol) { break; }
	}

	// Write the optimal value function and policy to files
	string aStr = to_string(int(fLamb*1000+1e-12));
	string aStr2 = to_string(int(fRho*1000 + 1e-12));
	string filename1 = "Sto_strategic_vertical_valuefn_rho" + aStr2 + "_lamb" + aStr + ".dat";
	string filename2 = "Sto_strategic_vertical_policy_rho" + aStr2 + "_lamb" + aStr + ".dat";
	io::writeToFile2D<double>(filename1, fMyValfn);
	io::writeToFile2D<bool>(filename2, fMyPolicy);
}

// This function writes the domain parameters to a file
void PDE_Solver::writeDomainToFile(const string aFilename) {

	array_1D params(boost::extents[11]); // Creating boostarray "row vector"

	params[0] = fNx;
	cout << "fNx = " << params[0] << endl;

	params[1] = fNy;
	cout << "fNy = " << params[1] << endl;

	params[2] = fDx;
	cout << "Delta x = " << params[2] << endl;

	params[3] = fDy;
	cout << "Delta y = " << params[3] << endl;

	params[4] = fEpsilon;
	cout << "Eps = " << params[4] << endl;

	params[5] = fRks;
	cout << "Rks = " << params[5] << endl;

	params[6] = fGamma;
	cout << "Gamma = " << params[6] << endl;

	params[7] = fLamb;
	cout << "Lamb = " << params[7] << endl;

	params[8] = fRho;
	cout << "Rho = " << params[8] << endl;

	params[9] = fVictory_thres;
	cout << "Victory thres = " << params[9] << endl;

	params[10] = fDefeat_thres;
	cout << "Defeat thres = " << params[10] << endl;

	// Writing vector of parameters to file
	io::writeToFile1D<double>(aFilename + "_DomainParameters.dat", params);
}
