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
 * Description: This file contains the implementations of member functions that compute
 * the population dynamics
 * and carry out time marching (backwards in time) to solve for the value function and recover the 
 * corresponding myopic-optimal policy based on a first-order semi-Lagrangian discretization.
 * The file also contains the implementation of the iterative map for the limiting behavior of
 * the recoverd myopic-optimal policy.
 *
 *============================================================================*/

 //----------------------Project specific header files---------------------------
#include "CSolver.h"
#include "WriteToFile.h"

//----------------------------Libraries------------------------------------------
#include <omp.h>

// This function computes the drift part of the population growth dynamics 
// based on the current state (aX, aY) (frequency of the killer: "f", total population: "N") and control (aCtrl; toxin-production rate "a").
// It returns an array of size 2, where the first element is the drift in the x-direction (frequency of the killer)
// and the second element is the drift in the y-direction (total population).
inline array<double, 2> PDE_Solver::f_drift(const double aX, const double aY, const int aCtrl)
{
	array<double, 2> drift;
	drift[0] = aX * (1 - aX) * ((1 - aY) * (fRks * (1 - fEpsilon * double(aCtrl)) - 1) + double(aCtrl) * fGamma * aX * aY);
	drift[1] = aY * (1 - aY) * (1 + aX * (fRks * (1 - fEpsilon * double(aCtrl)) - 1)) - double(aCtrl) * fGamma * pow(aY, 2) * aX * (1 - aX);
	return drift;
}

// This function updates the ode system based on the current state and control.
void PDE_Solver::ode_system(const state_type& state, state_type& dstate_dt, double t, int aCtrl) {
	array<double, 2> drift = f_drift(state[0], state[1], aCtrl);
	dstate_dt[0] = drift[0];    // dx/dt = fdot
	dstate_dt[1] = drift[1];   // dy/dt = Ndot
}

// Initialization of matrices which store value function values, optimal policies at a given time slice
void PDE_Solver::InitializeMat(int aTime_slice) {
	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
		for (int j = 0; j < fNy + 1; j++) {// y loop
			if (aTime_slice == fM - 1) { // terminal condition
				if (j == 0) {//zero population
					fMyValfn_previous[i][j] = 0;
				}
				else {
					fMyValfn_previous[i][j] = xi; //terminal condition for the value function
				}

				if (j == 0 || i == 0) {//zero population or zero killers
					fMyValfn_current[i][j] = 0;
				}
				else if (i == fNx) {
					fMyValfn_current[i][j] = 1; // all killers
				}
				else {
					fMyValfn_current[i][j] = 1e-6; //initialize to a small positive value
				}

			}
			else {// lower time slices

				fMyValfn_previous[i][j] = fMyValfn_current[i][j];

				if (j == 0 || i == 0) {//zero population or zero killers
					fMyValfn_current[i][j] = 0;
				}
				else if (i == fNx) {
					fMyValfn_current[i][j] = 1; // all killers
				}
				else {
					fMyValfn_current[i][j] = 1e-6; //initialize to a small positive value
				}
			}
			fMyPolicy_current[i][j] = 0; //initialize the policy to 0 (no toxin production)
		}
	}
}


// Initialization of matrices which store value function values at a given time slice (for N map)
void PDE_Solver::InitializeMat_N(int aTime_slice) {
	double start_time = double(aTime_slice) * fTau;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
		for (int j = 0; j < fNy + 1; j++) {// y loop
			double yj = fy0 + double(j) * fDy; // y coordinate at the j-th grid point
			if (aTime_slice == fM - 1) { // terminal condition
				if (j == 0) {//zero population at the previous slice
					fMyValfn_previous_N[i][j] = 0;
				}
				else {
					fMyValfn_previous_N[i][j] = yj; // terminal condition for N
				}

				if (j == 0) {//zero population at the current slice
					fMyValfn_current_N[i][j] = 0;
				}
				else if (i == 0 || i == fNx ) {//zero killer or zero sensitive; logistic growth
					std::vector<double> my_times;
					std::vector<double> my_f;
					std::vector<double> my_N;
					array<double, 2> xy{ xi,yj };
					// Lambda function as observer
					auto observer = [&my_f, &my_N, &my_times](const state_type& x, double t) {
						//my_states.push_back(x);
						my_f.push_back(x[0]);
						my_N.push_back(x[1]);
						my_times.push_back(t);
					};

					auto system_with_ctrl = [this](const state_type& state, state_type& dstate_dt, double t) {
						this->ode_system(state, dstate_dt, t, 0); };
					//integrate the ode
					integrate_adaptive(ctrl_rkck54(), system_with_ctrl, xy, start_time, fTf, fTau, observer);
					fMyValfn_current_N[i][j] = my_N.back();
				}
				else {
					fMyValfn_current_N[i][j] = 1e-6; //initialize to a small positive value
				}

			}
			else {// lower time slices
				fMyValfn_previous_N[i][j] = fMyValfn_current_N[i][j]; // update the previous value function
				if (j == 0) {//zero population or zero killers
					fMyValfn_current_N[i][j] = 0;
				}
				else if (i == 0 || i == fNx) {//zero killer or zero sensitive; logistic growth
					std::vector<double> my_times;
					std::vector<double> my_f;
					std::vector<double> my_N;
					array<double, 2> xy{ xi,yj };
					// Lambda function as observer
					auto observer = [&my_f, &my_N, &my_times](const state_type& x, double t) {
						//my_states.push_back(x);
						my_f.push_back(x[0]);
						my_N.push_back(x[1]);
						my_times.push_back(t);
					};

					auto system_with_ctrl = [this](const state_type& state, state_type& dstate_dt, double t) {
						this->ode_system(state, dstate_dt, t, 0); };
					//integrate the ode
					integrate_adaptive(ctrl_rkck54(), system_with_ctrl, xy, start_time, fTf, fTau, observer);
					fMyValfn_current_N[i][j] = my_N.back();
				}
				else {
					fMyValfn_current_N[i][j] = 1e-6; //initialize to a small positive value
				}


			}
		}
	}
}

// This function pre-computes all the coefficients and foot of characteristics of the population growth dynamics
// for each grid point (i, j) and control "a" (0 or 1).
// It fills the 4D fullstate_loc_mat with the future state locations for each (i,j,a) tuple, for each state pair (x,y) in the grid.
// and fullstate_ind_mat with the corresponding indices.
// 4D: i -> x ; j -> y ; a (0,1) -> ctrl ; jj (0,1) -> (x,y)
void PDE_Solver::precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat)
{
	for (int i = 1; i < fNx; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
		for (int j = 1; j < fNy + 1; j++) {// y loop
			double yj = fy0 + double(j) * fDy; // y coordinate at the j-th grid point

			array<double, 2> future_state; // future state after Delta t
			array<double, 2> drift; // drift vector (df, dN) for the current state (fi, Nj) and control a
			for (int a = 0; a < 2; a++) {

				// start with the current state (xi, yj)
				future_state[0] = xi;
				future_state[1] = yj;

				// compute the drift vector based on the current state and control
				drift = f_drift(xi, yj, a);

				for (int ii = 0; ii < 2; ii++) {
					future_state[ii] += fTau * drift[ii]; // update the future state using a 1st-order approximation
				}
				
				// compute the indices for the future state
				int x_indx = find_index(future_state[0], fDx, fx0, fxmax, fNx);
				int y_indx = find_index(future_state[1], fDy, fy0, fymax, fNy);

				array<int, 2> future_indx{ x_indx, y_indx };
				// Fill the fullstate_loc_mat and fullstate_ind_mat with the future state and indices
				for (int jj = 0; jj < 2; jj++) {
					fullstate_loc_mat[i][j][a][jj] = future_state[jj];
					fullstate_ind_mat[i][j][a][jj] = future_indx[jj];
				}
			}
		}
	}

}

// This function carries out the time marching (backwards in time) process of computing the value function
// for a single point in the grid (i, j) based on the pre-computed future states and indices.
// aI, aJ (inputs) : i-th and j-th grid points in the x and y directions, respectively.
// afullstate_loc_mat: 4D array containing the future state locations for each (i,j,a) tuple for each state pair
// afullstate_ind_mat: 4D array containing the corresponding indices for each (i,j,a) tuple for each state pair
// aGrid: CGrid_2D object that provides the bilinear interpolation functionality.
// Output: an array of size 2 containing the value function at the current grid point (i, j) for each control (0 and 1).
inline array<double, 2> PDE_Solver::TimeMarch_single_point(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	array<double, 2> returnVal;
	int counter = 0;
	for (int a = 0; a < 2; a++) {
		//extract loc, indx
		int x_indx = afullstate_ind_mat[aI][aJ][a][0];
		int y_indx = afullstate_ind_mat[aI][aJ][a][1];

		double x_loc = afullstate_loc_mat[aI][aJ][a][0];
		double y_loc = afullstate_loc_mat[aI][aJ][a][1];
		
		// Generate the 4-point stencil for bilinear interpolation
		array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn_previous, x_indx, y_indx);
		// Compute the value function at the current grid point (i, j) for control a using bilinear interpolation
		double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

		returnVal[a] = aVal;
	}
	return returnVal;
}

// This function carries out the time marching (backwards in time) process of computing the N map
// for a single point in the grid (i, j) based on the pre-computed future states and indices,
// and the recovered optimal policy from TimeMarch_single_point().
// aI, aJ (inputs) : i-th and j-th grid points in the x and y directions, respectively.
// afullstate_loc_mat: 4D array containing the future state locations for each (i,j,a) tuple for each state pair//
// afullstate_ind_mat: 4D array containing the corresponding indices for each (i,j,a) tuple for each state pair
// aGrid: CGrid_2D object that provides the bilinear interpolation functionality.
// Output: the value of N at the current grid point (i, j) based on the value function computed in TimeMarch_single_point().
inline double PDE_Solver::TimeMarch_single_point_N(const int aI, const int aJ,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid)
{
	// Extract the current policy (0 or 1) for the grid point (aI, aJ)
	bool policy = fMyPolicy_current[aI][aJ];

	// Extract pre-computed loc and indx for the given optimal policy
	int x_indx = afullstate_ind_mat[aI][aJ][int(policy)][0];
	int y_indx = afullstate_ind_mat[aI][aJ][int(policy)][1];

	double x_loc = afullstate_loc_mat[aI][aJ][int(policy)][0];
	double y_loc = afullstate_loc_mat[aI][aJ][int(policy)][1];

	// Generate the 4-point stencil for bilinear interpolation
	array<double, 4> aVal_sq = aGrid.Stencil_for_Bilinear_Interp(fMyValfn_previous_N, x_indx, y_indx);
	// Compute the value of N at the current grid point (i, j) using bilinear interpolation
	double aVal = aGrid.Bilinear_Interp(aVal_sq, x_indx, y_indx, x_loc, y_loc);

	return aVal;
}

// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
// to compute the new values of the killer frequency and total population size 
// at a single grid point (aI, aJ) with a grid object (aGrid) using the results from the PREVIOUS ITERATION.
// It updates 2^n cycles at once for the n-th iteration, where n is the number of iterations.
array<double,2> PDE_Solver::Ite_single_point(const int aI, const int aJ, CGrid_2D& aGrid)
{
	// Extract (f, N) values at the grid point (aI, aJ) from the previous iteration
	// and compute the new N value based on the rho factor.
	double N = fMyN_lim_old[aI][aJ];
	double f = fMyfreq_lim_old[aI][aJ];
	double newN = fRho * N;
	// Initialize the return values for f and N to be -1
	double aVal_f = -1;
	double aVal_N = -1;

	if (aI == 0 || aI == fNx){// indices for zero killers or all killers
		if (aI == 0){
			aVal_f = 0; // zero killers
		}
		else if (aI == fNx){
			aVal_f = 1; // all killers
		}
		// find the index for newN in the y-direction
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy); 
		// create the interval for linear interpolation in the y-direction USING the results of the previous iteration
		array<double, 2> aVal_interval{fMyN_lim_old[aI][N_indx - 1], fMyN_lim_old[aI][N_indx]};
		// Perform linear interpolation to find the value of N at the newN location
		aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}
	else{
		// For other indices, find the indices for f and newN in the x and y directions, respectively
		int f_indx = find_index(f, fDx, fx0, fxmax, fNx);
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);
		// Generate the 4-point stencil for bilinear interpolation in the 2D spatial grid USING the results of the previous iteration
		array<double, 4> aVal_sq_f = aGrid.Stencil_for_Bilinear_Interp(fMyfreq_lim_old, f_indx, N_indx);
		// Perform bilinear interpolation to find the value of f at the newN location
		aVal_f = aGrid.Bilinear_Interp(aVal_sq_f, f_indx, N_indx, f, newN);

		if (f_indx == fNx)
		{   // Project it back to f = 1 - dx for biological reasons
			array<double, 2> aVal_interval{fMyN_lim_old[fNx - 1][N_indx - 1], fMyN_lim_old[fNx - 1][N_indx]};
			aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		}
		else{
			// Generate the 4-point stencil for bilinear interpolation in the 2D spatial grid for N USING the results of the previous iteration
			array<double, 4> aVal_sq_N = aGrid.Stencil_for_Bilinear_Interp(fMyN_lim_old, f_indx, N_indx);
			// Perform bilinear interpolation to find the value of N at the newN location
			aVal_N = aGrid.Bilinear_Interp(aVal_sq_N, f_indx, N_indx, f, newN);
		}
		// Return the interpolated values of f and N as an array
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}	
}

// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
// to compute the new values of the killer frequency and total population size
// at a single grid point (aI, aJ) with a grid object (aGrid) using the results from the INITIAL MAP.
// It updates 1 cycle per iteration.
array<double,2> PDE_Solver::Ite_single_point_one_time_update(const int aI, const int aJ, CGrid_2D& aGrid)
{
	// Extract (f, N) values at the grid point (aI, aJ) from the previous iteration
	// and compute the new N value based on the rho factor.
	double N = fMyN_lim_old[aI][aJ];
	double f = fMyfreq_lim_old[aI][aJ];
	double newN = fRho * N;
	// Initialize the return values for f and N to be -1
	double aVal_f = -1;
	double aVal_N = -1;

	if (aI == 0 || aI == fNx){ // indices for zero killers or all killers
		if (aI == 0){
			aVal_f = 0; // zero killers
		}
		else if (aI == fNx){
			aVal_f = 1; // all killers
		}
		// find the index for newN in the y-direction
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);
		// create the interval for linear interpolation in the y-direction USING the initial N map
		array<double, 2> aVal_interval{fMyN_map[aI][N_indx - 1], fMyN_map[aI][N_indx]};
		// Perform linear interpolation to find the value of N at the newN location
		aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}
	else{
		// For other indices, find the indices for f and newN in the x and y directions, respectively
		int f_indx = find_index(f, fDx, fx0, fxmax, fNx);
		int N_indx = find_index(newN, fDy, fy0, fymax, fNy);
		// Generate the 4-point stencil for bilinear interpolation in the 2D spatial grid for f USING the initial frequency map
		array<double, 4> aVal_sq_f = aGrid.Stencil_for_Bilinear_Interp(fMyfreq_map, f_indx, N_indx);
		// Perform bilinear interpolation to find the value of f at the newN location
		double aVal_f = aGrid.Bilinear_Interp(aVal_sq_f, f_indx, N_indx, f, newN);

		if (f_indx == fNx)
		{   // project it back to f = 1 - dx for biological reasons (USING the initial N map)
			array<double, 2> aVal_interval{fMyN_map[fNx - 1][N_indx - 1], fMyN_map[fNx - 1][N_indx]};
			aVal_N = aGrid.Linear_Interp(aVal_interval, N_indx, newN);
		}
		else{
			// Generate the 4-point stencil for bilinear interpolation in the 2D spatial grid for N USING the initial N map
			array<double, 4> aVal_sq_N = aGrid.Stencil_for_Bilinear_Interp(fMyN_map, f_indx, N_indx);
			// Perform bilinear interpolation to find the value of N at the newN location
			double aVal_N = aGrid.Bilinear_Interp(aVal_sq_N, f_indx, N_indx, f, newN);
		}
		// Return the interpolated values of f and N as an array
		array<double,2> returnVal{aVal_f, aVal_N};
		return returnVal;
	}	
}

// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
// to perform a single step (iteration) of the iterative map for the limiting behavior of the myopic-optimal policy.
// It computes the new values of the killer frequency and total population size at each grid point (i, j)
// using the results from the previous iteration.
// It only takes the grid object (aGrid) as input,and outputs the maximum error encountered during the iteration.
double PDE_Solver::Ite_single_step(CGrid_2D& aGrid)
{
	double new_err = 0; // Initialize the error to 0
	array<double,2> Interped_val;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			// Call the accelerated one to compute the new values of f and N at the SINGLE grid point (i, j)
			Interped_val = Ite_single_point(i, j,aGrid);
			// Update the new values of f and N in the corresponding matrices
			fMyfreq_lim_new[i][j] = Interped_val[0];
			fMyN_lim_new[i][j] = Interped_val[1];
			// Compute the absolute change in the frequency value and update the error
			double this_change = abs(fMyfreq_lim_new[i][j] - fMyfreq_lim_old[i][j]);
			new_err = max(new_err, this_change);
		}
	}
	return new_err; // Return the maximum error encountered during the iteration
}

// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
// to perform a single step (iteration) of the iterative map for the limiting behavior of the myopic-optimal policy.
// It computes the new values of the killer frequency and total population size at each grid point (i, j)
// using the results from the initial map.
// It only takes the grid object (aGrid) as input,and outputs the maximum error encountered during the iteration.
double PDE_Solver::Ite_single_step_one_time_update(CGrid_2D& aGrid)
{
	double new_err = 0; // Initialize the error to 0
	array<double,2> Interped_val;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			// Call the standard one to compute the new values of f and N at the SINGLE grid point (i, j)
			Interped_val = Ite_single_point_one_time_update(i, j,aGrid);
			// Update the new values of f and N in the corresponding matrices
			fMyfreq_lim_new[i][j] = Interped_val[0];
			fMyN_lim_new[i][j] = Interped_val[1];
			// Compute the absolute change in the frequency value and update the error
			double this_change = abs(fMyfreq_lim_new[i][j] - fMyfreq_lim_old[i][j]);
			// Update the maximum error encountered during the iteration
			new_err = max(new_err, this_change);
		}
	}
	return new_err; // Return the maximum error encountered during the iteration
}

// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
// to perform a single step (iteration) of the iterative map for the limiting behavior of the myopic-optimal policy.
// It computes the new values of the killer frequency and total population size at each grid point (i, j)
// using the results from the initial map.
// It takes the grid object (aGrid) and the number of iterations performed (aIterNum) as input. 
// It outputs the maximum error encountered during the iteration
// and updates the map of competitive exclusion (fIterMap) and check matrix (fCheck) accordingly.
double PDE_Solver::Ite_single_step_one_time_update_OutputIter(CGrid_2D& aGrid, int aIterNum)
{
	double new_err = 0; // Initialize the error to 0
	array<double, 2> Interped_val;
	for (int i = 0; i < fNx + 1; i++) { // x loop
		for (int j = 1; j < fNy + 1; j++) {// y loop
			if (!fCheck[i][j]){ // If not convergent yet, we proceed with the update
				// Call the standard one to compute the new values of f and N at the SINGLE grid point (i, j)
				Interped_val = Ite_single_point_one_time_update(i, j, aGrid);
				// Update the new values of f and N in the corresponding matrices
				fMyfreq_lim_new[i][j] = Interped_val[0];
				fMyN_lim_new[i][j] = Interped_val[1];
			}
			// Compute the absolute change in the frequency value and update the error
			double this_change = abs(fMyfreq_lim_new[i][j] - fMyfreq_lim_old[i][j]);

			// If the change is below the tolerance level and the point is not checked yet
			if (this_change < fTol && !fCheck[i][j]) {
				// Update the map of competitive exclusion (fIterMap) and check matrix (fCheck)
				if (fMyfreq_lim_new[i][j] > 1-10*fTol) { // if the frequency is close to 1
					fIterMap[i][j] = aIterNum; // store the positive iteration number indicating the victory of the killers
					fCheck[i][j] = true; // Convergent mark is checked
				}
				else if (fMyfreq_lim_new[i][j] < 10*fTol) { // if the frequency is close to 0
					fIterMap[i][j] = -aIterNum; // store the negative iteration number indicating the victory of the sensitive
					fCheck[i][j] = true; // Convergent mark is also checked
				}
			}
			new_err = max(new_err, this_change); // Update the maximum error encountered during the iteration
		}
	}
	return new_err;
}


// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
// to perform the iterative map for the limiting behavior of the myopic-optimal policy.
// It iteratively updates the values of the killer frequency and total population size
// until the error is below a specified tolerance level (fTol).
void PDE_Solver::Iterative_map(CGrid_2D& aGrid) {
	double new_err = -1; // Initialize the error to -1
	long MaxIter = 10000; // Maximum number of iterations for the iterative map
	// Loop until the error is below the tolerance level
	for (long iter = 1; iter < MaxIter; ++iter)
	{
		// Call the accelerated one to perform a single step (iteration) of the iterative map
		new_err = Ite_single_step(aGrid); 
		// Update the old values of f and N with the new values
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyfreq_lim_new[i][j];
				fMyN_lim_old[i][j] = fMyN_lim_new[i][j];
			}
		}
		cout << "Iteration #: " << iter << " with err = " << new_err << "\n";\

		if (new_err < fTol) { break; }
	}
	
}

// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
// to perform the iterative map for the limiting behavior of the myopic-optimal policy.
// It iteratively updates the values of the killer frequency and total population size
// until the error is below a specified tolerance level (fTol).
void PDE_Solver::Iterative_map_one_time_update(CGrid_2D& aGrid) {
	double new_err = -1; // Initialize the error to -1
	long MaxIter = 10000; // Maximum number of iterations for the iterative map

	// Loop until the error is below the tolerance level
	for (long iter = 1; iter < MaxIter; ++iter)
	{
		// Call the standard one to perform a single step (iteration) of the iterative map
		new_err = Ite_single_step_one_time_update(aGrid);
		// Update the old values of f and N with the new values
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyfreq_lim_new[i][j];
				fMyN_lim_old[i][j] = fMyN_lim_new[i][j];
			}
		}
		
		cout << "Iteration #: " << iter << " with err = " << new_err << "\n";\

		if (new_err < fTol) { break; }
	}
}

// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
// to perform the iterative map for the limiting behavior of the myopic-optimal policy.
// It iteratively updates the values of the killer frequency and total population size
// until the error is below a specified tolerance level (fTol).
// This version updates the map of competitive exclusion (fIterMap) and check matrix (fCheck) accordingly.
void PDE_Solver::Iterative_map_one_time_update_OutputIter(CGrid_2D& aGrid) {
	double new_err = -1; // Initialize the error to -1
	int MaxIter = 10000; // Maximum number of iterations for the iterative map
	// Loop until the error is below the tolerance level
	for (int iter = 1; iter < MaxIter; ++iter)
	{
		// Call the standard one to perform a single step (iteration) of the iterative map
		// that also outputs the iteration number and error 
		// and updates the map of competitive exclusion (fIterMap) and check matrix (fCheck)
		new_err = Ite_single_step_one_time_update_OutputIter(aGrid, iter);
		// Update the old values of f and N with the new values
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyfreq_lim_new[i][j];
				fMyN_lim_old[i][j] = fMyN_lim_new[i][j]; 
			}
		}
		cout << "Iteration #: " << iter << " with err = " << new_err << "\n"; \
		fNum_dilution = iter; // Update the number of dilutions

		if (new_err < fTol) { break; }
	}
}

// This is the main myopic-optimal PDE solver
void PDE_Solver::Main_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);

	string filename1 = "Regular_Myopic_opt_valuefn.dat";
	string filename2 = "Regular_Myopic_opt_policy.dat";

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	//Main Time Marching Solver
	for (int t = fM - 1; t >= 0; t--) // (backwards in time)
	{	
		// Initialize the value function and policy matrices for the current time slice
		InitializeMat(t);
		// Initial write of the value function and policy matrices to files
		if (t == fM - 1) {
			io::writeToFile2D<double>(filename1, fMyValfn_previous);
			io::writeToFile2D<bool>(filename2, fMyPolicy_current);
		}
		// Parallel loop to compute the value function and policy for each grid point (i, j)
		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				// Call the time marching function for a SINGLE point (i, j) to compute the value function
				array<double, 2> InterpVals;
				InterpVals
					= TimeMarch_single_point(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);
				// Update the value function and policy based on the computed values
				if (InterpVals[0] >= InterpVals[1]) {
					fMyValfn_current[i][j] = InterpVals[0];
					fMyPolicy_current[i][j] = 0;
				}
				else {
					fMyValfn_current[i][j] = InterpVals[1];
					fMyPolicy_current[i][j] = 1;
				}
			}

		}
		// Append the computed value function and policy to the output files
		io::AppendToFile2D<double>(filename1, fMyValfn_current);
		io::AppendToFile2D<bool>(filename2, fMyPolicy_current);
		cout << "Time slice: " << t << endl;
	}
}


// This function computes the infinite limits based on the optimal PDE results
void PDE_Solver::Inf_limit_Solver() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);

	string filename1 = "Regular_Myopic_opt_limit_f.dat";
	string filename2 = "Regular_Myopic_opt_limit_N.dat";

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	//Main Time Marching Solver
	for (int t = fM - 1; t >= 0; t--) // (backwards in time)
	{
		// Initialize the value function and policy matrices for the current time slice
		InitializeMat(t);
		// Initial write of the value functions to files
		if (t == fM - 1) {
			io::writeToFile2D<double>(filename1, fMyValfn_previous);
			io::writeToFile2D<double>(filename2, fMyValfn_previous_N);
		}
		// Parallel loop to compute the value function and policy for each grid point (i, j)
		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				// Call the time marching function for a SINGLE point (i, j) to compute the value function
				array<double, 2> InterpVals;
				InterpVals
					= TimeMarch_single_point(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);

				// Update the value function and policy based on the computed values
				if (InterpVals[0] >= InterpVals[1]) {
					fMyValfn_current[i][j] = InterpVals[0];
					fMyPolicy_current[i][j] = 0;
				}
				else {
					fMyValfn_current[i][j] = InterpVals[1];
					fMyPolicy_current[i][j] = 1;
				}

			}

		}
		// Initialize the N map for the current time slice
		InitializeMat_N(t);
		// Parallel loop to compute the N map for each grid point (i, j)
		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				// Call the time marching function for a SINGLE point (i, j) to compute the N map
				fMyValfn_current_N[i][j]
					= TimeMarch_single_point_N(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);
			}
		}

		// Update the f and N limit matrices for the current time slice
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyValfn_current[i][j];
				fMyfreq_lim_new[i][j] = fMyfreq_lim_old[i][j];
				fMyN_lim_old[i][j] = fMyValfn_current_N[i][j];
				fMyN_lim_new[i][j] = fMyN_lim_old[i][j];

				fMyfreq_map[i][j] = fMyValfn_current[i][j];
				fMyN_map[i][j] = fMyValfn_current_N[i][j];
			}
		}

		// Call the Accelerated iterative map function to compute the limiting behavior of the myopic-optimal policy
		Iterative_map(myGrid);

		// Uncomment the following line to use the Standard iterative map with one-time update
		// Iterative_map_one_time_update(myGrid);

		// Append the computed value function and policy to the output files
		io::AppendToFile2D<double>(filename1, fMyfreq_lim_new);
		io::AppendToFile2D<double>(filename2, fMyN_lim_new);

		cout << "Time slice: " << t << endl;
	}
}


// This function computes the infinite limits based on the optimal PDE results
// and outputs the results for each iteration step.
void PDE_Solver::Inf_limit_Solver_OutputIter() {
	// Initialization
	array_4D myfull_loc_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	array_4D_int myfull_ind_mat(boost::extents[fNx + 1][fNy + 1][2][2]);
	CGrid_2D myGrid(fDx, fDy, fx0, fy0, fNx, fNy);

	auto start = std::chrono::steady_clock::now();
	// precomute all the coefficients
	precompute_coeff(myfull_loc_mat, myfull_ind_mat);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time for computing coefficients: " << elapsed_seconds.count() << "s\n";

	//Main Time Marching Solver
	for (int t = fM - 1; t >= 0; t--) // (backwards in time)
	{
		// Initialize the value function for the current time slice
		InitializeMat(t);

		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				// Call the time marching function for a SINGLE point (i, j) to compute the value function
				array<double, 2> InterpVals;
				InterpVals
					= TimeMarch_single_point(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);

				// Update the value function and policy based on the computed values
				if (InterpVals[0] >= InterpVals[1]) {
					fMyValfn_current[i][j] = InterpVals[0];
					fMyPolicy_current[i][j] = 0;
				}
				else {
					fMyValfn_current[i][j] = InterpVals[1];
					fMyPolicy_current[i][j] = 1;
				}

			}

		}

		// Initialize the N map for the current time slice
		InitializeMat_N(t);
		#pragma omp parallel for schedule(static,1)
		for (int i = 1; i < fNx; i++) { // x loop
			for (int j = 1; j < fNy + 1; j++) {// y loop
				// Call the time marching function for a SINGLE point (i, j) to compute the N map
				fMyValfn_current_N[i][j] 
					= TimeMarch_single_point_N(i, j, myfull_loc_mat, myfull_ind_mat, myGrid);
			}
		}

		// Update the f and N limit matrices for the current time slice
		for (int i = 0; i < fNx + 1; i++) { // x loop
			for (int j = 0; j < fNy + 1; j++) {// y loop
				fMyfreq_lim_old[i][j] = fMyValfn_current[i][j];
				fMyfreq_lim_new[i][j] = fMyfreq_lim_old[i][j];
				fMyN_lim_old[i][j] = fMyValfn_current_N[i][j];
				fMyN_lim_new[i][j] = fMyN_lim_old[i][j];

				fMyfreq_map[i][j] = fMyValfn_current[i][j];
				fMyN_map[i][j] = fMyValfn_current_N[i][j];

				// Initialize the map of competitive exclusion and check matrix for convergence
				fCheck[i][j] = 0;
				fIterMap[i][j] = 0;
			}
		}

		cout << "Time slice: " << t << endl;
	}

	// Call the standard iterative map function to compute the limiting behavior, 
	// which also updates the map of competitive exclusion and check matrix
	Iterative_map_one_time_update_OutputIter(myGrid);

	// Write the results to files
	// The results include the map of competitive exclusion and the limiting frequency.
	// The filenames are constructed based on the rho factor and time slice.
	string aStr = to_string(int(fTf));
	string aStr2 = to_string(int(fRho*1000000 + 1e-12));
	string filename1 = "Regular_myopic_opt_limit_iter_rho" + aStr2 + "_T" + aStr + ".dat";
	io::writeToFile2D<int>(filename1, fIterMap);
	string filename2 = "Regular_myoptic_opt_limit_f_rho" + aStr2 + "_T" + aStr + ".dat";
	io::writeToFile2D<double>(filename2, fMyfreq_lim_new);
}

// This function writes the domain parameters to a file for later reference
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

	params[7] = fTf;
	cout << "Tf = " << params[7] << endl;

	params[8] = fM;
	cout << "fM = " << params[8] << endl;

	params[9] = fTau;
	cout << "Tau = " << params[9] << endl;

	params[10] = fNum_dilution;
	cout << "Num of dilution = " << params[10] << endl;


	// Writing vector of parameters to file
	io::writeToFile1D<double>(aFilename + "_DomainParameters.dat", params);
}
