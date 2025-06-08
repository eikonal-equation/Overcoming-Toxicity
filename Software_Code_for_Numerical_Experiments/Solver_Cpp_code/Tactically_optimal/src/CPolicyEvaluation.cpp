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
  * File: CPolicyEvaluation.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the actual implementation of the member functions of
  * class "CPolicyEvaluation". 
  * It contains functions executing the policy evaluation for the intermediate step of 
  * the Value-Policy Iteration (VPI) algorithm that boosts the convergence rate.
  * These functions evaluate the value function for a given policy slice
  * by solving a linear system of equations with a semi-Lagrangian discretization
  * under the condition of a tactically-optimal control problem.
  * It uses the Eigen library for both direct sparse LU solver and iterative approximate solver.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "CPolicyEvaluation.h"
#include "WriteToFile.h"

// This function evaluates the value function for a given policy slice
// by solving a linear system of equations with a semi-Lagrangian discretization,
// using the Eigen library for a direct sparse LU solver.
//
// aPolicySlice (input): 2D array of booleans representing the policy slice
// afullstate_loc_mat (input): 4D array containing the locations of the full state
// afullstate_ind_mat (input): 4D array containing the indices of the full state
// aValfn (input): 2D array to store the evaluated value function
// The function constructs a sparse matrix A and a vector b, then solves the system Av = b
// and updates the value function aValfn with the solution v.
void CPolicyEvaluation::PolicyEvaluation(const array_2D_bool& aPolicySlice,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, array_2D& aValfn) 
{
	int mSize = (fNx + 1) * (fNy + 1); // total number of grid points; size of the linear system

	VectorXd b(mSize); // right-hand side vector for the linear system
	VectorXd v(mSize); // solution vector for the linear system
	SparseMatrix<double> A(mSize, mSize); // sparse matrix storing the coefficients of the linear system
	vector<Triplet<double>> entries; // vector to store non-zero entries of the sparse matrix

	int aXindex; // right index for the x-coordinate in the full state
	int aYindex; // right index for the y-coordinate in the full state
	double aXloc; // location for the x-coordinate in the full state
	double aYloc; // location for the y-coordinate in the full state
	double x_pos_tr; // x position corresponding to aXindex
	double y_pos_tr; // y position corresponding to aYindex


	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
 		for (int j = 0; j < fNy + 1; j++) { // y loop
			// Add diagonal entry for the current grid point (i, j)
			entries.push_back(Triplet<double>(i * (fNx + 1) + j, i * (fNx + 1) + j, 1.0));
			
			if (j > 0 && i > 0 && i < fNx && aValfn[i][j] != 1 && aValfn[i][j] != 0)
			{
				int current_policy = int(aPolicySlice[i][j]); // current policy for the grid point (i, j)

				b(i * (fNx + 1) + j) = fProb_arrival*xi; // RHS: probability of arrival times the x-coordinate

				// Get the indices and locations for the current policy
				aXindex = afullstate_ind_mat[i][j][current_policy][0];
				aYindex = afullstate_ind_mat[i][j][current_policy] [1];
				aXloc = afullstate_loc_mat[i][j][current_policy][0];
				aYloc = afullstate_loc_mat[i][j][current_policy][1];

				// Calculate the x and y positions corresponding to the indices
				x_pos_tr = aXindex * fDx + fx0;
				y_pos_tr = aYindex * fDy + fy0;

				// Add entries to the sparse matrix A for the semi-Lagrangian discretization
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (x_pos_tr - aXloc) / fDx)); //bottom left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //bottom right
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (x_pos_tr - aXloc) / fDx)); //top left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //top right
			}
			else {
				b(i * (fNx + 1) + j) = aValfn[i][j]; // if the policy is not applicable, set the RHS to the current value function
			}
		}
	}

	A.setFromTriplets(entries.begin(), entries.end()); // Construct the sparse matrix A from the triplets

	//-------------- solving the linear system Av = b by LU factorization ----------------------
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	solver.compute(A);

	if (solver.info() != Success) {
		cout << "LU Decomposition Failed" << endl;
		myAssert(solver.info() == Success);
	}
	v = solver.solve(b);

	//------------ update the top slice of our value function-----------------
	for (int i = 0; i < fNx + 1; i++) {
		for (int j = 0; j < fNy + 1; j++) {
			aValfn[i][j] = v(i * (fNx + 1) + j);
		}
	}
}

// This function evaluates the value function for a given policy slice
// by solving a linear system of equations with a semi-Lagrangian discretization,
// using the Eigen library for an iterative approximate solver with a given tolerance.
//
// aPolicySlice (input): 2D array of booleans representing the policy slice
// afullstate_loc_mat (input): 4D array containing the locations of the full state
// afullstate_ind_mat (input): 4D array containing the indices of the full state
// aValfn (input): 2D array to store the evaluated value function
// tolerance (input): the tolerance for the iterative solver
// The function constructs a sparse matrix A and a vector b, then solves the system Av = b
// and updates the value function aValfn with the solution v.
void CPolicyEvaluation::ApproximatePolicyEvaluation(const array_2D_bool& aPolicySlice,
	const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat,
	array_2D& aValfn, double tolerance) 
{
	int mSize = (fNx+ 1) * (fNy + 1); // total number of grid points; size of the linear system

	VectorXd b(mSize); // right-hand side vector for the linear system
	VectorXd v(mSize); // solution vector for the linear system
	VectorXd InitialGuess(mSize); // initial guess for the solution vector
	SparseMatrix<double> A(mSize, mSize); // sparse matrix storing the coefficients of the linear system
	vector<Triplet<double>> entries; // vector to store non-zero entries of the sparse matrix

	int aXindex; // right index for the x-coordinate in the full state
	int aYindex; //	right index for the y-coordinate in the full state
	double aXloc; // location for the x-coordinate in the full state
	double aYloc; // location for the y-coordinate in the full state
	double x_pos_tr; // x position corresponding to aXindex
	double y_pos_tr; // y position corresponding to aYindex

	for (int i = 0; i < fNx + 1; i++) { // x loop
		double xi = fx0 + double(i) * fDx; // x coordinate at the i-th grid point
		for (int j = 0; j < fNy + 1; j++) { // y loop
			// Add diagonal entry for the current grid point (i, j)
			entries.push_back(Triplet<double>(i * (fNx + 1) + j, i * (fNx + 1) + j, 1.0));
			// Initialize the initial guess for the value function to be its current value
			InitialGuess(i * (fNx + 1) + j) = aValfn[i][j];

			if (j > 0 && i > 0 && i < fNx && aValfn[i][j] != 1 && aValfn[i][j] != 0)
			{
				int current_policy = int(aPolicySlice[i][j]); // current policy for the grid point (i, j)
				b(i * (fNx + 1) + j) = fProb_arrival*xi; // RHS: probability of arrival times the x-coordinate

				// Get the indices and locations for the current policy
				aXindex = afullstate_ind_mat[i][j][current_policy][0];
				aYindex = afullstate_ind_mat[i][j][current_policy][1];
				aXloc = afullstate_loc_mat[i][j][current_policy][0];
				aYloc = afullstate_loc_mat[i][j][current_policy][1];

				// Calculate the x and y positions corresponding to the indices
				x_pos_tr = aXindex * fDx + fx0;
				y_pos_tr = aYindex * fDy + fy0;

				// Add entries to the sparse matrix A for the semi-Lagrangian discretization
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex - 1, -fProb_not_arrival*(y_pos_tr - aYloc) / fDy * (x_pos_tr - aXloc) / fDx)); //bottom left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex - 1, -fProb_not_arrival * (y_pos_tr - aYloc) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //bottom right
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex - 1) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (x_pos_tr - aXloc) / fDx)); //top left
				entries.push_back(Triplet<double>(i * (fNx + 1) + j, (aXindex) * (fNx + 1) + aYindex, -fProb_not_arrival * (aYloc - (y_pos_tr - fDy)) / fDy * (aXloc - (x_pos_tr - fDx)) / fDx)); //top right
			}
			else {
				b(i * (fNx + 1) + j) = aValfn[i][j]; // if the policy is not applicable, set the RHS to the current value function
			}
		}
	}

	A.setFromTriplets(entries.begin(), entries.end()); // Construct the sparse matrix A from the triplets

	//-------------- solving the linear system Av = b by iterative approximate solver ----------------------
	BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double> > approxSolver;
	approxSolver.compute(A);
	approxSolver.setTolerance(tolerance);
	v = approxSolver.solveWithGuess(b, InitialGuess);
	cout << "Approx Solver Iterations:" << approxSolver.iterations() << "\n";

	//------------ update the top slice of our value function-----------------
	for (int i = 0; i < fNx + 1; i++) {
		for (int j = 0; j < fNy + 1; j++) {
			aValfn[i][j] = v(i * (fNx + 1) + j);
		}
	}
}
