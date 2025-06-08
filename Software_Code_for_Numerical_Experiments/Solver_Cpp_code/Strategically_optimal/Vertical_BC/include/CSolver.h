#pragma once
#ifndef CSOLVER_H
#define CSOLVER_H
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
  * File: CSolver.h
  *
  * Author: MingYi Wang
  * 
  * Description: This file contains the declarations of the class "CPDE_Solver"
  * which constructs the PDE solver for the strategically-optimal control problem 
  * with vertical boundary conditions.
  *	It carrys out a Gauss–Seidel Value-Policy-Iteration (GS-VPI) algorithm to 
  * solve for the value function and recover the corresponding strategically-optimal policy 
  * based on a first-order semi-Lagrangian discretization.
  * Here, we specifically handles the case of a vertical boundary condition.
  *
  *
  * Details of all of these functions are found in CSolver.cpp.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "MyDefinitions.h"
#include "CGrid_2D.h"
#include "CPolicyEvaluation.h"


class PDE_Solver
{
public:
	// constructors
	PDE_Solver() = default;
	PDE_Solver(int aRefinement_factor, double aMaxPopSize, double aMaxFrequency,
		double aTol, double aDefeat_thres, double aVictory_thres,
		double aEpsilon, double aRks, double aGamma, double aLamb, double aRho)
	{
		fNx = aRefinement_factor * 100; // Number of grid points for in the x-direction (frequency of the killer) (with baseline 100)
		fNy = aRefinement_factor * 100; // Number of grid points for in the y-direction (total population) (with baseline 100)
		fymax = aMaxPopSize; // Max y (total population; default 1.0)
		fxmax = aMaxFrequency; // Max x (frequency of the killer; default 1.0)
		fx0 = 0.0; // Min x (frequency of the killer; default 0.0)
		fy0 = 0.0; // Min y (total population; default 0.0)
		fDx = (fxmax - fx0) / fNx; // Delta x, grid spacing in x-direction
		fDy = (fymax - fy0) / fNy; // Delta y, grid spacing in y-direction
		
		fDefeat_thres = aDefeat_thres; // Threshold for the defeat (default 0.1)
		fVictory_thres = aVictory_thres; // Threshold for the victory (default 0.9)
		fDefeat_indx = int(fDefeat_thres / fDx); // Index for the defeat threshold
		fVic_indx = int(fVictory_thres / fDx); // Index for the victory threshold

		fEpsilon = aEpsilon; // Cost of constitutive toxin production
		fRks = aRks; // Intrinsic cell division rates ratio between the killer and the sensitive
		fGamma = aGamma; // Rescaled killing rate
		fLamb = aLamb; // Arrival rate of the toxin at the next time step
		fRho = aRho; // Rescaled survival rate
		fDelta = 0.7; // Percentage of reduction of error to enter the policy evaluation step

		fMyValfn.resize(boost::extents[fNx + 1][fNy + 1]); // Storage of the value function for all (x,y) grid points
		fMyPolicy.resize(boost::extents[fNx + 1][fNy + 1]); // Storage of the optimal policy for all (x,y) grid points

		fTau = sqrt(fDx); // Delta t (time step size)
		fMaxIter = 100000; // Maximum number of iterations for the iterative solver
		fTol = aTol; // Tolerance for the convergence of the iterative solver

		fProb_not_arrival = exp(-fLamb * fTau); // Probability of dilution not arriving at the next time step
		fProb_arrival = 1.0 - fProb_not_arrival; // Probability of dilution arriving at the next time step
	}

	// Compute the drift part of the population growth dynamics based on the current state and control
	array<double, 2> f_drift(const double aX, const double aY, const int aCtrl);

	// Initialization of matrices which store the value function and the optimal policy
	void InitializeMat();

	// This function pre-computes all the coefficients and foot of characteristics of the population growth dynamics
	void precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat);

	// This function computes the Gauss–Seidel (GS) Value iterations for a single gridpt
	array<double,4> ValueIte_single_point(const int aI, const int aJ,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	// This function computes the Gauss–Seidel (GS) Value iterations for a single step (iteration)
	// I.e., looping over all grid points (i,j) once
	double ValueIte_single_step(const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	// This is the main strategically-optimal control problem solver 
	// with Gauss–Seidel Value-Policy-Iteration (GS-VPI) algorithm
	void Main_VPI_Solver();
	
	// This function writes the value function and the optimal policy to a file
	void writeDomainToFile(const string aFilename);
	
	// destructor
	~PDE_Solver() { cout << "Destructor called." << "\n"; };

private:
	//static paramters
	double fEpsilon; // Cost of constitutive toxin production
	double fRks; // Intrinsic cell division rates ratio between the killer and the sensitive
	double fGamma; // Rescaled killing rate
	double fRho; // Rescaled survival rate
	double fDelta; // Percentage of reduction of error to enter the policy evaluation step
	double fProb_not_arrival; // Probability of dilution not arriving at the next time step
	double fProb_arrival; // Probability of dilution arriving at the next time step
	double fDefeat_thres; // Threshold for the defeat
	double fVictory_thres; // Threshold for the victory
	int fDefeat_indx; // Index for the defeat threshold
	int fVic_indx; // Index for the victory threshold
	

	//Discretization and definitions
	int fNx; //number of grid points for x
	int fNy; //number of grid points for y
	double fymax; //Max y
	double fxmax; //Max x
	double fDx; // Delta x, grid spacing in x-direction
	double fDy; // Delta y, grid spacing in y-direction
	double fx0; // Min x (default: 0)
	double fy0; // Min y (default: 0)
	double fTau; // Delta t (time step size)
	long fMaxIter; // Maximum number of iterations for the iterative solver
	double fTol; // Tolerance for the convergence of the iterative solver
	array_2D fMyValfn; // Storage of the value function for all (x,y) grid points
	array_2D_bool fMyPolicy; // Storage of the optimal policy for all (x,y) grid points
};

#endif // !CSOLVER_H

