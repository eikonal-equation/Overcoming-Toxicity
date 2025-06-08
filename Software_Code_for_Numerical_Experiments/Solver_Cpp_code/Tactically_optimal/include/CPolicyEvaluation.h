#pragma once
#ifndef CPOLICYEVALUATION_H
#define CPOLICYEVALUATION_H
/*=============================================================================
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
  * File: CPolicyEvaluation.h
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the declarations of the class "CPolicyEvaluation"
  * which implements the policy evaluation for the intermediate step of 
  * the Value-Policy Iteration (VPI) algorithm that boosts the convergence rate.
  * It contains functions to evaluate the value function for a given policy slice
  * by solving a linear system of equations with a semi-Lagrangian discretization
  * under the condition of a tactically-optimal control problem.
  * It uses the Eigen library for both direct sparse LU solver and iterative approximate solver.
  *
  * Details of all of these functions are found in CPolicyEvaluation.cpp.
  * 
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "MyDefinitions.h"
//----------------------------Libraries------------------------------
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iomanip>

using namespace Eigen;

class CPolicyEvaluation
{
public:
	// constructors
	CPolicyEvaluation() = default;
	CPolicyEvaluation(double aXDim, double aYDim, double aX0, double aY0, double aDx, double aDy,double aProb_not_arrival, double aProb_arrival)
	{
		fNx = aXDim; // number of grid points for x (x: frequency of the killer in this project)
		fNy = aYDim; // number of grid points for y (y: total population in this project)
		fx0 = aX0; // Xmin (default: 0 frequency)
		fy0 = aY0; // Ymin (default: 0 population)

		fDx = aDx; // Delta x, grid spacing in x-direction
		fDy = aDy; // Delta y, grid spacing in y-direction

		fProb_arrival = aProb_arrival; // probability of dilution arriving at the next time step
		fProb_not_arrival = aProb_not_arrival; // probability of dilution not arriving at the next time step

	};

	// This function evaluates the value function for a given policy slice
	// by solving a linear system of equations with a semi-Lagrangian discretization.
	// It uses the Eigen library for a direct sparse LU solver.
	void PolicyEvaluation(const array_2D_bool& aPolicySlice,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, array_2D& aValfn);

	// This function also evaluates the value function for a given policy slice
	// by solving a linear system of equations with a semi-Lagrangian discretization.
	// But this version uses the Eigen library for a iterative approximate solver with a given tolerance.
	void ApproximatePolicyEvaluation(const array_2D_bool& aPolicySlice,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat,
		array_2D& aValfn, double tolerance);

	// destructor
	~CPolicyEvaluation() {};

private:
	//static paramters
	double fProb_not_arrival; // probability of dilution not arriving at the next time step 
	double fProb_arrival; // probability of dilution arriving at the next time step

	//Discretization and definitions
	int fNx; //number of grid points for x
	int fNy; //number of grid points for y
	double fDx; // Delta x, grid spacing in x-direction
	double fDy; // Delta y, grid spacing in y-direction
	double fx0; // Xmin (default: 0)
	double fy0; // Ymin (default: 0)


};
#endif // !CPOLICYEVALUATION_H

