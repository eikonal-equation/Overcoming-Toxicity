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
  * which implements the linear time-dependent PDE solver for constitutive killers with regular dilutions. 
  * It carrys out time marching (backwards in time) to solve for the value function 
  * based on a first-order semi-Lagrangian discretization.
  * The file also contains the declarations of the iterative map for the limiting behavior of
  * constitutive killers with infinitely many regular dilutions.
  *
  *
  * Details of all of these functions are found in CSolver.cpp.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "MyDefinitions.h"
#include "CGrid_2D.h"

//----------------------------Libraries------------------------------------------
#include <functional>
#include <boost/numeric/odeint.hpp>

//----------------------------Definitions------------------------------------------
using namespace boost::numeric::odeint;
typedef std::array< double, 2 > state_type;

class PDE_Solver
{
public:
	PDE_Solver() = default;
	PDE_Solver(int aRefinement_factor, double aMaxPopSize, double aMaxFrequency,
		double aTol, double aEpsilon, double aRks, double aGamma, double aTf, double aRho)
	{
		fNx = aRefinement_factor * 100; // Number of grid points for in the x-direction (frequency of the killer) (with baseline 100)
		fNy = aRefinement_factor * 100; // Number of grid points for in the y-direction (total population) (with baseline 100)
		fymax = aMaxPopSize; // Max y (total population; default 1.0)
		fxmax = aMaxFrequency; // Max x (frequency of the killer; default 1.0)
		fx0 = 0.0; // Min x (frequency of the killer; default 0.0)
		fy0 = 0.0; // Min y (total population; default 0.0)
		fDx = (fxmax - fx0) / fNx; // Delta x, grid spacing in x-direction (frequency of the killer)
		fDy = (fymax - fy0) / fNy; // Delta y, grid spacing in y-direction (total population)

		fEpsilon = aEpsilon; // Cost of constitutive toxin production
		fRks = aRks; // Intrinsic cell division rates ratio between the killer and the sensitive
		fGamma = aGamma; // Rescaled killing rate
		fRho = aRho; // Rescaled survival rate
		fTf = aTf; // Rescaled inter-dilution time (E.g., T = 1 the inter-dilution time equal to the inverse growth rate of sensitive cells)
		fM = int(10 * fTf * aRefinement_factor); // Number of time slices (with baseline 10*Tf)

		fMyValfn_current.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of  value function at the current time slice
		fMyPolicy_current.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the optimal policy at the current time slice
		fMyValfn_current_N.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of value-function-based N map  at the current time slice
		fMyValfn_previous.resize(boost::extents[fNx + 1][fNy + 1]);	 //Initialize the storage of  value function at the previous time slice
		fMyValfn_previous_N.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of value-function-based N map at the previous time slice

		fMyfreq_lim_old.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the limiting frequency map at the previous iteration
		fMyfreq_lim_new.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the limiting frequency map at the current iteration
		fMyN_lim_old.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the limiting total population at the previous iteration
		fMyN_lim_new.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the limiting total population at the current iteration

		fMyfreq_map.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the frequency map (recoverd from the 0-th slice of the value function)
		fMyN_map.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the value-function-based N map 

		fIterMap.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of map of competitive exclusion for the limiting behavior of constitutive killers
		fCheck.resize(boost::extents[fNx + 1][fNy + 1]); //Initialize the storage of the check matrix for the above map

		fTau = fTf / double(fM); // Delta t (time step in the time direction)
		fTol = aTol; // Tolerance for the convergence of the iterative solver
		fNum_dilution = -1; // Number of dilutions (initialized to -1, will be updated in the iterative solver)
	}

	// Basic stepper:
	// follows given timestep size "dt"
	typedef runge_kutta4<state_type> rk4;

	// Error stepper, used to create the controlled stepper
	typedef runge_kutta_cash_karp54< state_type > rkck54;

	// Controlled stepper:
	// it's built on an error stepper and allows us to have the output at each
	// internally defined (refined) timestep, via integrate_adaptive call
	typedef controlled_runge_kutta< rkck54 > ctrl_rkck54;

	// Update an ode system based on the current state and control with the above steppers
	void ode_system(const state_type& state, state_type& dstate_dt, double t, int aCtrl);

	// Compute the drift part of the population growth dynamics based on the current state and control
	array<double, 2> f_drift(const double aX, const double aY, const int aCtrl);

	// Initialization of matrices which store value function values at a given time slice 
	// (the policy is fixed at "a = 1" for constitutive killers)
	void InitializeMat(int aTime_slice);

	// Initialization of matrices which store value function values at a given time slice (for N map)
	// (the toxin-production rate is fixed at "a = 1" for constitutive killers)
	void InitializeMat_N(int aTime_slice);

	// This function pre-computes all the coefficients and foot of characteristics of the population growth dynamics
	void precompute_coeff(array_3D& fullstate_loc_mat, array_3D_int& fullstate_ind_mat);

	// This function carrys out the time marching (backwards in time) process of computing the value function
	// for a single point in the grid 
	array<double,2> TimeMarch_single_point(const int aI, const int aJ,
		const array_3D& afullstate_loc_mat, const array_3D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	// This function carrys out the time marching (backwards in time) process of computing the N map
	// for a single point in the grid based on the value function computed in TimeMarch_single_point()
	double TimeMarch_single_point_N(const int aI, const int aJ,
		const array_3D& afullstate_loc_mat, const array_3D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	// This is the main linear time-dependent PDE solver using the semi-Lagrangian discretization
	void Main_Solver();

	// This function computes the infinite limits based on the linear PDE results
	void Inf_limit_Solver();

	// This function computes the infinite limits based on the linear PDE results
	// and outputs the number of iterations till competitive exclusion
	void Inf_limit_Solver_OutputIter();

	// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
	// to compute the new values of the killer frequency and total population size 
	// at a single grid point (aI, aJ) with a grid object (aGrid) using the results from the PREVIOUS ITERATION.
	// It updates 2^n cycles at once for the n-th iteration, where n is the number of iterations.
	array<double,2> Ite_single_point(const int aI, const int aJ,CGrid_2D& aGrid);
	
	// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
	// to compute the new values of the killer frequency and total population size
	// at a single grid point (aI, aJ) with a grid object (aGrid) using the results from the INITIAL MAP.
	// It updates 1 cycle per iteration.
	array<double,2> Ite_single_point_one_time_update(const int aI, const int aJ,CGrid_2D& aGrid);

	// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
	// to perform a single step (iteration) of the iterative map for the limiting behavior of constitutive killers.
	// It computes the new values of the killer frequency and total population size at each grid point (i, j)
	// using the results from the previous iteration.
	double Ite_single_step(CGrid_2D& aGrid);

	// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
	// to perform a single step (iteration) of the iterative map for the limiting behavior of constitutive killers.
	// It computes the new values of the killer frequency and total population size at each grid point (i, j)
	// using the results from the initial map.
	double Ite_single_step_one_time_update(CGrid_2D& aGrid);

	// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
	// to perform a single step (iteration) of the iterative map for the limiting behavior of constitutive killers.
	// It computes the new values of the killer frequency and total population size at each grid point (i, j)
	// using the results from the initial map.
	// In addition, it updates the number of iterations till for a single step (iteration).
	double Ite_single_step_one_time_update_OutputIter(CGrid_2D& aGrid, int aIterNum);

	// This function applies our Accelearated Algorithm (Algo. S3 in the Supportiing Information)
	// to perform the iterative map for the limiting behavior of constitutive killers.
	void Iterative_map(CGrid_2D& aGrid);
	
	// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
	// to perform the iterative map for the limiting behavior of constitutive killers.
	void Iterative_map_one_time_update(CGrid_2D& aGrid);

	// This function applies our Standard Algorithm (Algo. S2 in the Supportiing Information)
	// to perform the iterative map for the limiting behavior of constitutive killers.
	//In addition, it updates the number of iterations till competitive exclusion.
	void Iterative_map_one_time_update_OutputIter(CGrid_2D& aGrid);

	// This function writes the domain parameters to a file
	void writeDomainToFile(const string aFilename);

	// destructor
	~PDE_Solver() { cout << "Destructor called." << "\n"; };

private:
	//static paramters
	double fEpsilon; // Cost of constitutive toxin production
	double fRks; // Intrinsic cell division rates ratio between the killer and the sensitive
	double fGamma; // Rescaled killing rate
	double fTf; // Rescaled inter-dilution time (E.g., T = 1 the inter-dilution time equal to the inverse growth rate of sensitive cells)
	double fRho; // Rescaled survival rate
	double fTol; // Tolerance for the convergence of the iterative solver

	//Discretization and definitions
	int fNx; // number of grid points for x // (frequency of the killer)
	int fNy; //number of grid points for y // (total population)
	double fymax; //Max y (total population)
	double fxmax; // Max x (frequency of the killer)
	double fDx; // Delta x, grid spacing in x-direction (frequency of the killer)
	double fDy; // Delta y, grid spacing in y-direction (total population)
	double fx0; // Min x (frequency of the killer)
	double fy0; // Min y (total population)
	double fTau; // Delta t (time step in the time direction)
	int fM; //number of time slices
	int fNum_dilution; // Number of dilutions

	// Value function and policy matrices
	array_2D fMyValfn_current; // Value function at the current time slice
	array_2D fMyValfn_previous; // Value function at the previous time slice
	array_2D fMyValfn_current_N; // Value-function-based N map at the current time slice
	array_2D fMyValfn_previous_N; // Value-function-based N map at the previous time slice
	array_2D_bool fMyPolicy_current; // Optimal policy at the current time slice

	array_2D fMyfreq_lim_old; // Limiting frequency map at the previous iteration
	array_2D fMyfreq_lim_new; // Limiting frequency map at the current iteration
	array_2D fMyN_lim_old; // Limiting total population at the previous iteration
	array_2D fMyN_lim_new; // Limiting total population at the current iteration

	array_2D fMyN_map; // Value-function-based N map (recoverd from the 0-th slice of the value function)
	array_2D fMyfreq_map; // Frequency map (recoverd from the 0-th slice of the value function)

	array_2D_bool fCheck; // Check matrix to store the convergence status of each grid point
	array_2D_int fIterMap; // Map of competitive exclusion for each grid point
};

#endif // !CSOLVER_H

