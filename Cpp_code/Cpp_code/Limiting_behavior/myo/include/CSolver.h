#pragma once
#ifndef CSOLVER_H
#define CSOLVER_H
#define _CRT_SECURE_NO_WARNINGS
#include "MyDefinitions.h"
#include "CGrid_2D.h"
#include <functional>
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
typedef std::array< double, 2 > state_type;

class PDE_Solver
{
public:
	PDE_Solver() = default;
	PDE_Solver(int aRefinement_factor, double aMaxPopSize, double aMaxFrequency,
		double aTol, double aEpsilon, double aRks, double aGamma, double aTf, double aRho)
	{
		fNx = aRefinement_factor * 100;
		fNy = aRefinement_factor * 100;
		fymax = aMaxPopSize;
		fxmax = aMaxFrequency;
		fx0 = 0.0;
		fy0 = 0.0;
		fDx = (fxmax - fx0) / fNx;
		fDy = (fymax - fy0) / fNy;

		fEpsilon = aEpsilon;
		fRks = aRks;
		fGamma = aGamma;
		fRho = aRho;
		fTf = aTf;
		fM = int(10 * fTf * aRefinement_factor);

		fMyValfn_current.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyPolicy_current.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyValfn_current_N.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyValfn_previous.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyValfn_previous_N.resize(boost::extents[fNx + 1][fNy + 1]);

		fMyfreq_lim_old.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyfreq_lim_new.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyN_lim_old.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyN_lim_new.resize(boost::extents[fNx + 1][fNy + 1]);

		fMyN_map.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyfreq_map.resize(boost::extents[fNx + 1][fNy + 1]);


		fIterMap.resize(boost::extents[fNx + 1][fNy + 1]);
		fCheck.resize(boost::extents[fNx + 1][fNy + 1]);


		fTau = fTf / double(fM);
		fTol = aTol;
		fNum_dilution = -1;
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

	// compute the ode system
	void ode_system(const state_type& state, state_type& dstate_dt, double t, int aCtrl);

	// drift function f
	array<double, 2> f_drift(const double aX, const double aY, const int aCtrl);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	void InitializeMat(int aTime_slice);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	void InitializeMat_N(int aTime_slice);

	// This function pre-compute all the coefficients and foot of characteristics
	void precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat);

	//This function computes the GS Value iterations for a single gridpt
	array<double,2> TimeMarch_single_point(const int aI, const int aJ,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	double TimeMarch_single_point_N(const int aI, const int aJ,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	// This is the main solver with GS VI
	void Main_Solver();

	void Inf_limit_Solver();

	void Inf_limit_Solver_OutputIter();

	array<double, 2> Ite_single_point(const int aI, const int aJ, CGrid_2D& aGrid);

	array<double, 2> Ite_single_point_one_time_update(const int aI, const int aJ, CGrid_2D& aGrid);


	double Ite_single_step(CGrid_2D& aGrid);

	double Ite_single_step_one_time_update(CGrid_2D& aGrid);

	double Ite_single_step_one_time_update_OutputIter(CGrid_2D& aGrid, int aIterNum);

	void Iterative_map(CGrid_2D& aGrid);

	void Iterative_map_one_time_update(CGrid_2D& aGrid);

	void Iterative_map_one_time_update_OutputIter(CGrid_2D& aGrid);

	void writeDomainToFile(const string aFilename);

	~PDE_Solver() { cout << "Destructor called." << "\n"; };

private:
	//static paramters
	double fEpsilon;
	double fRks;
	double fGamma;
	double fTf;
	double fRho;
	double fTol;

	//Discretization and definitions
	int fNx; //number of grid points for x
	int fNy; //number of grid points for y
	double fymax; //Max y
	double fxmax; //Max x
	double fDx;
	double fDy;
	double fx0;
	double fy0;
	double fTau;
	int fM; //number of time slices
	array_2D fMyValfn_current;
	array_2D fMyValfn_previous;
	array_2D fMyValfn_current_N;
	array_2D fMyValfn_previous_N;
	array_2D_bool fMyPolicy_current;

	array_2D fMyfreq_lim_old;
	array_2D fMyfreq_lim_new;
	array_2D fMyN_lim_old;
	array_2D fMyN_lim_new;

	array_2D fMyN_map;
	array_2D fMyfreq_map;

	array_2D_bool fCheck;
	array_2D_int fIterMap;

	int fNum_dilution;
};

#endif // !CSOLVER_H

