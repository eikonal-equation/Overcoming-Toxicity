#pragma once
#ifndef CSOLVER_H
#define CSOLVER_H
#include "MyDefinitions.h"
#include "CGrid_2D.h"
#include "CPolicyEvaluation.h"
class PDE_Solver
{
public:
	PDE_Solver() = default;
	PDE_Solver(int aRefinement_factor, double aMaxPopSize, double aMaxFrequency,
		double aTol, double aDefeat_thres, double aVictory_thres,
		double aEpsilon, double aRks, double aGamma, double aLamb, double aRho)
	{
		fNx = aRefinement_factor * 100;
		fNy = aRefinement_factor * 100;
		fymax = aMaxPopSize;
		fxmax = aMaxFrequency;
		fx0 = 0.0;
		fy0 = 0.0;
		fDx = (fxmax - fx0) / fNx;
		fDy = (fymax - fy0) / fNy;

		fDefeat_thres = aDefeat_thres;
		fVictory_thres = aVictory_thres;

		fEpsilon = aEpsilon;
		fRks = aRks;
		fGamma = aGamma;
		fLamb = aLamb;
		fRho = aRho;
		fDelta = 0.7;

		fMyValfn.resize(boost::extents[fNx + 1][fNy + 1]);
		fMyPolicy.resize(boost::extents[fNx + 1][fNy + 1]);

		fTau = sqrt(fDx);
		fMaxIter = 100000;
		fTol = aTol;

		fProb_not_arrival = exp(-fLamb * fTau);
		fProb_arrival = 1.0 - fProb_not_arrival;
	}

	// drift function f
	array<double, 2> f_drift(const double aX, const double aY, const int aCtrl);

	// Initialization of matrices which store value function values, optimal policies, and switchgrids
	void InitializeMat();

	// This function pre-compute all the coefficients and foot of characteristics
	void precompute_coeff(array_4D& fullstate_loc_mat, array_4D_int& fullstate_ind_mat);

	//This function computes the GS Value iterations for a single gridpt
	array<double,4> ValueIte_single_point(const int aI, const int aJ,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);


	//This function computes the GS Value iterations for a single gridpt
	double ValueIte_single_step(const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, CGrid_2D& aGrid);

	//This is the main solver with GS VI
	void Main_VPI_Solver();

	void writeDomainToFile(const string aFilename);

	~PDE_Solver() { cout << "Destructor called." << "\n"; };

private:
	//static paramters
	double fEpsilon;
	double fRks;
	double fGamma;
	double fLamb;
	double fRho;
	double fDelta;
	double fDefeat_thres;
	double fVictory_thres;
	double fProb_not_arrival;
	double fProb_arrival;

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
	long fMaxIter;
	double fTol;
	array_2D fMyValfn;
	array_2D_bool fMyPolicy;
};

#endif // !CSOLVER_H

