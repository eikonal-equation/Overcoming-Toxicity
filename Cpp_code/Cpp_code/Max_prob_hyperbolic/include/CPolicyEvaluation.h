#pragma once
#ifndef CPOLICYEVALUATION_H
#define CPOLICYEVALUATION_H


//#include "CBarycentricLagrange.h"

//----------------------------Libraries------------------------------
#include "MyDefinitions.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iomanip>

using namespace Eigen;
class CPolicyEvaluation
{
public:
	// constructors
	CPolicyEvaluation() = default;
	CPolicyEvaluation(double aXDim, double aYDim, double aX0, double aY0, double aDx, double aDy,
		 double aDefeat_thres, double aVictory_thres, double aProb_not_arrival, double aProb_arrival)
	{
		fNx = aXDim;
		fNy = aYDim;
		fx0 = aX0;
		fy0 = aY0;

		fDx = aDx;
		fDy = aDy;

		fDefeat_thres = aDefeat_thres;
		fVictory_thres = aVictory_thres;


		fProb_arrival = aProb_arrival;
		fProb_not_arrival = aProb_not_arrival;

	};


	void PolicyEvaluation(const array_2D_bool& aPolicySlice,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat, array_2D& aValfn);


	void ApproximatePolicyEvaluation(const array_2D_bool& aPolicySlice,
		const array_4D& afullstate_loc_mat, const array_4D_int& afullstate_ind_mat,
		array_2D& aValfn, double tolerance);

	~CPolicyEvaluation() {};

private:
	//static paramters
	double fDefeat_thres;
	double fVictory_thres;
	double fProb_not_arrival;
	double fProb_arrival;

	//Discretization and definitions
	int fNx; //number of grid points for x
	int fNy; //number of grid points for y
	double fDx;
	double fDy;
	double fx0;
	double fy0;


};
#endif // !CPOLICYEVALUATION_H

