#pragma once
#ifndef GRID_2D_H
#define CGRID_2D_H

/*==============================================================================
 * File: CGrid_2D.h
 *
 * Author: MingYi, Lesley, Vishal, Nicolas
 *
 * Description:
 *
 *============================================================================*/
#include "MyDefinitions.h"

class CGrid_2D
{
public:
	//constructors
	CGrid_2D() = default;
	CGrid_2D(double aDx, double aDy, double aXmin, double aYmin, int aXDim, int aYDim)
	{
		fDx = aDx;
		fDy = aDy;
		fx0 = aXmin;
		fy0 = aYmin;
		fNx = aXDim;
		fNy = aYDim;
		f2h = 2.0 * fDx;
		f3h = 3.0 * fDx;
	}

	//---------------------------------(Bi)-Linear-Interpolation--------------------------------------------------------
	array<double, 2> Stencil_for_Linear_Interp(const array_2D& aValfnSlice,
		const int aXindex, const int aYindex);

	double Linear_Interp(const array<double, 2>& aValArray, const int aYindex, const double aYloc);

	array<double, 4> Stencil_for_Bilinear_Interp(const array_2D& aValfnSlice,
		const int aXindex, const int aYindex);

	//Bilinear Interpolation in 2D
	double Bilinear_Interp(const array<double, 4>& aValArray, const int aXindex, const int aYindex,
		const double aXloc, const double aYloc);


	//-------------------------------------------Cubic-Interpolation---------------------------------------------------
	// This function computes the second divided difference
	template <typename ValType, size_t Dim>
	tuple<double, double, double> Second_Divided_Difference(const double D2, const double D12, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray);

	// This function computes the third divided difference
	template <typename ValType, size_t Dim>
	double Third_Divided_Difference(const double D3, const double D23, const double D123, const array<ValType, Dim>& aValueArray,const array<ValType, Dim>& aStencilArray);

	// This function is the cubic interpolation in Newton form given a 4-point stencil
	template <typename ValType, size_t Dim>
	double NewtonInterp(const array<ValType, Dim>& aStencilArray, const array<ValType, Dim>& aCoefficientArray, const double aXloc);

	// This function computes the cubic interpolation for a given stencil
	double CubicInterp(const array<double, 4>& aStencilArray, const array<double, 4>& aValueArray, const double aXloc);

	// This function computes the quadratic interpolation for a given stencil
	double QuadInterp(const array<double, 3>& aStencilArray, const array<double, 3>& aValueArray, const double aXloc);

	//destructor
	~CGrid_2D() {};

	//member variables of the class
protected:
	double fDx; // Delta x or h
	double fDy;
	double fx0; // starting position
	double fy0;
	int fNx; // numeber of points of the position array, assume fN_x
	int fNy;
	double f2h;
	double f3h;
};


//linear interpolation function for boundary condition j=0 and j=fN
inline double CGrid_2D::Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc) {

	return aValArray[0] + (aValArray[1] - aValArray[0]) / fDx * (aXloc - (fx0 + (double(aIndex) - 1.0) * fDx));
}
// reference: return aValArray[0] + (aValArray[1] - aValArray[0]) / fDx * (aXloc - (fx0 + (double(aIndex) - 1.0) * fDx));
//aValArray[0] is the "bottom"


//aXindex & aYindex are the index for the topright grid for bilinear interpolation
inline double CGrid_2D::Bilinear_Interp(const array<double, 4>& aValArray, const int aXindex, const int aYindex,
	const double aXloc, const double aYloc) {
	double x_pos_tr = aXindex * fDx + fx0;
	double y_pos_tr = aYindex * fDy + fy0;
	return (y_pos_tr - aYloc) / fDy * ((x_pos_tr - aXloc) / fDx * aValArray[0] + (aXloc - (x_pos_tr - fDx)) / fDx * aValArray[1])
		+ (aYloc - (y_pos_tr - fDy)) / fDy * ((x_pos_tr - aXloc) / fDx * aValArray[2] + (aXloc - (x_pos_tr - fDx)) / fDx * aValArray[3]);
}

#endif // ! CGrid_2D_H
