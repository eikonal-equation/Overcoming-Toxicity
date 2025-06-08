#pragma once
#ifndef GRID_2D_H
#define CGRID_2D_H

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
  * File: CGrid_2D.h
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the declarations of the class "CGrid_2D"
  * that sets up a unifrom grid in 2D as well as functions for
  * bilinear interpolation and linear interpolation.
  *
  * Details of all of these functions are found in CGrid_2D.cpp.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "MyDefinitions.h"

class CGrid_2D
{
public:
	// constructors
	CGrid_2D() = default;
	CGrid_2D(double aDx, double aDy, double aXmin, double aYmin, int aXDim, int aYDim)
	{
		fDx = aDx; // Delta x (x: frequency of the killer in this project)
		fDy = aDy; // Delta y (y: total population in this project)
		fx0 = aXmin; // Xmin (default: 0 frequency)
		fy0 = aYmin; // Ymin (default: 0 population)
		fNx = aXDim; // Number of points in the x-direction
		fNy = aYDim; // Number of points in the y-direction
		f2h = 2.0 * fDx; // 2.0 * Delta x
		f3h = 3.0 * fDx; // 3.0 * Delta x
	}

	//---------------------------------(Bi)-Linear-Interpolation--------------------------------------------------------
	// Construct an array of length 2 for linear Interpolation (from a 2D time-slice of the value function
	// Returns array of length 2 containing the value of the 2 gridpoints for linear interpolation (in the x-direction)
	array<double, 2> Stencil_for_Linear_Interp(const array_2D& aValfnSlice,
		const int aXindex, const int aYindex);

	// Linear Interpolation in 1D
	// aValArray is the output array from Stencil_for_Linear_Interp
	double Linear_Interp(const array<double, 2>& aValArray, const int aYindex, const double aYloc);

	// Construct an array of length 4 for Bilinear Interpolation
	// Returns array of length 4 containing the value of the four gridpoints for bilinear interpolation
	array<double, 4> Stencil_for_Bilinear_Interp(const array_2D& aValfnSlice,
		const int aXindex, const int aYindex);

	// Bilinear Interpolation in 2D
	// aValArray is the output array from Stencil_for_Bilinear_Interp
	double Bilinear_Interp(const array<double, 4>& aValArray, const int aXindex, const int aYindex,
		const double aXloc, const double aYloc);

	// destructor
	~CGrid_2D() {};

//member variables of the class
protected:
	double fDx; // Delta x (x: frequency of the killer)
	double fDy; // Delta y (y: total population)
	double fx0; // Xmin (default: 0 frequency)
	double fy0; // Ymin (default: 0 population)
	int fNx; // Number of points in the x-direction
	int fNy; // Number of points in the y-direction
	double f2h; // 2.0 * Delta x
	double f3h; // 3.0 * Delta x
};


// This function computes a linear interpolation on a uniform 1D grid
// aValueArray(input) : a 2x1 double array of values of the 2-point stencil to be used for linear interpolation
// aXloc(input) : x-coordinate of the query point
// aXindex(input) : the first index such that the x-coordinate of the query point <= x, where x is the gridpoints in the x-direction
//
// Output: interpolated value
inline double CGrid_2D::Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc) {

	return aValArray[0] + (aValArray[1] - aValArray[0]) / fDx * (aXloc - (fx0 + (double(aIndex) - 1.0) * fDx));
}

// This function computes a bilinear interpolation on a uniform 2D grid
// aValueArray(input) : a 4x1 double array of values of the 4-point stencil to be used for bilinear interpolation
// aXloc(input) : x-coordinate of the query point
// aYloc(input) : y-coordinate of the query point
// aXindex(input) : the first index such that the x-coordinate of the query point <= x, where x is the gridpoints in the x-direction
// aYindex(input) : the first index such that the y-coordinate of the query point <= y, where y is the gridpoints in the y-direction
//
// Output: interpolated value
inline double CGrid_2D::Bilinear_Interp(const array<double, 4>& aValArray, const int aXindex, const int aYindex,
	const double aXloc, const double aYloc) {
	double x_pos_tr = aXindex * fDx + fx0;
	double y_pos_tr = aYindex * fDy + fy0;
	return (y_pos_tr - aYloc) / fDy * ((x_pos_tr - aXloc) / fDx * aValArray[0] + (aXloc - (x_pos_tr - fDx)) / fDx * aValArray[1])
		+ (aYloc - (y_pos_tr - fDy)) / fDy * ((x_pos_tr - aXloc) / fDx * aValArray[2] + (aXloc - (x_pos_tr - fDx)) / fDx * aValArray[3]);
}

#endif // ! CGrid_2D_H
