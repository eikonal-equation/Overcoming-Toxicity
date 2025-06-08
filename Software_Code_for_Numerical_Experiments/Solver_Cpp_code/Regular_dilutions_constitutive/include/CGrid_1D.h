#pragma once
#ifndef GRID_1D_H
#define CGRID_1D_H

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
  * File: CGrid_1D.h
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the declarations of the class "CGrid_1D"
  * that sets up a unifrom grid in 1D as well as functions for linear interpolation.
  *
  * Details of all of these functions are found in CGrid_1D.cpp.
  *
  *============================================================================*/

//----------------------Project specific header files---------------------------
#include "MyDefinitions.h"
 //----------------------Project specific header files---------------------------
#include <iostream>


class CGrid_1D
{
public:
	//constructors
	CGrid_1D() = default;
	CGrid_1D(double aDx, double aXmin, int aDim)
	{
		fDx = aDx; // Delta x
		fx0 = aXmin; //Xmin (default: 0)
		fN = aDim; // Number of points in the x-direction
	}

	// Construct an array of length 2 for linear Interpolation from a 1D data array
	// Returns array of length 2 containing the value of the 2 gridpoints for linear interpolation 
	array<double, 2> Stencil_for_Linear_Interp(const array_1D& aValfnArray, const int aIndex);

	// Linear Interpolation in 1D
	// aValArray is the output array from Stencil_for_Linear_Interp
	double Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc);

	//member variables of the class
protected:
	double fDx; // Delta x
	double fx0; // Xmin (default: 0)
	int fN; // Number of points in the x-direction

};
#endif // ! CGRID_1D_H
