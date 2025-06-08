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
  * File: CGrid_1D.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the implementation of functions that extract
  * the stencil for linear interpolation, and implement linear interpolation
  *
  *============================================================================*/

  //----------------------Project specific header files---------------------------
#include "CGrid_1D.h"


// This function returns a stencil (array of size 2) for linear interpolation
// aValfnArray (input) : a 1D boost array of actual values on each grid point
// aIndex (input): the first index such that (the query point) xloc <= the gridpoint on x-axis
//
// theValArray (output): a std::array of size 2 that contains (the query point) xloc
array<double, 2> CGrid_1D::Stencil_for_Linear_Interp(const array_1D& aValfnArray, const int aIndex) 
{
	array<double, 2> theValArray;
	theValArray[0] = aValfnArray[aIndex - 1]; //left
	theValArray[1] = aValfnArray[aIndex]; //right

	return theValArray;
}


// This function implements linear interpolation in 1D
// aValArray (input) : a (size 2) std::array of actual values on the stencil
// aIndex (input): the first index such that (the query point) xloc <= the gridpoint on x-axis
// aXloc (input): the query point that we want to apply linear interpolation to
//
// Output: a (double) scalar that returns the interpolated value at xloc
double CGrid_1D::Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc) 
{
	return aValArray[0] + (aValArray[1] - aValArray[0]) / fDx * (aXloc - (fx0 + (double(aIndex) - 1.0) * fDx));
}
