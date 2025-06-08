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
  * File: CGrid_2D.cpp
  *
  * Author: MingYi Wang
  *
  * Description: This file contains the acutal implementation of member functions
  * that sets up a unifrom grid in 2D as well as functions for
  * bilinear interpolation and linear interpolation.
  *
  *============================================================================*/

 //----------------------Project specific header files---------------------------
#include "CGrid_2D.h"


//---------------------------------(Bi)-Linear-Interpolation--------------------------------------------------------
//
// Extarct the 1D stencil from a 2D time-slice of the value function 
// for linear interpolation in the x-direction.
// It takes a 2D time-slice of the value function, 
// the x and y-index (in the spatial grid) of the query point 
// and returns an array of 2 double that stores the value at
// left and right gridpoints with the same x-index.
array<double, 2> CGrid_2D::Stencil_for_Linear_Interp(const array_2D& aValfnSlice,
	const int aXindex, const int aYindex) {
	array<double, 2> theValArray;
	theValArray[1] = aValfnSlice[aXindex][aYindex]; // left
	theValArray[0] = aValfnSlice[aXindex - 1][aYindex]; // right
	return theValArray;
}



// Create the 2D rectangular stencil from a 2D time-slice of the value function 
// for bi-linear interpolation.
// It takes a 2D time-slice of the value function, 
// the x and y-index (in the spatial grid) of the query point 
// and returns an array of 4 double that stores the value at
// topright, topleft, bottomright, and bottomleft gridpoints.
// Order of Dimensions for aValfnSlice: x-dimension (frequency of the killer), y-dimension (total population)
array<double, 4> CGrid_2D::Stencil_for_Bilinear_Interp(const array_2D& aValfnSlice,
	const int aXindex, const int aYindex) {
	array<double, 4> theValArray;

	theValArray[0] = aValfnSlice[aXindex - 1][aYindex - 1];// bottom left
	theValArray[1] = aValfnSlice[aXindex][aYindex - 1];// bottom right
	theValArray[2] = aValfnSlice[aXindex - 1][aYindex];//top left
	theValArray[3] = aValfnSlice[aXindex][aYindex];//top right

	return theValArray;
}
