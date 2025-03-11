/*==============================================================================
 * File: CGrid_1D.cpp
 * Author: MingYi, Lesley, Vishal, Nicolas
 *
 * Description: This file contains the implementations of functions that set up a
 * stencil for linear interpolation as well as the actual computation of
 * linear interpolation.
 *============================================================================*/
#include "CGrid_1D.h"


/*==============================================================================================
* This function returns a stencil (array of size 2) for linear interpolation
* aValfnArray (input) : a 1D boost array of actual values on each grid point
* aIndex (input): the first index such that (the query point) xloc <= the gridpoint on x-axis
*
* theValArray (output): a std::array of size 2 that contains (the query point) xloc
*==============================================================================================*/
array<double, 2> CGrid_1D::Stencil_for_Linear_Interp(const array_1D& aValfnArray, const int aIndex) {
	array<double, 2> theValArray;
	theValArray[0] = aValfnArray[aIndex - 1];
	theValArray[1] = aValfnArray[aIndex];

	return theValArray;
}


/*==============================================================================================
* This function implements linear interpolation in 1D
* aValArray (input) : a (size 2) std::array of actual values on the stencil
* aIndex (input): the first index such that (the query point) xloc <= the gridpoint on x-axis
* aXloc (input): the query point that we want to apply linear interpolation to
*
* Output: a (double) scalar that returns the interpolated value at xloc
*==============================================================================================*/
double CGrid_1D::Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc) {
	//return max(((fx0 + (aIndex)*fDx - aXloc) * aValArray[0] + (aXloc - (fx0 + (aIndex - 1.0) * fDx)) * aValArray[1])/fDx, 0.0);
	return aValArray[0] + (aValArray[1] - aValArray[0]) / fDx * (aXloc - (fx0 + (double(aIndex) - 1.0) * fDx));
}
