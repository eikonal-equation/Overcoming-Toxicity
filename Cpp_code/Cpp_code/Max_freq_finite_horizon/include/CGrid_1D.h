#pragma once
#ifndef GRID_1D_H
#define CGRID_1D_H

/*==============================================================================
 * File: CGrid_1D.h
 *
 * Author: MingYi, Lesley, Vishal, Nicolas
 *
 * Description: This file contains the declarations of a class that sets up
 * a unifrom grid in 1D as well as a function for linear interpolation.
 *
 * Details of all of these functions are found in CGrid_1D.cpp.
 *
 *============================================================================*/

 //---------------------------Libraries-----------------------------------------
#include <iostream>
#include "MyDefinitions.h"


class CGrid_1D
{
public:
	//constructors
	CGrid_1D() = default;
	CGrid_1D(double aDx, double aXmin, int aDim)
	{
		fDx = aDx;
		fx0 = aXmin;
		fN = aDim;
	}

	//construct an array with length 2 for Liner Interpolation
	array<double, 2> Stencil_for_Linear_Interp(const array_1D& aValfnArray, const int aIndex);

	//Linear Interpolation in 1D
	double Linear_Interp(const array<double, 2>& aValArray, const int aIndex, const double aXloc);

	//member variables of the class
protected:
	double fDx; // Delta x or h
	double fx0; // starting position
	int fN; // numeber of points of the position array

};
#endif // ! CGRID_1D_H
