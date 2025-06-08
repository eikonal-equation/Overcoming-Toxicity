#pragma once
#ifndef MYDEFINITIONS_H
#define MYDEFINITIONS_H
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
 * File: MyDefinitions.h
 *
 * Author: MingYi Wang
 *
 * Description: This file contains definitions and helper functions 
 *
 *============================================================================*/

//-------------------------Libraires-------------------------------------------
#include <algorithm>
#include <numeric>
#include <tuple>
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <chrono>
#include <boost/multi_array.hpp>

//---------------------------Definitions---------------------------------------
using namespace std;
typedef boost::multi_array<double, 4> array_4D;
typedef boost::multi_array<int, 4> array_4D_int;
typedef boost::multi_array<double, 3> array_3D;
typedef boost::multi_array<int, 3> array_3D_int;
typedef boost::multi_array<bool, 3> array_3D_bool;
typedef boost::multi_array<double, 2> array_2D;
typedef boost::multi_array<bool, 2> array_2D_bool;
typedef boost::multi_array<int, 2> array_2D_int;
typedef boost::multi_array<double, 1> array_1D;


// ---------------------------Function Declarations--------------------------------
//My own assertion function
inline void myAssert(bool b) {
	if (!b)
	{
		std::cerr << "My Assertion Failed!" << "\n";
		std::abort();
	}
}


// This function finds the first index such that (the query point) xloc <= the uniform gridpoint on x-axis
// aXloc (input): the coordiante of the query point
// aDeltaX (input): the uniform grid spacing on x-axis
// aX0 (input): the left boundary (Xmin) of the uniform grid on x-axis
// aXmax (input): the right boundary (Xmax) of the uniform grid on x-axis
// aDim (input): the number of grid points on x-axis
//
// Output: the first integer-valued index such that xloc <= gridpt at this index
__inline int find_index(const double aXloc, const double aDeltaX, const double aX0, const double aXmax, const int aDim)
{
	//assert(aXloc <= aXmax);
	//assert(aXloc >= aX0);
	if (aXloc >= aXmax) {
		return aDim;
	}
	else {
		return int(std::ceil((aXloc + 1e-20 - aX0) / aDeltaX));
	}
}


// This function uses binary search to return the first index such that xloc <= a sorted grid
// aValfnArray (input) : a 1D boost array of a sorted uniform grid
// aXloc (input): the query point that we want to find the index of which gridpoint is >= it.
//
// right (output): the first integer-valued index such that aXloc <= gridpt at this index
inline int find_index(const array_1D& aArray, const double aLoc, const int aN) {
	int left = -1; //aArray[left] < x
	int right = aN + 1;
	while (right > left + 1) {
		int mid = (right + left) / 2;
		if (aArray[mid] < aLoc) {
			left = mid;
		}
		else {
			right = mid;
		}
	}
	return right;
}

#endif // !MYDEFINITIONS_H

