#pragma once
#ifndef MYDEFINITIONS_H
#define MYDEFINITIONS_H


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


// my own assertion function in release mode
inline void myAssert(bool b) {
	if (!b)
	{
		std::cerr << "My Assertion Failed!" << "\n";
		std::abort();
	}
}

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

/*==============================================================================================
* This function uses binary search to return the first index such that xloc <= a sorted grid
* aValfnArray (input) : a 1D boost array of a sorted uniform grid
* aXloc (input): the query point that we want to find the index of which gridpoint is >= it.
*
* Right (output): the first integer-valued index such that aXloc <= gridpt at this index
*==============================================================================================*/
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

