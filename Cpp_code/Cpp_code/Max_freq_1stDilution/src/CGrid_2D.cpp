/*==============================================================================
 * File: CGrid_2D.cpp
 * Author: MingYi, Lesley, Vishal, Nicolas
 *
 * Description:
 *============================================================================*/
#include "CGrid_2D.h"
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>
using boost::math::interpolators::makima;
using boost::math::interpolators::pchip;

//---------------------------------(Bi)-Linear-Interpolation--------------------------------------------------------
//stencil for boundary condition j=0 and j=fN
array<double, 2> CGrid_2D::Stencil_for_Linear_Interp(const array_2D& aValfnSlice,
	const int aXindex, const int aYindex) {
	array<double, 2> theValArray;
	theValArray[1] = aValfnSlice[aXindex][aYindex];// top
	theValArray[0] = aValfnSlice[aXindex - 1][aYindex];//bottom
	return theValArray;
}




 //Stencil_for_Bilinear_Interp function returns an array of 4 double that stores the value at
 // topright,topleft, bottomright, and bottomleft gridpoints.
 // taking aIndex as the index of gridpoints on the right.
 //aIndex is the index of i in marching
array<double, 4> CGrid_2D::Stencil_for_Bilinear_Interp(const array_2D& aValfnSlice,
	const int aXindex, const int aYindex) {
	array<double, 4> theValArray;

	theValArray[0] = aValfnSlice[aXindex - 1][aYindex - 1];// bottom left
	theValArray[1] = aValfnSlice[aXindex][aYindex - 1];// bottom right
	theValArray[2] = aValfnSlice[aXindex - 1][aYindex];//top left
	theValArray[3] = aValfnSlice[aXindex][aYindex];//top right

	//based on this definition, if aValmat[i][j], i is the column, j is the row?
	return theValArray;
}


// This function computes the second divided difference
// D2(input) : the second 0-th order divided difference
// D12(input) : the first 1st order divided difference
// aValueArray(input) : 3x1 vector of corresponding values of sample points
//
// D3(output) : the third 0-th order divided difference
// D23(output) : the second 1st order divided difference
// D123(output) : the first 2nd order divided difference
//
template <typename ValType, size_t Dim>
tuple<double, double, double>
CGrid_2D::Second_Divided_Difference(const double D2, const double D12, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray)
{
	// D123 : = f[x1, x2, x3]
	//        = (f[x2, x3] - f[x1, x2]) / (x3 - x1).
	double D3 = aValueArray[2];
	double D23 = (D3 - D2) / (aStencilArray[2] - aStencilArray[1]);
	double D123 = (D23 - D12) / (aStencilArray[2] - aStencilArray[0]);
	return make_tuple(D3, D23, D123);
}

// This function computes the third divided difference
// D3(input) : the third 0-th order divided difference
// D23(input) : the second 1st order divided difference
// D123(input) : the first 2nd order divided difference
// aValueArray(input) : 4x1 vector of corresponding values of sample points
//
// D1234(output) : the first 3rd order divided difference
//
template <typename ValType, size_t Dim>
double CGrid_2D::Third_Divided_Difference(const double D3, const double D23, const double D123, const array<ValType, Dim>& aValueArray, const array<ValType, Dim>& aStencilArray)
{
	// D1234 : = f[x1, x2, x3, x4]
	//         = (f[x2, x3, x4] - f[x1, x2, x3]) / (x4 - x1).
	double D4 = aValueArray[3];
	double D34 = (D4 - D3) / (aStencilArray[3] - aStencilArray[2]);
	double D234 = (D34 - D23) / (aStencilArray[3] - aStencilArray[1]);
	double D1234 = (D234 - D123) / (aStencilArray[3] - aStencilArray[0]);
	return D1234;
}

// This function computes cubic interpolation in Newton form
// aStencilArray(input) : 4x1 vector of sample points
// aCoefficientArray(input) : 4x1 coefficient vector of x computed from Newton divided difference.
// xloc(input) : coordinate of the query point
// val(output) : interpolated value
template <typename ValType, size_t Dim>
double CGrid_2D::NewtonInterp(const array<ValType, Dim>& aStencilArray, const array<ValType, Dim>& aCoefficientArray, const double aXloc)
{
	const int n = aCoefficientArray.size() - 1; //degree of the interpolating polynomial
	double val = aCoefficientArray[n];
	for (int i = 0; i < n; i++)
	{
		val = aCoefficientArray[n - i - 1] + (aXloc - aStencilArray[n - i - 1]) * val;
	}
	return val;
}

// This function computes cubic interpolation for a given stencil
// aStencilArray(input) : 4x1 vector of sample points
// aValueArray(input) : 4x1 vector of values at sample points
// xloc(input) : coordinate of the query point
// val(output) : interpolated value
double CGrid_2D::CubicInterp(const array<double, 4>& aStencilArray, const array<double, 4>& aValueArray, const double aXloc)
{
	const int degree = aStencilArray.size() - 1; //degree of the interpolating polynomial

	//compute 1st DD
	double D1 = aValueArray[0];
	double D2 = aValueArray[1];
	double D12 = (D2 - D1) / fDx;

	//compute 2nd DD
	double D3, D23, D123;
	std::tie(D3, D23, D123) = Second_Divided_Difference(D2, D12, aValueArray, aStencilArray);

	//compute 3rd DD
	double D1234 = Third_Divided_Difference(D3, D23, D123, aValueArray, aStencilArray);

	array<double, 4> aCoefficientArray{ D1,D12,D123,D1234 };

	// Newton Interpolation
	double val = NewtonInterp(aStencilArray, aCoefficientArray, aXloc);

	return val;
}

// This function computes quadratic interpolation for a given stencil
// aStencilArray(input) : 4x1 vector of sample points
// aValueArray(input) : 4x1 vector of values at sample points
// xloc(input) : coordinate of the query point
// val(output) : interpolated value
double CGrid_2D::QuadInterp(const array<double, 3>& aStencilArray, const array<double, 3>& aValueArray, const double aXloc)
{
	const int degree = aStencilArray.size() - 1; //degree of the interpolating polynomial

	//compute 1st DD
	double D1 = aValueArray[0];
	double D2 = aValueArray[1];
	double D12 = (D2 - D1) / fDx;

	//compute 2nd DD
	double D3, D23, D123;
	std::tie(D3, D23, D123) = Second_Divided_Difference(D2, D12, aValueArray, aStencilArray);

	array<double, 3> aCoefficientArray{ D1,D12,D123 };

	// Newton Interpolation
	double val = NewtonInterp(aStencilArray, aCoefficientArray, aXloc);

	return val;
}