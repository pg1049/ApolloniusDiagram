#pragma once
#include <vector>
#include <Eigen\sparse>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
using namespace std;
struct SingleValuedFunction
{
	virtual K::FT operator()(const vector<K::FT>& X) const = 0;
};

struct SingleValuedBinaryFunction
{
	virtual K::FT operator()(const Point_2& pt1, const Point_2& pt2) const = 0;
	virtual K::FT operator()(const vector<K::FT>& X1, const vector<K::FT>& X2) const
	{
		return (*this)(Point_2(X1[0], X1[1]), Point_2(X2[0], X2[1]));
	}
};
//
//struct SingleValuedFunction_Smooth_1Order : public SingleValuedFunction
//{
//	virtual vector<double> GetPrime(const vector<double>& X) const = 0;
//};
//
//struct SingleValuedFunction_Smooth_2Order_Dense : public SingleValuedFunction
//{
//	virtual Eigen::MatrixXd GetHessian(const vector<double>& X) const = 0;
//};
//
//struct SingleValuedFunction_Smooth_2Order_Sparse : public SingleValuedFunction
//{
//	virtual Eigen::SparseMatrix<double> GetHessian(const vector<double>& X) const = 0;
//};
//
//struct VectorValuedFunction
//{
//	virtual vector<double> operator()(const vector<double>& X) const = 0;
//};
//
//struct VectorValuedFunction_Smooth_1Order_Dense : public VectorValuedFunction
//{
//	virtual Eigen::MatrixXd GetPrime(const vector<double>& X) const = 0;
//};
//
//struct VectorValuedFunction_Smooth_1Order_Sparse : public VectorValuedFunction
//{
//	virtual Eigen::SparseMatrix<double> GetPrime(const vector<double>& X) const = 0;
//};

