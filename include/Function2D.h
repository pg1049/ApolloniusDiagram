#pragma once
//#include "Polygon_2_Operation.h"
#include <boost/shared_ptr.hpp>
#include "Function.h"
#include <cmath>
typedef  K::Vector_2  Vector_2;
using namespace std;
struct SingleValuedFunction2D : public SingleValuedFunction
{
	virtual K::FT operator()(const Point_2& pt) const = 0;
	virtual K::FT operator()(const vector<K::FT>& X) const
	{
		return (*this)(Point_2(X[0], X[1]));
	}
};

struct SingleValuedBinaryFunction2D : public SingleValuedBinaryFunction
{
	virtual K::FT operator()(const Point_2& pt1, const Point_2& pt2) const = 0;
	virtual K::FT operator()(const vector<K::FT>& X1, const vector<K::FT>& X2) const
	{
		return (*this)(Point_2(X1[0], X1[1]), Point_2(X2[0], X2[1]));
	}
};

//struct Distortion : public SingleValuedFunction2D
//{
//	virtual K::FT operator()(const Point_2& pt) const
//	{
//		double x = pt.x() + 0.2 * abs(pt.y()) * pt.y();
//		double y = pt.y() + 0.2 * abs(pt.x()) * pt.x();
//		return abs(x) - abs(y);
//	}
//};

//struct SingleValuedFunction_Smooth_1Order2D : public SingleValuedFunction2D
//{
//	virtual vector<K::FT> GetPrime(const Point_2& pt) const = 0;
//	vector<K::FT> GetPrime(const vector<double>& X) const
//	{
//		return this->GetPrime(Point_2(X[0], X[1]));
//	}
//};
//
//struct SingleValuedFunction_Smooth_2Order_Dense2D : public SingleValuedFunction_Smooth_1Order2D
//{
//	virtual Eigen::MatrixXd GetHessian(const Point_2& pt) const = 0;
//	virtual Eigen::MatrixXd GetHessian(const vector<double>& X) const
//	{
//			return this->GetHessian(Point_2(X[0], X[1]));
//	}
//};
//
//struct SingleValuedFunction_Smooth_2Order_Sparse2D : public SingleValuedFunction_Smooth_1Order2D
//{
//	virtual Eigen::SparseMatrix<double> GetHessian(const Point_2& pt) const = 0;
//	virtual Eigen::SparseMatrix<double> GetHessian(const vector<double>& X) const
//	{
//			return this->GetHessian(Point_2(X[0], X[1]));
//	}
//};
//
//struct VectorValuedFunction2D : public VectorValuedFunction
//{
//	virtual vector<K::FT> operator()(const Point_2& X) const = 0;
//	virtual vector<K::FT> operator()(const vector<double>& X) const
//	{
//		return (*this)(Point_2(X[0], X[1]));
//	}
//};
//
//struct VectorValuedFunction_Smooth_1Order_Dense2D : public VectorValuedFunction2D
//{
//	virtual Eigen::MatrixXd GetPrime(const Point_2& X) const = 0;
//	virtual Eigen::MatrixXd GetPrime(const vector<double>& X) const
//	{
//		return this->GetPrime(Point_2(X[0], X[1]));
//	}
//};
//
//struct VectorValuedFunction_Smooth_1Order_Sparse2D : public VectorValuedFunction2D
//{
//	virtual Eigen::SparseMatrix<K::FT> GetPrime(const Point_2& X) const = 0;
//	virtual Eigen::SparseMatrix<K::FT> GetPrime(const vector<double>& X) const
//	{
//		return this->GetPrime(Point_2(X[0], X[1]));
//	}
//};

struct Rho2D_Uniform : public SingleValuedFunction2D
{
	Rho2D_Uniform()
	{
	}
	virtual K::FT operator()(const Point_2& pt) const
	{
		return 1.0;
	}
};
struct SingleValuedBinaryFunction_Differentiable : public SingleValuedBinaryFunction
{
	//differential at pt_var
	virtual K::Vector_2 GetParitial(const Point_2& pt_var, const Point_2& pt) const = 0;
};

struct SquaredEuclideanDistance2D : public SingleValuedBinaryFunction_Differentiable
{
	virtual K::FT operator()(const Point_2& pt1, const Point_2& pt2) const
	{
		return (pt1 - pt2).squared_length();
	}
	virtual K::Vector_2 GetParitial(const Point_2& pt_var, const Point_2& pt) const
	{
		return 2 * K::Vector_2(pt_var.x() - pt.x(), pt_var.y() - pt.y());
	}
	virtual K::Vector_2 GetParitial_Without_DetaX(const Point_2& pt_var, const Point_2& pt) const
	{
		return K::Vector_2(2, 2);
	}
};


typedef SquaredEuclideanDistance2D L2Distance2D;

struct StraightLineDistance : public SingleValuedBinaryFunction_Differentiable
{
	virtual K::FT operator()(const Point_2& pt1, const Point_2& pt2) const
	{
		K::FT detaX = pt1.x() - pt2.x();
		K::FT detaY = pt1.y() - pt2.y();
		detaX *= detaX;
		detaY *= detaY;
		return sqrt(detaX + detaY);
	}
	virtual K::Vector_2 GetParitial(const Point_2& pt_var, const Point_2& pt) const
	{
		K::FT detaX = pt_var.x() - pt.x();
		K::FT detaY = pt_var.y() - pt.y();
		K::FT norm_inf = max(abs(detaX), abs(detaY));
		detaX /= norm_inf;
		detaY /= norm_inf;
		K::FT denom = sqrt(detaX * detaX + detaY * detaY);
		return K::Vector_2(detaX / denom, detaY / denom);
	}
	virtual K::Vector_2 GetParitial_Without_DetaX(const Point_2& pt_var, const Point_2& pt) const
	{
		K::FT detaX = pt_var.x() - pt.x();
		K::FT detaY = pt_var.y() - pt.y();
		K::FT _sqrt = sqrt(detaX * detaX + detaY * detaY);

		return K::Vector_2(1.0 / _sqrt, 1.0 / _sqrt);
	}
};

struct X_4_Y_4 : public SingleValuedBinaryFunction_Differentiable
{
	virtual K::FT operator()(const Point_2& pt1, const Point_2& pt2) const
	{
		K::FT detaX = pt1.x() - pt2.x();
		K::FT detaY = pt1.y() - pt2.y();
		detaX *= detaX;
		detaX *= detaX;
		detaY *= detaY;
		detaY *= detaY;
		return detaX + detaY;
	}
	virtual Vector_2 GetParitial(const Point_2& pt_var, const Point_2& pt) const
	{
		K::FT detaX = pt_var.x() - pt.x();
		K::FT detaY = pt_var.y() - pt.y();
		detaX *= detaX * detaX;
		detaY *= detaY * detaY;
		return 4 * Vector_2(detaX, detaY);
	}
	virtual Vector_2 GetParitial_Without_DetaX(const Point_2& pt_var, const Point_2& pt) const
	{
		K::FT detaX = pt_var.x() - pt.x();
		K::FT detaY = pt_var.y() - pt.y();
		detaX *= detaX;
		detaY *= detaY;
		return 4 * Vector_2(detaX, detaY);
	}
};


//struct Rho2D_Normal_Distribution : public SingleValuedFunction2D
//{
//	double m_deviation;
//	double m_scale;
//	double m_maxDisToBoundary;
//	Point m_center; 
//	Rho2D_Normal_Distribution()
//	{
//		m_center = CGAL::ORIGIN;
//		m_deviation = 0;
//		m_scale = 0;
//		m_maxDisToBoundary = FLT_MAX;
//	}
//	Rho2D_Normal_Distribution(const Polygon_2& boundary,  double deviation) 
//		: m_deviation(deviation)
//	{
//		m_center = GetCentroid(boundary);
//		m_maxDisToBoundary = FindClosestPoint(boundary, m_center).first;
//		m_scale = 1;
//		UpdateScale(boundary);
//	}
//	virtual double operator()(const Point& pt) const 
//	{
//		double dis = sqrt((m_center - pt).squared_length()) / m_maxDisToBoundary;
//		return exp(-dis * m_deviation) * m_scale;
//	}
//protected:	
//	void UpdateScale(const Polygon_2& boundary);
//};
