#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Site.h"

#include <iostream>
#include <vector>
#include <queue>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include "Eigen/Dense"
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SVD>




#define PI 3.1415926
using namespace std;
using namespace Eigen;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon2;

typedef Site_3DApollonius::Site_3D Site;

	struct Event_poly
	{
		Polygon2 poly;
		double area;
		friend bool operator<(const Event_poly& event1, const Event_poly& event2)
		{
			return event1.area < event2.area;
		}

	};
	Point_3 get_Fun_face_normal(Point_3 s1, Point_3 s2, Point_3 p);
	Point_3 get_normal(vector<Point_3> face_one);
	double dot_product(Point_3 P1, Point_3 P2);
	vector<double> cross_product(vector<double> v1,vector<double> v2);
	vector<double> solve_equation(vector<Site> m_site, Point_3 p);
	vector<double> solve_equation(vector<Site> m_site);
	vector<double> solve_equation(vector<Site> m_site, Point_3 p,vector<double> range);
	vector<int> GetInsectOfVector(vector<int> v1,vector<int> v2);
	Point_3 codeRotateByX(Point_3 p, double rx);
	Point_3 codeRotateByY(Point_3 p, double ry);
	Point_3 codeRotateByZ(Point_3 p, double rz);
	void RotateXYZ(vector<double> v1, vector<double> &theta, vector<Point_3>& boundary_point);
	void ReRotateXYZ(vector<double> v1, vector<double> v2, vector<Point_3>& boundary_point);
	vector<Point_2> RandomPointInPolygon(Polygon2 poly,int num);
	Point_2 compute2DPolygonCentroid(Polygon2 poly);
	double Distance_3D(Point_3 p1, Point_3 p2);
	double Distance_2D(Point_2 p1, Point_2 p2);
	double MaxLenthOfTriangle(Point_2 p1,Point_2 p2,Point_2 p3);
	Point_3 RotateByVector(Point_3 old_point, Point_3  vn, double theta);
	void delay(int seconds);//ÑÓ³Ùº¯Êý
