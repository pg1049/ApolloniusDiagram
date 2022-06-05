#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
using namespace std;

class Site_3DApollonius {
public:
	struct Site_3D
	{
		Point site_position;
		double weight;
		int id;
		Site_3D() {};
		Site_3D(Point p, double w = 0) :site_position(p), weight(w) { id = -1; }
		bool operator<(const Site_3D& other) const
		{
			return make_pair(weight, site_position) < make_pair(other.weight, other.site_position);
		}
		bool operator>(const Site_3D& other) const
		{
			return make_pair(weight, site_position) > make_pair(other.weight, other.site_position);
		}

		bool operator==(const Site_3D& other) const
		{
			return make_pair(weight, site_position) == make_pair(other.weight, other.site_position);
		}
		double Weighted_Distance(const Point& pt) const
		{
			return sqrt((site_position - pt).squared_length()) - weight;
		}
	};
	int num_virture_site = 0;

public:
	int virturre_site_size;
	vector<Site_3D> m_sites;
	vector<double> m_boundary;
	vector<double> m_cube_boundary;
	double weight_max;
//	double ebsion = 1e-5;
	void Receive_sites(const char* filename);
	//friend bool operator==(vector<Site_3D> site1, vector<Site_3D> site2);
	//friend bool operator==(const Site_3D site1, const Site_3D site2);
	//friend bool operator==(const Point  P1, const Point P2);
	void add_virtual_site(double boundary);
	void GetIntersectOfSitegroup(vector<Site_3D> &sg1, vector<Site_3D> &sg2, vector<Site_3D> &sg);
	bool IfSubset(vector<Site_3D> &sg1, vector<Site_3D> &sg2);
	bool IfEquality(vector<Site_3D> &sg1, vector<Site_3D> &sg2);
	bool IfCoplanarity(vector<Site_3D> sg);
	bool IfCollineation(vector<Site_3D> sg);
};