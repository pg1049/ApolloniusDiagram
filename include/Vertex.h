#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include "Site.h"
#include "Cube.h"
#include "tool.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <list>
#include <iomanip>
#include <string>
#include <algorithm>
#include <time.h>
#include <functional>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef Site_3DApollonius::Site_3D    Site;
typedef Cube_3DApollonius::Cube       Cube;
using namespace std;
class Vertex_3DApollonius: public Cube_3DApollonius
{
public:
	struct range_vertex
	{
		double x_start;
		double x_end;
		double y_start;
		double y_end;
		double z_start;
		double z_end;
		range_vertex() {};
		range_vertex(double a1, double b1, double a2, double b2,double a3,double b3)
		{
			x_start = a1;
			x_end   = b1;
			y_start = a2;
			y_end = b2;
			z_start = a3;
			z_end = b3;
		}
	};
struct Vertex_3D
{
	int id;
	int cube_id;
	Point vertex_3d;
	vector<Site> site_group;
	range_vertex vertex_range;
	bool flag;
};

struct VertexAndSite
{
	vector<Vertex_3D> vg;
	vector<Site> sg;
	bool flag = 1;
};
struct Bound_vertex
{
	Vertex_3D v;
	vector<vector<Site>> sgg;
};
struct Event
{
	Cube cube_one;
	double distance;
	friend bool operator<(const Event& event1, const Event& event2)
	{
		return event1.distance > event2.distance;       //优先级按照有大到小排列
	}

};
protected:
vector<Cube> m_cube_result;
//vector<Cube> m_box;
vector<Cube> m_cube_bound;
vector<Cube> m_cube_bound_result;
vector<Vertex_3D> m_vertex;

//friend ostream& operator<<(ostream& out, Vertex_3D& vertex);
//friend bool operator==(const Vertex_3D  V1, const Vertex_3D V2);
void vertex_3D_get();
void vertex_remove();
Cube update_cube_one(Site m_site, Cube cube,bool& flag);
void first_update_cube();
void update_cube();
void compute_vertex();
void remove_cube();
void refine_cube(double &lenth_flag);
void confirm_vertex();
vector<double> solve(Cube &cube);
vector<double> solve1(Cube &cube);
vector<Vertex_3D> GetIntersectionOfVertex(vector<Vertex_3D> vg1,vector<Vertex_3D> vg2);
public:
	void out_vertex( const char* filename3);
};