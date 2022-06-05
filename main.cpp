#pragma comment(lib, "libmpfr-4.lib")
#pragma comment(lib, "libgmp-10.lib")
#include <Windows.h>
#include <fstream>
#include <MA\quadrature.hpp>
#include <CGAL\Polygon_2.h>
#include <Eigen\core>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL\Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Hyperbola_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Hyperbola_2.h>
#include <vector>
//#include "..\\PQP\Distance_OBB.hpp"
//#include "..\\CVTOMTCPD\ConstructSparseCDT.hpp"
//#include "..\\Model3D\RichModel.h"
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\Model3D.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\Model3D.lib")
//#endif
//#include "..\\Geodesic\Xin_Wang.h"
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\Geodesic.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\Geodesic.lib")
//#endif
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\LBFGS.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\LBFGS.lib")
//#endif
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\PQP.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\PQP.lib")
//#endif
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\SolveAxb.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\SolveAxb.lib")
//#endif
//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\CVTOMTCPD.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\CVTOMTCPD.lib")
//#endif
#include <iostream>
#include <vector>
using namespace std;


#include <cgal/Hyperbola_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL\Apollonius_graph_traits_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;

struct Site_3D
{
	double weight;
	Point_3 site_position;
	Site_3D() {};
	Site_3D(Point_3 p, double w = 0) :site_position(p), weight(w) {}
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
	double Weighted_Distance(const Point_3& pt) const
	{
		return sqrt((site_position - pt).squared_length()) - weight;
	}
};

struct Cube
{
	Point_3 center;
	Point_3 corners[8];
	double sideLength;
	Cube() {};
	Cube(Point_3 center, double sidelen)
		: center(center), sideLength(sidelen)
	{
		corners[0] = center + 0.5 *sideLength* Vector_3(-1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, -1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, -1);
		corners[1] = center + 0.5 *sideLength* Vector_3(1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, -1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, -1);
		corners[2] = center + 0.5 *sideLength* Vector_3(-1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, 1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, -1);
		corners[3] = center + 0.5 *sideLength* Vector_3(-1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, -1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, 1);
		corners[4] = center + 0.5 *sideLength* Vector_3(1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, 1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, -1);
		corners[5] = center + 0.5 *sideLength* Vector_3(-1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, 1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, 1);
		corners[6] = center + 0.5 *sideLength* Vector_3(1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, -1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, 1);
		corners[7] = center + 0.5 *sideLength* Vector_3(1, 0, 0)
			+ 0.5 *sideLength* Vector_3(0, 1, 0)
			+ 0.5 *sideLength* Vector_3(0, 0, 1);
	}

	bool operator<(const Cube& other) const
	{
		return make_pair(center, sideLength) < make_pair(other.center, other.sideLength);
	}

	bool operator>(const Cube& other) const
	{
		return make_pair(center, sideLength) > make_pair(other.center, other.sideLength);
	}

	bool operator==(const Cube& other) const
	{
		return make_pair(center, sideLength) == make_pair(other.center, other.sideLength);
	}
};

bool DefeatTheOther(const Site_3D& f_large, const Site_3D& f_small)
{
	double distance = sqrt((f_large.site_position - f_small.site_position).squared_length());
	return f_large.weight - f_small.weight >= distance;
}

Point_3 MiddlePoint_without_hidden_sites(const Site_3D& f1, const Site_3D& f2)
{
	double distance = sqrt((f1.site_position - f2.site_position).squared_length());
	double d1 = (distance + f1.weight - f2.weight) * 0.5;
	double d2 = (distance - f1.weight + f2.weight) * 0.5;
	return CGAL::ORIGIN + d2 / (d1 + d2) * (f1.site_position - CGAL::ORIGIN) + d1 / (d1 + d2) * (f2.site_position - CGAL::ORIGIN);
}

Plane_3 MiddlePlane_without_hidden_sites(const Site_3D& f1, const Site_3D& f2)
{
	return Plane_3(MiddlePoint_without_hidden_sites(f1, f2), f2.site_position - f1.site_position);
}

bool Colinear(const Point_3& pt1, const Point_3& pt2, const Point_3& pt3)
{
	return CGAL::cross_product(pt2 - pt1, pt3 - pt2).squared_length() < 1e-8;
}

Plane_3 MiddlePlane_without_hidden_sites(const Site_3D& f_small, const Site_3D& f_large, const Point_3& pt)
{
	auto cross = CGAL::cross_product(f_large.site_position - f_small.site_position, pt - f_large.site_position);
	if (cross.squared_length() < 1e-8)
		return Plane_3(MiddlePoint_without_hidden_sites(f_small, f_large), f_large.site_position - f_small.site_position);
	Vector_3 dir1 = f_large.site_position - f_small.site_position;
	dir1 = 1.0 / sqrt(dir1.squared_length()) * dir1;
	Vector_3 dir2 = CGAL::cross_product(cross, dir1);
	dir2 = 1.0 / sqrt(dir2.squared_length()) * dir2;
	//check dir1 and dir2...
	K::Point_2 site2d_small((f_small.site_position - CGAL::ORIGIN) * dir1, (f_small.site_position - CGAL::ORIGIN) * dir2);
	K::Point_2 site2d_large((f_large.site_position - CGAL::ORIGIN) * dir1, (f_large.site_position - CGAL::ORIGIN) * dir2);
	K::Point_2 pt2d((pt - CGAL::ORIGIN) * dir1, (pt - CGAL::ORIGIN) * dir2);
	CGAL::Hyperbola_2<CGAL::Apollonius_graph_traits_2<K>> hy(CGAL::Apollonius_graph_traits_2<K>::Site_2(site2d_small, f_small.weight), CGAL::Apollonius_graph_traits_2<K>::Site_2(site2d_large, f_large.weight));
	auto proj = hy.f(hy.t(pt2d)); //right?
	auto proj_3d = CGAL::ORIGIN + proj.x() * dir1 + proj.y() * dir2;
	auto dir_large = f_large.site_position - proj_3d;
	dir_large = 1.0 / sqrt(dir_large.squared_length()) * dir_large;
	auto dir_small = f_small.site_position - proj_3d;
	dir_small = 1.0 / sqrt(dir_small.squared_length()) * dir_small;
	Plane_3 plane(proj_3d, dir_large - dir_small);
	return plane;
}

int PKatCube_large_small(const Site_3D& f_large, const Site_3D& f_small, const Cube& cube, bool strict = false)
{
	static const double sqrt3 = sqrt(3.0);
	double weighted_distance_large = f_large.Weighted_Distance(cube.center);
	double weighted_distance_small = f_small.Weighted_Distance(cube.center);
	if (weighted_distance_large < weighted_distance_small - sqrt3 * cube.sideLength)
		return 1;
	if (weighted_distance_small < weighted_distance_large - sqrt3 * cube.sideLength)
		return -1;
	int couningOrentation(0);
	for (int i = 0; i < 8; ++i)
	{
		if (f_small.Weighted_Distance(cube.corners[i]) < f_large.Weighted_Distance(cube.corners[i]))
			couningOrentation++;
	}
	if (couningOrentation == 8)
		return -1;
	if (couningOrentation > 0)
		return 0;

	Plane_3 plane;
	if (strict)
		plane = MiddlePlane_without_hidden_sites(f_small, f_large, cube.center);
	else
		plane = MiddlePlane_without_hidden_sites(f_small, f_large);
	couningOrentation = 0;
	for (int i = 0; i < 8; ++i)
	{
		auto side = plane.oriented_side(cube.corners[i]);
		if (side == CGAL::ON_POSITIVE_SIDE
			|| side == CGAL::ON_ORIENTED_BOUNDARY)
			couningOrentation++;
	}
	if (couningOrentation == 8)
		return 1;
	return 0;
}

int PKatCube(const Site_3D& f1, const Site_3D& f2, const Cube& cube, bool strict = false)
{
	if (f1.weight > f2.weight)
		return PKatCube_large_small(f1, f2, cube, strict);
	return -PKatCube_large_small(f2, f1, cube, strict);
}

set<Site_3D> AddOneMoreSiteIntoCube(const set<Site_3D> &sites, const Site_3D& newSite, const Cube& cube, bool strict = false)
{
	set<Site_3D> winners(sites);
	set<Site_3D> to_be_removed;
	bool flag_add(false);
	for (auto winner : winners)
	{
		auto res = PKatCube(newSite, winner, cube, strict);
		if (res == 0)
		{
			flag_add = true;
		}
		else if (res == 1)
		{
			flag_add = true;
			to_be_removed.insert(winner);
		}
		else
		{
			flag_add = false;
			break;
		}
	}

	for (auto item : to_be_removed)
		winners.erase(item);
	if (flag_add)
		winners.insert(newSite);

	return winners;
}

set<Site_3D> PKatCube(const set<Site_3D> &sites, const Cube& cube, bool strict = false)
{
	if (sites.empty())
		return set<Site_3D>();
	set<Site_3D> winners;
	auto ptr = sites.begin();
	winners.insert(*ptr);
	ptr++;
	for (; ptr != sites.end(); ptr++)
	{
		auto newSite = *ptr;
		winners = AddOneMoreSiteIntoCube(winners, newSite, cube, strict);
	}
	return winners;
}

int main(int argc, char** argv)
{
	ifstream in("sites_cube.txt");
	char buf[256];
	Cube cube;
	set<Site_3D> sites;
	while (in.getline(buf, sizeof buf))
	{
		istringstream line(buf);
		string word;
		line >> word;
		if (word == "c")
		{
			Point_3 center;
			double sideLen;
			line >> center >> sideLen;
			cube = Cube(center, sideLen);
		}
		else if (word == "v")
		{
			Point_3 site;
			double weight;
			line >> site >> weight;
			sites.insert(Site_3D(site, weight));
		}
	}
	in.close();

	double t = GetTickCount();
	set<Site_3D> winners;
	//for (int i = 0; i < 1000; ++i)
	{
		winners = PKatCube(sites, cube, true);
	}
	t = GetTickCount() - t;
	t /= 1000.0;
	//t /= 1000;
	cerr << "total time: " << t << " seconds.\n";
	cerr << "=============winners:=============\n";
	for (auto winner : winners)
		cerr << winner.site_position << "\n";
	cerr << endl;
	//cerr << "=============losers:=============\n";
	//for (auto site : sites)
	//	if (winners.find(site) == winners.end())
	//		cerr << site .site_position << "\n";
	//cerr << endl;
}