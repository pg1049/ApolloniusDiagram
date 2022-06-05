#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
#include "Site.h"
#include <cgal/Hyperbola_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL\Apollonius_graph_traits_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
typedef Site_3DApollonius::Site_3D  Site_3D;
using namespace std;
class Cube_3DApollonius: public Site_3DApollonius
{
public:
	struct Cube
	{
		Point indexOfcube;
		Point pointOfcube_cente;
		double lenthOfcube_sides;
		vector<pair<Site_3D, double>> sites;
		double min_distance;
		vector<Point> cube_vertice;
		int id;
		bool inOrbound;
		bool empty;
		bool operator<(const Cube& other) const
		{
			return make_pair(pointOfcube_cente, lenthOfcube_sides) < make_pair(other.pointOfcube_cente, other.lenthOfcube_sides);
		}

		bool operator>(const Cube& other) const
		{
			return make_pair(pointOfcube_cente, lenthOfcube_sides) > make_pair(other.pointOfcube_cente, other.lenthOfcube_sides);
		}

		bool operator==(const Cube& other) const
		{
			return make_pair(pointOfcube_cente, lenthOfcube_sides) == make_pair(other.pointOfcube_cente, other.lenthOfcube_sides);
		}
	};
	int numOfcube;
	int x_max_cube;
	int y_max_cube;
	int z_max_cube;
	Point source_point;
	vector<Cube> m_cube;
	int cmp(const pair<Site_3D, double>& x, const pair<Site_3D, double>& y)
	{
		return x.second < y.second;
	}
	friend ostream& operator<<(ostream& out, Cube& cube);
	double Distance_point2point(Point P1, Point P2);
	Point get_SiteInCube(Cube &cube_one, Point p_site);
	vector<Point> get_CubeNeighbor(Cube &cube_one);
	vector<Cube> segment_box(Cube &box);
	void Generate_Cube(vector<double>& boundary,int sites_size);
	void get_vertice_Cube(Cube &cube);
	Cube get_CubeByIndex(vector<Cube>& m_cube, Point index);
	void Cube_Inorbound(Cube &cube);

	bool DefeatTheOther(const Site_3D& f_large, const Site_3D& f_small)
	{
		double distance = sqrt((f_large.site_position - f_small.site_position).squared_length());
		return f_large.weight - f_small.weight >= distance;
	}

	Point MiddlePoint_without_hidden_sites(const Site_3D& f1, const Site_3D& f2)
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

	bool Colinear(const Point& pt1, const Point& pt2, const Point& pt3)
	{
		return CGAL::cross_product(pt2 - pt1, pt3 - pt2).squared_length() < 1e-8;
	}

	Plane_3 MiddlePlane_without_hidden_sites(const Site_3D& f_small, const Site_3D& f_large, const Point& pt)
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
		double weighted_distance_large = f_large.Weighted_Distance(cube.pointOfcube_cente);
		double weighted_distance_small = f_small.Weighted_Distance(cube.pointOfcube_cente);
		if (weighted_distance_large < weighted_distance_small - sqrt3 * cube.lenthOfcube_sides)
			return 1;
		if (weighted_distance_small < weighted_distance_large - sqrt3 * cube.lenthOfcube_sides)
			return -1;
		int couningOrentation(0);
		for (int i = 0; i < 8; ++i)
		{
			if (f_small.Weighted_Distance(cube.cube_vertice[i]) < f_large.Weighted_Distance(cube.cube_vertice[i]))
				couningOrentation++;
		}
		if (couningOrentation == 8)
			return -1;
		if (couningOrentation > 0)
			return 0;

		Plane_3 plane;
		if (strict)
			plane = MiddlePlane_without_hidden_sites(f_small, f_large, cube.pointOfcube_cente);
		else
			plane = MiddlePlane_without_hidden_sites(f_small, f_large);
		couningOrentation = 0;
		for (int i = 0; i < 8; ++i)
		{
			auto side = plane.oriented_side(cube.cube_vertice[i]);
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
};