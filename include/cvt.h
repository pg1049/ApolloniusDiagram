#pragma once
#define _CRT_SECURE_NO_WARNINGS
// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include <MA/functions.hpp>
#include <CGAL/Vector_2.h>
using namespace std;
// typedefs for defining the adaptor
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                    VD;
typedef CGAL::Polygon_2<K> Polygon_2;

typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;

// typedef for the result type of the point location
typedef AT::Site_2                    CVTSite_2;
typedef AT::Point_2                   Point_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;

typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef VD::Face_iterator  Face_iterator;
using namespace std;
class CVT
{
	struct  boundary
	{
		Polygon_2 in_boundary;
		Polygon_2 out_boundary;
	};
protected:
	vector<Point_2> m_site;
	boundary m_boundary;
	map<Point_2, Polygon_2> m_partition;
public:	
	VD m_new_vd;
	VD m_old_vd;
	vector<Point_2>  m_old_site;
	vector<Point_2>  m_new_site;
	vector<Point_2>  m_boundary_site;
	map<Point_2, Polygon_2> m_new_partition;
	map<Point_2, Polygon_2> m_boundary_partition;
	CVT() {};
	CVT(vector<Point_2> m_point, Polygon_2 _boundary)
	{
		VD vd;
		for (int i = 0; i < m_point.size(); ++i)
		{
			vd.insert(m_point.at(i));
		}
	
		for (int i = 0; i < _boundary.size(); ++i)
		{
			vd.insert(_boundary[i]);
		}
		m_new_vd = vd;
		m_old_vd = vd;
		m_boundary.out_boundary = _boundary;
		add_virture_vertice();
		for (VD::Bounded_faces_iterator fb = m_new_vd.bounded_faces_begin(); fb != m_new_vd.bounded_faces_end(); fb++)
		{
			auto site = fb->dual()->point();
			Ccb_halfedge_circulator ec_start = fb->ccb();
			Ccb_halfedge_circulator ec = ec_start;
			Polygon_2 partition_one;
			do {
				partition_one.push_back(ec->source()->point());
			} while (++ec != ec_start);
			Polygon_2 CVT_partition_one;
			CVT_partition_one = Intersect(partition_one, m_boundary.out_boundary);
			//CVT_partition_one = partition_one;
			m_partition[site] = CVT_partition_one;
			Point_2 p(site.hx(),site.hy());
			if(CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(m_boundary.out_boundary.vertices_begin(), m_boundary.out_boundary.vertices_end(), p, K()))
			{
				m_new_site.push_back(p);
				m_new_partition[site] = CVT_partition_one;
			}
			else
			{
				m_boundary_site.push_back(p);
				m_boundary_partition[site] = CVT_partition_one;
			}
		}
	
	};
	CVT(vector<Point_2> m_point, Polygon_2 _boundary, Polygon_2 _in_boundary)
	{
		VD vd;
		for (int i = 0; i < m_point.size(); ++i)
		{
			vd.insert(m_point.at(i));
		}

		for (int i = 0; i < _boundary.size(); ++i)
		{
			vd.insert(_boundary[i]);
		}
		for (int i = 0; i < _in_boundary.size(); ++i)
		{
			vd.insert(_in_boundary[i]);
		}
		m_new_vd = vd;
		m_old_vd = vd;
		m_boundary.out_boundary = _boundary;
		m_boundary.in_boundary = _in_boundary;
		add_virture_vertice();
		for (VD::Bounded_faces_iterator fb = m_new_vd.bounded_faces_begin(); fb != m_new_vd.bounded_faces_end(); fb++)
		{
			auto site = fb->dual()->point();
			Ccb_halfedge_circulator ec_start = fb->ccb();
			Ccb_halfedge_circulator ec = ec_start;
			Polygon_2 partition_one;
			do {
				partition_one.push_back(ec->source()->point());
			} while (++ec != ec_start);
			Polygon_2 CVT_partition_one;
			CVT_partition_one = Intersect(partition_one, m_boundary.out_boundary);
			//CVT_partition_one = partition_one;
			m_partition[site] = CVT_partition_one;
			Point_2 p(site.hx(), site.hy());
			if (CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(m_boundary.out_boundary.vertices_begin(), m_boundary.out_boundary.vertices_end(), p, K()))
			{
				m_new_site.push_back(p);
				m_new_partition[site] = CVT_partition_one;
			}
			if (CGAL::ON_BOUNDARY == CGAL::bounded_side_2(m_boundary.out_boundary.vertices_begin(), m_boundary.out_boundary.vertices_end(), p, K())
			|| CGAL::ON_BOUNDARY == CGAL::bounded_side_2(m_boundary.in_boundary.vertices_begin(), m_boundary.in_boundary.vertices_end(), p, K()))
			{
				m_new_site.push_back(p);
				m_new_partition[site] = CVT_partition_one;
			}
			{
				m_boundary_site.push_back(p);
				m_boundary_partition[site] = CVT_partition_one;
			}
		}

	};
	vector<Point_2> get_m_site()const { return m_site; };
	int get_m_site_size()const { return m_site.size(); };
    pair<double, vector<double>> ComputeEnergyAndGradient_CVT( );
	void print_endpoint(Halfedge_handle e, bool is_src);
	void add_virture_vertice();
	void draw_bounded_voronoi_to_matlab(char* filename);
	void draw_CVT_to_matlab(string filename);
	void draw_CVT_boundary_to_matlab(string filename);
	void draw_CVT_dual_to_matlab(string filename);
	void update_partition();
	Polygon_2 Intersect(const Polygon_2& poly1, const Polygon_2& poly2);
	Point_2 compute2DPolygonCentroid(Polygon_2 poly);
	void draw_CVT_to_obj(string filename);
	void draw_site_to_obj(string filename);
};

