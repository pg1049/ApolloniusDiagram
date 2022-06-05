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
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include "gpc.h"
#include <CGAL/Vector_2.h>
using namespace std;
// typedefs for defining the adaptor
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                    VD;
typedef CGAL::Polygon_2<K> Polygon_2;

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

class CVT_lloyd
{
	struct  boundary
	{
		Polygon_2  in_boundary;
		Polygon_2 out_boundary;
	};
protected:
	vector<Point_2> m_site;
	boundary m_boundary;
public:
	VD m_new_vd;
	VD m_old_vd;
	vector<Point_2>  m_original_site;
	vector<Point_2>  m_old_site;
	vector<Point_2>  m_new_site;
	vector<Point_2>  m_boundary_site;
	map<Point_2, Polygon_2> m_boundary_partition;
	map<Point_2, Polygon_2> m_new_cell;
	map<Point_2, Polygon_2> m_original_cell;
	vector<Point_2> m_vertex;
	CVT_lloyd() {};
	CVT_lloyd(vector<Point_2> m_sites, Polygon_2 _boundary, int num);
	CVT_lloyd(VD vd, Polygon_2 _boundary, int num);
	vector<Point_2> get_m_site()const { return m_site; };
	int get_m_site_size()const { return m_site.size(); };
	void print_endpoint(Halfedge_handle e, bool is_src);
	void add_virture_vertice(VD &vd);
	void draw_domain(string filename);
	void draw_CVT_to_matlab(string filename);
	void draw_CVT_to_obj(string filename);
	void draw_CVT_old_site_obj(string filename);
	void draw_CVT_site_to_obj(string filename);
	void draw_orginal_CVT_to_obj(string filename);
	void draw_original_site_to_obj(string filename);
	void update_partition();
	Polygon_2 Intersect(const Polygon_2& poly1, const Polygon_2& poly2);
	Point_2 compute2DPolygonCentroid(Polygon_2 poly);
	void GetVertex();
};
