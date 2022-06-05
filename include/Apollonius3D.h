#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include "Edge.h"
#include "cvt.h"
#include "cvtlloyd.h"
#include "CVT_LBFGS.h"
#include<random>
#include<ctime>
#include<cmath>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>


typedef CVT CVT;
typedef Site_3DApollonius::Site_3D Site;
typedef Edge_3DApollonius::Edge_3D Edge;
typedef Vertex_3DApollonius::Vertex_3D Vertex;
typedef CGAL::Polygon_2<K> Polygon2;
typedef K::Point_2 Point_2;


using namespace std;
using namespace  CGAL::internal;
using namespace CGAL::VoronoiDiagram_2;
class Apollonius3D : public Edge_3DApollonius
{
private:
	void ComputeFace();
	void MeshTheFace();
	void OrderTheFaceBoundary();
	void MapPointToOneFace(int i);
	void AddThePointToOneFace(int i);
	void ComputerCell();
	bool IfIntheFace(int faceid, Point p);
	
public:
	struct Face_3D
	{
		int id;
		int site1;
		int site2;
		//int site_id1;
		//int site_id2;
		vector<int> face_boundary;
		vector<Point> boundary;
		vector<int> boundary_index;
		vector<int> in_boundary_index;
		vector<Point>in_boundary;
	
		vector<Point> Points;
		vector<tuple<int, int, int>> connect;
	//	CDT cdt;
		Point normal;
	};
	struct Cell_3D
	{
		int id;
		int site_id;
		vector<int> face_ids;
		vector<int> site_ids;
		Point center;
	};
	Apollonius3D(char *inputfilename1);
	Apollonius3D(vector<Site_3D> sites);
	void run();
	Site GetSite(int id);
	Edge GetEdge(int id);
	Face_3D GetFace(int id);
	Vertex GetVertex(int id);
	Face_3D GetFace(int site1, int site2);
	void out_EdgeObj(char* file_edge);
	void out_FaceObj(string file_face);
	void out_faceone(string faceone, int id);
	void out_edgeone(string edgeone, int id);
	void out_PointObj(string file_point);
	void out_VertexObj(string file_vertex);
	void out_SiteObj(string file_site);
	void out_SiteText(string file_site);
	void out_SitePointObj(string file_site);
	void out_cell_one(string file_cell,int id);
	vector<Face_3D> m_faces;
	vector<Cell_3D> m_cells;
	map<int, vector<int>> cellboundary;
	vector<Point> m_points;
	int edge_point_endid;
	bool open;
	double  step;
	int getnumOfSite() { return m_sites.size(); }
	int getnumOfEdge() { return m_edges.size();}
	int getnumOfFace() { return m_faces.size(); }
	int getnumOfVertex() { return m_vertex.size();}
	int getPointId(Point p);
	
};