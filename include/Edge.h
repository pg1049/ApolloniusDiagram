#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Site.h"
#include "Vertex.h"
#include "tool.h"
#include "VectorOperations.h"
#include<iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef Site_3DApollonius::Site_3D  Site;
typedef Vertex_3DApollonius::Vertex_3D Vertex;
typedef CGAL::Polygon_2<K> Polygon_2;

typedef pair<int, double> PAIR;

using namespace std;
class Edge_3DApollonius :public Vertex_3DApollonius

{
public:
	struct Edge_3D
	{
		vector<Site> focus;
		Site site;
		Point start;
		Point end;
		int start_id;
		int end_id;
		vector<Point>  line;
		vector<Point>  reline;
		Point normal;
		int id;

	};
	bool cmp_by_value(const PAIR& lhs, const PAIR& rhs) {
		return lhs.second < rhs.second;
	}

	struct CmpByValue {
		bool operator()(const PAIR& lhs, const PAIR& rhs) {
			return lhs.second < rhs.second;
		}
	};
	vector<Edge_3D> m_edges;
protected:

	map<vector<int>, vector<int>> site_vertex;
	friend ostream& operator<<(ostream& out, Edge_3D edge);
	friend bool operator==( Edge_3D hyperbola1,  Edge_3D hyperbola2);
	void GetEdgeEndPoint();
	void insert_poinBetweenEndpoint(Edge_3D& hyperbola1, double step);
	void get_hyperbola(double step);
	void remove_pointofhyperbola(Edge_3D&  edge, vector<double>& boundary);
	Edge_3D creatEdge(vector<Site> sg, vector<Vertex> vg);
	map<int, vector<int>> GetCellboundary( );
};
