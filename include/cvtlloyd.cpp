#include "cvtlloyd.h"

CVT_lloyd::CVT_lloyd(vector<Point_2> m_sites, Polygon_2 _boundary, int num)
{
	for (int i = 0; i < m_sites.size(); ++i)
	{
		m_new_vd.insert(m_sites.at(i));
	}
	for (int i = 0; i < _boundary.size(); ++i)
	{
		m_boundary_site.push_back(_boundary[i]);
		m_new_vd.insert(_boundary[i]);
	}
	m_boundary.out_boundary = _boundary;
	add_virture_vertice(m_new_vd);
	for (VD::Bounded_faces_iterator fb = m_new_vd.bounded_faces_begin(); fb != m_new_vd.bounded_faces_end(); fb++)
	{
		auto site = fb->dual()->point();
		m_original_site.push_back(site);
		Ccb_halfedge_circulator ec_start = fb->ccb();
		Ccb_halfedge_circulator ec = ec_start;
		Polygon_2 partition_one;
		do {
			partition_one.push_back(ec->source()->point());
		} while (++ec != ec_start);
		Polygon_2 CVT_partition_one;
		CVT_partition_one = Intersect(partition_one, m_boundary.out_boundary);
		Point_2 p(site.hx(), site.hy());
		bool flag = true;
		for (int i = 0; i < m_boundary_site.size(); i++)
		{
			Point_2 v = m_boundary_site.at(i);
			if (p == v)
				flag = false;
		}
		if (flag == true)
		{
			Point_2 center = compute2DPolygonCentroid(CVT_partition_one);
			m_new_site.push_back(center);
			m_original_cell[site] = CVT_partition_one;
			m_new_cell[center] = CVT_partition_one;
		}
		else
		{

			m_boundary_partition[site] = CVT_partition_one;
		}

	}
	//draw_domain("E://code//Apollonius_new//result//domain.obj");
	//draw_CVT_to_obj("E://code//Apollonius_new//result//cvt.obj");
	//getchar();
	for (int i = 0; i < num; ++i)
	{
		update_partition();
	}

};


CVT_lloyd::CVT_lloyd(VD vd, Polygon_2 _boundary, int num)
{
	m_boundary.out_boundary = _boundary;
	add_virture_vertice(vd);
	for (VD::Bounded_faces_iterator fb = vd.bounded_faces_begin(); fb != vd.bounded_faces_end(); fb++)
	{
		auto site = fb->dual()->point();
		m_original_site.push_back(site);
		Ccb_halfedge_circulator ec_start = fb->ccb();
		Ccb_halfedge_circulator ec = ec_start;
		Polygon_2 partition_one;
		do {
			partition_one.push_back(ec->source()->point());
		} while (++ec != ec_start);
		Polygon_2 CVT_partition_one;
		CVT_partition_one = Intersect(partition_one, m_boundary.out_boundary);
		Point_2 center = compute2DPolygonCentroid(CVT_partition_one);
		m_new_site.push_back(center);
		m_original_cell[site] = CVT_partition_one;
		m_new_cell[center] = CVT_partition_one;
	}
	for (int i = 0; i < num; ++i)
	{
		update_partition();
	}

};
void CVT_lloyd::print_endpoint(Halfedge_handle e, bool is_src)
{
	std::cout << "\t";
	if (is_src) {
		if (e->has_source())  std::cout << e->source()->point() << std::endl;
		else  std::cout << "point at infinity" << std::endl;
	}
	else {
		if (e->has_target())  std::cout << e->target()->point() << std::endl;
		else  std::cout << "point at infinity" << std::endl;
	}
}

void CVT_lloyd::add_virture_vertice(VD &vd)
{

	double x_min = m_boundary.out_boundary.bbox().xmin();
	double x_max = m_boundary.out_boundary.bbox().xmax();
	double y_min = m_boundary.out_boundary.bbox().ymin();
	double y_max = m_boundary.out_boundary.bbox().ymax();
	double x_lenth = x_max - x_min;
	double y_lenth = y_max - y_min;
	Point_2 p1(x_min - x_lenth, y_min - y_lenth);
	Point_2 p2(x_min - x_lenth, y_max + y_lenth);
	Point_2 p3(x_max + x_lenth, y_min - y_lenth);
	Point_2 p4(x_max + x_lenth, y_max + y_lenth);
	vd.insert(p1);
	vd.insert(p2);
	vd.insert(p3);
	vd.insert(p4);

}

Polygon_2 CVT_lloyd::Intersect(const Polygon_2& poly1, const Polygon_2& poly2)
{
	gpc_vertex_list vertextList1;
	vertextList1.num_vertices = poly1.size();
	std::vector<gpc_vertex> vecVertex1;
	for (int i = 0; i < vertextList1.num_vertices; i++)
	{
		gpc_vertex tm;
		tm.x = poly1[i].x();
		tm.y = poly1[i].y();
		vecVertex1.push_back(tm);
	}
	vertextList1.vertex = &(vecVertex1[0]);
	gpc_polygon gpc_polygon1;
	gpc_polygon1.num_contours = 1;
	gpc_polygon1.hole = NULL;
	gpc_polygon1.contour = &vertextList1;

	gpc_vertex_list vertextList2;
	vertextList2.num_vertices = poly2.size();
	std::vector<gpc_vertex> vecVertex2;
	for (int i = 0; i < vertextList2.num_vertices; i++)
	{
		gpc_vertex tm;
		tm.x = poly2[i].x();
		tm.y = poly2[i].y();
		vecVertex2.push_back(tm);
	}
	vertextList2.vertex = &(vecVertex2[0]);
	gpc_polygon gpc_polygon2;
	gpc_polygon2.num_contours = 1;
	gpc_polygon2.hole = NULL;
	gpc_polygon2.contour = &vertextList2;

	gpc_polygon result;
	gpc_polygon_clip(GPC_INT, &gpc_polygon1, &gpc_polygon2, &result);
	if (result.contour != NULL)
	{
		int numver = result.contour[0].num_vertices;
		vector<Point_2> resultpoint;
		if (numver > 0)
		{
			for (int m = numver - 1; m >= 0; m--)
			{
				Point_2 t(result.contour[0].vertex[m].x, result.contour[0].vertex[m].y);
				resultpoint.push_back(t);
			}
			Polygon_2 tmp(resultpoint.begin(), resultpoint.end());
			return tmp;
		}
		else
			return poly2;
	}
	else
		return poly2;

}


void CVT_lloyd::draw_CVT_to_matlab(string filename)
{

	ofstream out(filename);
	out << "clc;clear;close;\n";
	out << "figure(1);\n";
	map<Point_2, Polygon_2>::iterator  it;
	for (it = m_new_cell.begin(); it != m_new_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		out << "plot(" << site.x() << "," << site.y() << ",\'b*\');" << endl;
		out << "hold on;" << endl;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			out << "x=[" << eit->source().hx() << "," << eit->target().hx() << "];\n";
			out << "y=[" << eit->source().hy() << "," << eit->target().hy() << "];\n";
			out << "plot(x,y,\'r.-\');\nhold on;\n";
			//	out << "x=[" << eit->source().hx() << "," << site.x() << "];\n";
			//	out << "s=[" << eit->source().hy() << "," << site.y() << "];\n";
			//	out << "plot(x,s,\'b.-\');\nhold on;\n";
		}
	}
	out << "axis equal" << endl;
	out << "axis off;" << endl;
	out.close();
};
void CVT_lloyd::update_partition()
{
	VD vd;
	for (int i = 0; i < m_boundary_site.size(); ++i)
	{
		vd.insert(m_boundary_site[i]);
	}
	m_old_site.clear();
	for (int i = 0; i < m_new_site.size(); ++i)
	{
		double  x = m_new_site.at(i).x();
		double  y = m_new_site.at(i).y();
		if (x > -DBL_MAX&&x<DBL_MAX&&y>-DBL_MAX&&y < DBL_MAX)
		{
			CVTSite_2 t(m_new_site.at(i).x(), m_new_site.at(i).y());
			vd.insert(t);
			m_old_site.push_back(t);
		}
	}
	m_new_vd = vd;
	add_virture_vertice(m_new_vd);
	m_new_cell.clear();
	m_new_site.clear();
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
		Point_2 p(site.hx(), site.hy());
		bool flag = true;
		for (int i = 0; i < m_boundary_site.size(); i++)
		{
			Point_2 v = m_boundary_site.at(i);
			if (p == v)
				flag = false;
		}
		if (flag == true)
		{
			Point_2 center = compute2DPolygonCentroid(CVT_partition_one);
			m_new_cell[center] = CVT_partition_one;
			m_new_site.push_back(center);
		}
		else
		{
			m_boundary_partition[site] = CVT_partition_one;
		}

		//	m_new_partition[site] = partition_one;
	}

	//cout << endl;
}

Point_2 CVT_lloyd::compute2DPolygonCentroid(Polygon_2 poly)
{
	double ctx = 0;
	double cty = 0;
	double signedArea = 0.0;
	double x0 = 0.0; // Current vertex X
	double y0 = 0.0; // Current vertex Y
	double x1 = 0.0; // Next vertex X
	double y1 = 0.0; // Next vertex Y
	double a = 0.0;  // Partial signed area

					 // For all vertices
	int vertexCount = poly.size();
	for (int i = 0; i < vertexCount; ++i)
	{
		x0 = poly[i].x();
		y0 = poly[i].y();
		x1 = poly[(i + 1) % vertexCount].x();
		y1 = poly[(i + 1) % vertexCount].y();
		a = x0*y1 - x1*y0;
		signedArea += a;
		ctx += (x0 + x1)*a;
		cty += (y0 + y1)*a;
	}

	signedArea *= 0.5;
	ctx /= (6.0*signedArea);
	cty /= (6.0*signedArea);
	Point_2 centroid(ctx, cty);
	return centroid;
}
void CVT_lloyd::GetVertex()
{
	map<Point_2, Polygon_2>::iterator  it;
	for (it = m_new_cell.begin(); it != m_new_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		//
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			Point_2 p1(eit->source().hx(), eit->source().hy());
			bool flag1 = true;
			for (int i = 0; i < m_vertex.size(); ++i)
			{
				Point_2 p = m_vertex.at(i);
				if (p == p1)
				{
					flag1 = false;
					break;
				}
			}

			if (flag1 == true)
			{
				m_vertex.push_back(p1);
			}

		}
	}
}

void CVT_lloyd::draw_CVT_to_obj(string filename)
{
	set<Point_2> vertex;
	ofstream out(filename);
	vector<int>  line;
	map<Point_2, int> vertex_to_id;
	map<Point_2, Polygon_2>::iterator  it;
	set<Point_2>::iterator it_set;
	for (it = m_new_cell.begin(); it != m_new_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			if (vertex.find(eit->source()) == vertex.end())
			{
				vertex.insert(eit->source());
			}
		}
	}
	int id = 0;
	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		Point_2 p(it1->hx(), it1->hy());
		vertex_to_id[p] = id;
		id++;
	}
	for (it = m_new_cell.begin(); it != m_new_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			int index1 = vertex_to_id.find(eit->source())->second;
			int index2 = vertex_to_id.find(eit->target())->second;
			line.push_back(index1);
			line.push_back(index2);
		}
	}

	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	for (int i = 0; i < line.size() / 2; i++)
	{
		out << "l " << line.at(i * 2) + 1 << "  " << line.at(i * 2 + 1) + 1 << endl;
	}


	out.close();
};

void CVT_lloyd::draw_CVT_site_to_obj(string filename)
{
	ofstream out(filename);
	for (auto it1 = m_new_site.begin(); it1 != m_new_site.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	out.close();
};

void CVT_lloyd::draw_CVT_old_site_obj(string filename)
{
	ofstream out(filename);
	for (auto it1 = m_old_site.begin(); it1 != m_old_site.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	out.close();
};
void CVT_lloyd::draw_orginal_CVT_to_obj(string filename)
{
	set<Point_2> vertex;
	ofstream out(filename);
	vector<int>  line;
	map<Point_2, int> vertex_to_id;
	map<Point_2, Polygon_2>::iterator  it;
	set<Point_2>::iterator it_set;
	for (it = m_original_cell.begin(); it != m_original_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			if (vertex.find(eit->source()) == vertex.end())
			{
				vertex.insert(eit->source());
			}
		}
	}
	int id = 0;
	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		Point_2 p(it1->hx(), it1->hy());
		vertex_to_id[p] = id;
		id++;
	}
	for (it = m_original_cell.begin(); it != m_original_cell.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			int index1 = vertex_to_id.find(eit->source())->second;
			int index2 = vertex_to_id.find(eit->target())->second;
			line.push_back(index1);
			line.push_back(index2);
		}
	}

	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	for (int i = 0; i < line.size() / 2; i++)
	{
		out << "l " << line.at(i * 2) + 1 << "  " << line.at(i * 2 + 1) + 1 << endl;
	}


	out.close();

}
void CVT_lloyd::draw_original_site_to_obj(string filename)
{
	set<Point_2> vertex;
	ofstream out(filename);
	vector<int>  line;
	map<Point_2, int> vertex_to_id;
	map<Point_2, Polygon_2>::iterator  it;
	set<Point_2>::iterator it_set;
	for (it = m_original_cell.begin(); it != m_original_cell.end(); it++)
	{
		Point_2 site = it->first;
		if (vertex.find(it->first) == vertex.end())
		{
			vertex.insert(it->first);
		}
	}

	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	out.close();

}

void CVT_lloyd::draw_domain(string filename)
{
	ofstream out(filename);
	vector<Point_2> vertex;
	vector<int>  line;
	for (auto eit = m_boundary.out_boundary.edges_begin(); eit != m_boundary.out_boundary.edges_end(); eit++)
	{
		auto point = eit->source();
		vertex.push_back(point);
	}
	for (auto eit = m_boundary.out_boundary.edges_begin(); eit != m_boundary.out_boundary.edges_end(); eit++)
	{
		vector<Point_2>::iterator iter1 = std::find(vertex.begin(), vertex.end(), eit->source());
		vector<Point_2>::iterator iter2 = std::find(vertex.begin(), vertex.end(), eit->target());
		int index1 = std::distance(vertex.begin(), iter1);
		int index2 = std::distance(vertex.begin(), iter2);
		line.push_back(index1);
		line.push_back(index2);
	}
	for (auto it1 = vertex.begin(); it1 != vertex.end(); it1++)
	{
		out << "v " << it1->hx() << " " << it1->hy() << " " << 0 << endl;
	}
	for (int i = 0; i < line.size() / 2; i++)
	{
		out << "l " << line.at(i * 2) + 1 << "  " << line.at(i * 2 + 1) + 1 << endl;
	}
	out.close();
}

