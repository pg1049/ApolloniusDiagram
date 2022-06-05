#include "cvt.h"
#include "gpc.h"
#include <MA/quadrature.hpp>
#include "Function2D.h"
//#include "stdafx.h"
void CVT::print_endpoint(Halfedge_handle e, bool is_src)
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

void CVT::add_virture_vertice()
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
	m_new_vd.insert(p1);
	m_new_vd.insert(p2);
	m_new_vd.insert(p3);
	m_new_vd.insert(p4);



}

Polygon_2 CVT::Intersect(const Polygon_2& poly1, const Polygon_2& poly2)
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
	if (result.contour == NULL)
		return poly2;
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

void CVT::draw_bounded_voronoi_to_matlab(char* filename)
{
	ofstream out(filename);
	out << "clc;clear;close;\n";
	out << "figure(1);\n";
	map<Point_2, Polygon_2>::iterator  it;
	for (it = m_partition.begin(); it != m_partition.end(); it++)
	{
		Point_2 site = it->first;
		Polygon_2 cell = it->second;
		out << "plot(" << site.x() << "," << site.y() << ",\'b*\');" << endl;
		out << "hold on;" << endl;
		for (auto eit = cell.edges_begin(); eit != cell.edges_end(); eit++)
		{
			out << "x=[" << eit->source().hx() << "," << eit->target().hx() << "];\n";
			out << "y=[" << eit->source().hy() << "," << eit->target().hy() << "];\n";
			out << "plot(x,y,\'r*-\');\nhold on;\n";
		}
	}
	out << "axis equal" << endl;
	out << "axis off;" << endl;
	out.close();
};

void CVT::draw_CVT_boundary_to_matlab(string filename)
{
	ofstream out(filename);
	out << "clc;clear;close;\n";
	out << "figure(1);\n";
	for (int i = 0; i < m_boundary_site.size(); ++i)
	{
		int size = m_boundary_site.size();
		out << "x=[" << m_boundary_site.at(i).hx() << "," << m_boundary_site.at((i + 1) % size).hx() << "];\n";
		out << "y=[" << m_boundary_site.at(i).hy() << "," << m_boundary_site.at((i + 1) % size).hy() << "];\n";
		out << "plot(x,y,\'r.-\');\nhold on;\n";
	}
	out << "axis equal" << endl;
	out << "axis off;" << endl;
	out.close();
}

void CVT::draw_CVT_to_matlab(string filename)
{

	ofstream out(filename);
	out << "clc;clear;close;\n";
	out << "figure(1);\n";
	map<Point_2, Polygon_2>::iterator  it;
	for (it = m_new_partition.begin(); it != m_new_partition.end(); it++)
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
		}
	}
	for (it = m_boundary_partition.begin(); it != m_boundary_partition.end(); it++)
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
		}
	}
	out << "axis equal" << endl;
	out << "axis off;" << endl;
	out.close();
};

pair<double, vector<double>> CVT::ComputeEnergyAndGradient_CVT()
{
	vector<double> gradients(m_new_site.size() * 2, 0);
	double F(0);
	SquaredEuclideanDistance2D  dist;
	for (int i = 0; i < m_new_site.size(); i++)
	{
		Polygon_2 _boundary_one = m_new_partition.find(m_new_site.at(i))->second;
		F += MA::integrate_3<K::FT>(_boundary_one, K::FT(0), [&](Point_2 p)
		{
			return dist(p, m_new_site.at(i));
		});
		double f = MA::integrate_3<K::FT>(_boundary_one, K::FT(0), [&](Point_2 p)
		{
			return 1;
		});
		Point_2 centroid = compute2DPolygonCentroid(_boundary_one);
		K::Vector_2 diff = 2 * f*dist.GetParitial(m_new_site.at(i), centroid);
		gradients[i] += diff.x();
		gradients[i + m_new_site.size()] += diff.y();
	}

	for (int i = 0; i < m_boundary_site.size(); i++)
	{
		Polygon_2 _boundary_one = m_boundary_partition.find(m_boundary_site.at(i))->second;
		F += MA::integrate_3<K::FT>(_boundary_one, K::FT(0), [&](Point_2 p)
		{
			return dist(p, m_boundary_site.at(i));
		});
	}
	return make_pair(F, gradients);
}
void CVT::update_partition()
{

	VD vd;
	for (int i = 0; i < m_new_site.size(); ++i)
	{
		while (!(CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(m_boundary.out_boundary.vertices_begin(), m_boundary.out_boundary.vertices_end(), m_new_site.at(i), K())))
		{
			double x1 = m_new_site.at(i).hx();
			double x2 = m_old_site.at(i).hx();
			double y1 = m_new_site.at(i).hy();
			double y2 = m_old_site.at(i).hy();
		//	Point_2 p((x1 + x2) / 2, (y1 + y2) / 2);
			m_new_site.at(i) = m_old_site.at(i);
		}
		CVTSite_2 t(m_new_site.at(i).x(), m_new_site.at(i).y());
		vd.insert(t);
	}

	for (int i = 0; i < m_boundary_site.size(); ++i)
	{
		vd.insert(m_boundary_site.at(i));
	}
	m_new_vd = vd;
	add_virture_vertice();
	//cout << "number_of_faces()"<<vd.number_of_faces() << endl;
	m_new_partition.clear();
	m_boundary_partition.clear();
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
			m_new_partition[site] = CVT_partition_one;
		}
		else
		{
			m_boundary_partition[site] = CVT_partition_one;
		}
	}
	//cout << endl;
}

Point_2 CVT::compute2DPolygonCentroid(Polygon_2 poly)
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
void CVT::draw_CVT_dual_to_matlab(string filename)
{
	ofstream out(filename);
	out << "clc;clear;close;\n";
	out << "figure(1);\n";
	map<Point_2, Polygon_2>::iterator  it;
	for (it = m_new_partition.begin(); it != m_new_partition.end(); it++)
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
		}
	}
	for (it = m_boundary_partition.begin(); it != m_boundary_partition.end(); it++)
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
		}
	}
	Polygon_2 bound = m_boundary.out_boundary;
	/*VD m_dual_vd;
	for (int i = 0; i < m_boundary_site.size(); ++i)
	{
		m_dual_vd.insert(m_boundary_site.at(i));
	}
	for (int i = 0; i < m_old_site.size(); ++i)
	{
		m_dual_vd.insert(m_old_site.at(i));
	}
	for (auto vit = m_dual_vd.vertices_begin(); vit != m_dual_vd.vertices_end(); vit++)*/
	//	for (auto vit = m_new_vd.vertices_begin(); vit != m_new_vd.vertices_end(); vit++)
	//	{        
	//		Point_2 v = vit->point();
	////		out << "plot(" << v.x() << "," << v.y() << ",\'r*\');" << endl;
	//			auto fit = vit->dual();
	//			vector<Point_2> polypoints;			
	//			for (int i = 0; i < 3; i++)
	//			{
	//				polypoints.push_back(fit->vertex(i)->point());
	//				out << "x=[" << fit->vertex(i)->point().x() << "," << fit->vertex((i + 1) % 3)->point().x() << "];\n";
	//				out << "y=[" << fit->vertex(i)->point().y() << "," << fit->vertex((i + 1) % 3)->point().y() << "];\n";
	//				out << "plot(x,y,'color','k','LineWidth',1);\nhold on;\n";
	//			}	
	//	}	
	CDT cdt;
	for (int i = 0; i < m_boundary_site.size() - 1; ++i)
	{
		cdt.insert_constraint(m_boundary_site.at(i), m_boundary_site.at(i + 1));
		//	cdt.insert(m_boundary_site.at(i));
	}
	//	cdt.insert(m_boundary_site.at(m_boundary_site.size()-1));
	for (int i = 0; i < m_old_site.size(); ++i)
	{
		cdt.insert(m_old_site.at(i));
	}
	for (auto fit = cdt.faces_begin(); fit != cdt.faces_end(); fit++)
	{
		for (int i = 0; i < 3; i++)
		{
			out << "x=[" << fit->vertex(i)->point().x() << "," << fit->vertex((i + 1) % 3)->point().x() << "];\n";
			out << "y=[" << fit->vertex(i)->point().y() << "," << fit->vertex((i + 1) % 3)->point().y() << "];\n";
			out << "plot(x,y,'color','k','LineWidth',1);\nhold on;\n";
		}
	}
	out << "axis equal" << endl;
	out << "axis off;" << endl;
	out.close();
}
void CVT::draw_CVT_to_obj(string filename)
{
	set<Point_2> vertex;
	ofstream out(filename);
	vector<int>  line;
	map<Point_2, int> vertex_to_id;
	map<Point_2, Polygon_2>::iterator  it;
	set<Point_2>::iterator it_set;
	for (it = m_new_partition.begin(); it != m_new_partition.end(); it++)
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
	for (it = m_boundary_partition.begin(); it != m_boundary_partition.end(); it++)
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
	for (it = m_new_partition.begin(); it != m_new_partition.end(); it++)
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
	for (it = m_boundary_partition.begin(); it != m_boundary_partition.end(); it++)
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
	for (int i = 0; i<line.size() / 2; i++)
	{
		out << "l " << line.at(i * 2) + 1 << "  " << line.at(i * 2 + 1) + 1 << endl;
	}


	out.close();
};
void CVT::draw_site_to_obj(string filename)
{
	ofstream out(filename);
	for (int i = 0; i < m_old_site.size();++i )
	{
		Point_2 p = m_old_site.at(i);
		out << "v " << p.x() << " " << p.y() << " " << 0 << endl;
	}
	for (int i = 0; i< m_boundary.out_boundary.size(); ++i)
	{
		Point_2 p = m_boundary.out_boundary[i];
		out << "v " << p.x() << " " << p.y() << " " << 0 << endl;
	}
		out.close();
}