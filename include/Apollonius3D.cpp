#define _CRT_SECURE_NO_WARNINGS
#include "Apollonius3D.h"

typedef Apollonius3D::Face_3D Face;
using namespace  std;

Apollonius3D::Apollonius3D(char *inputfilename1)
{
	Receive_sites(inputfilename1);
	if (m_sites.size() <= 0)
	{
		cout << "the input sites is NULL" << endl;
		exit(0);
	}
	add_virtual_site(1);
}
Apollonius3D::Apollonius3D(vector<Site_3D> sites)
{
	m_sites = sites;
	for (int i = 0; i < m_sites.size(); i++)
	{
		m_sites[i].id = i;
	}
	if (m_sites.size() <= 0)
	{
		cout << "the input sites is NULL" << endl;
		exit(0);
	}
	add_virtual_site(1);
}
void Apollonius3D::run()
{
	
	clock_t starttime, endtime;
	starttime = clock();
	vertex_3D_get();
	endtime = clock();
	cout << "vertex time:" << (double)(endtime - starttime) / CLOCKS_PER_SEC << "s" << endl;
	
	int m_vertex_size = m_vertex.size();
	cout << m_vertex_size << endl;
	if (m_vertex_size > 1)
	{
		//for (int i = 0; i < m_vertex.size(); ++i)
		//{
		//	Point p = m_vertex.at(i).vertex_3d;
		//	if(p.hx()>=m_boundary.at(0)&& p.hx() <= m_boundary.at(1)
		//		&& p.hy() >= m_boundary.at(2) && p.hy() <= m_boundary.at(3)
		//		&& p.hz() >= m_boundary.at(4) && p.hz() <= m_boundary.at(5))
		//	m_points.push_back(p);
		//}
		clock_t starttime1, endtime1;
		starttime1 = clock();
		double x_lenth = m_boundary.at(1) - m_boundary.at(0);
		double y_lenth = m_boundary.at(3) - m_boundary.at(2);
		double z_lenth = m_boundary.at(5) - m_boundary.at(4);
		double _step = sqrt(x_lenth * x_lenth + y_lenth * y_lenth + z_lenth * z_lenth) /20;
		step = _step;
		get_hyperbola(step);
		endtime1 = clock();
		cout << "edge time:" << (double)(endtime1 - starttime1) / CLOCKS_PER_SEC << "s" << endl;
		int m_edges_size = m_edges.size();
		cout << "m_edges_size"<<m_edges_size << endl;

//	 m_edges_size = 0;
		if (m_edges_size > 1)
		{
			for (int i = 0; i < m_edges.size(); ++i)
			{
				vector<Point> line = m_edges.at(i).line;
				for (int j = 0; j < line.size(); ++j)
					m_points.push_back(line.at(j));
			}
			edge_point_endid = m_points.size() - 1;
			clock_t starttime2, endtime2;
			starttime2 = clock();
			cellboundary = GetCellboundary();
			ComputeFace();
			MeshTheFace();
			endtime2 = clock();
			cout << "face time:" << (double)(endtime2 - starttime2) / CLOCKS_PER_SEC << "s" << endl;
			ComputerCell();
		}
	}
}
void Apollonius3D::out_EdgeObj(char* file_edge)
{
	ofstream out(file_edge);
	vector<int> flag;
	flag.push_back(0);

	for (int i = 0; i < m_edges.size(); ++i)
	{
		Edge edge = m_edges.at(i);
		vector<Site> ids = edge.focus;
		int count = 0;
		//for (int j = 0; j < ids.size(); ++j)
		//{
		//	if (ids.at(j).id == 2 || ids.at(j).id == 3 || ids.at(j).id == 4)
		//		count++;
		//}
		if (count < 2)
		{
			for (int j = 0; j < m_edges.at(i).line.size(); ++j)
			{
				out << "v" << " " << m_edges.at(i).line.at(j) << endl;
			}
			int end = flag.back() + m_edges.at(i).line.size();
			flag.push_back(end);
		}
		
	}
	for (int i = 0; i < flag.size() - 1; ++i)
	{
		int start = flag.at(i) + 1;
		int end = flag.at(i + 1);
		for (int j = start; j < end; j++)
			out << "l" << " " << j << " " << j + 1 << endl;
	}
	out.close();
}
void Apollonius3D::out_FaceObj(string file_face)
{
	ofstream out(file_face);
	int num = 0;
	if (m_faces.size() > 0)
	{
		for (int i = 0; i < edge_point_endid + 1; ++i)
		{
			out << "v " << m_points.at(i) << endl;
		}
		for (int i = 0; i < m_faces.size(); ++i)
		{
			Face face = m_faces.at(i);
			int id1 = face.site1;
			int id2 = face.site2;
			int count = 0;
			if (id1 >=0 && id2 >=0)
			{
				vector<Point> points = m_faces.at(i).Points;
				vector<tuple<int, int, int>> connect = m_faces.at(i).connect;
				for (int j = 0; j < points.size(); ++j)
				{
					out << "v " << points.at(j) << endl;
				}
				for (int j = 0; j < connect.size(); j++)
				{
					auto index0 = get<0>(connect.at(j));
					auto index1 = get<1>(connect.at(j));
					auto index2 = get<2>(connect.at(j));
					if (index0 > edge_point_endid)
						index0 = index0 + num;
					if (index1 > edge_point_endid)
						index1 = index1 + num;
					if (index2 > edge_point_endid)
						index2 = index2 + num;
					out << "f " << index0 + 1 << " " << index1 + 1 << " " << index2 + 1 << endl;
				}
				num = num + points.size();
			}
		}
	}
	
	out.close();
}
void Apollonius3D::out_faceone(string faceone, int id)
{
	ofstream out(faceone);
	int num = 0;
	for (int i = 0; i < edge_point_endid + 1; ++i)
	{
		out << "v " << m_points.at(i) << endl;
	}
	for (int i = 0; i < m_faces.size(); ++i)
	{
		Face face = m_faces.at(i);
		if (face.id == id)
		{
			vector<Point> points = m_faces.at(i).Points;
			vector<tuple<int, int, int>> connect = m_faces.at(i).connect;
			for (int j = 0; j < points.size(); ++j)
			{
				out << "v " << points.at(j) << endl;
			}
			for (int j = 0; j < connect.size(); j++)
			{
				auto index0 = get<0>(connect.at(j));
				auto index1 = get<1>(connect.at(j));
				auto index2 = get<2>(connect.at(j));
				out << "f " << index0 + 1 << " " << index1 + 1 << " " << index2 + 1 << endl;
			}
			break;
		}
	}
	out.close();
}
void Apollonius3D::out_edgeone(string edgeone, int id)
{
	ofstream out(edgeone);

	for (int i = 0; i < m_edges.size(); ++i)
	{
		Edge edge = m_edges.at(i);
		vector<Site_3D> ids = edge.focus;
		if (edge.id == id )
		{
			for (int j = 0; j < m_edges.at(i).line.size(); ++j)
			{
				out << "v" << " " << m_edges.at(i).line.at(j) << endl;
			}
			for (int j = 0; j < m_edges.at(i).line.size()-1; ++j)
			{
				out << "l " << j + 1 << " " << j + 2 << endl;
			}
		}	
	}
	out.close();
}
void Apollonius3D::out_PointObj(string file_face)
{
	ofstream out(file_face);
	for (int i = 0; i < m_points.size(); ++i)
	{
		out << "v " << m_points.at(i) << endl;
	}
	out.close();
}
void Apollonius3D::out_SitePointObj(string file_site)
{
	ofstream out(file_site);
	for (int i = 0; i < m_sites.size(); ++i)
	{
	//	if (m_sites.at(i).id >= 0)
		{
			out << "v " << m_sites.at(i).site_position << endl;
		}
	}
	out.close();
}
void Apollonius3D::out_SiteText(string file_site)
{
	ofstream out(file_site);
	for (int i = 0; i < m_sites.size(); ++i)
	{
		{
			out <<m_sites.at(i).site_position<<" "<< m_sites.at(i).weight<< endl;
		}
	}
	out.close();
}
void Apollonius3D::out_SiteObj(string file_site)
{

}
void Apollonius3D::out_VertexObj(string file_vertex)
{
	ofstream out(file_vertex);
	for (int i = 0; i < m_vertex.size(); ++i)
	{
		out << "v " << m_vertex.at(i).vertex_3d << endl;
	}
	out.close();
}
Site Apollonius3D::GetSite(int id)
{
	if (id < 0)
	{
		for (int i = m_sites.size()-1; i >= 0; i--)
		{
			if (id == m_sites.at(i).id)
				return m_sites.at(i);
		}
	}

	if (id == m_sites.at(id).id)
		return m_sites.at(id);
	else
	{
		int index1 = id - 1;
		int index2 = id + 1;
		while (index1 >= 0 || index2 < m_sites.size())
		{
			if (id == m_sites.at(index1).id)
				return m_sites.at(index1);
			if (id == m_sites.at(index2).id)
				return m_sites.at(index2);
			if (index1 > 0)
				index1--;
			if (index2 < m_sites.size() - 1)
				index2++;
		}
	}

}
Edge Apollonius3D::GetEdge(int id)
{
	if (id == m_edges.at(id).id)
		return m_edges.at(id);
	if (id < 0)
	{
		for (int i = m_edges.size(); i >= 0; --i)
		{
			if (id == m_edges.at(i).id)
				return m_edges.at(id);
		}
	}
	int index1 = id - 1;
	int index2 = id + 1;
	while (index1 >= 0 || index2 < m_edges.size())
	{
		if (id == m_edges.at(index1).id)
			return m_edges.at(index1);
		if (id == m_edges.at(index2).id)
			return m_edges.at(index2);
		if (index1 > 0)
			index1--;
		if (index2 < m_edges.size() - 1)
			index2++;
	}
}
Face Apollonius3D::GetFace(int site1, int site2)
{

}

void Apollonius3D::OrderTheFaceBoundary()
{
	for (int i = 0; i < m_faces.size(); ++i)
	{
		Face face = m_faces.at(i);
		vector<int> face_boundary = face.face_boundary;
		vector<vector<Point>> boundarys;
		while (face_boundary.size() > 0)
		{
		//	cout << face_boundary.size() << endl;
			vector<Point> boundary;
			boundary.push_back(GetEdge(face_boundary.at(0)).start);
			int size = face_boundary.size();
			bool flag = false;
			while (size > 0)
			{
				int  a = face_boundary.size();
				for (int j = 0; j < face_boundary.size(); ++j)
				{
					Point start = GetEdge(face_boundary.at(j)).start;
					Point end = GetEdge(face_boundary.at(j)).end;
					vector<Point> line = GetEdge(face_boundary.at(j)).line;
					if (start == boundary.at(boundary.size() - 1))
					{
						for (int k = 1; k < line.size(); ++k)
						{
							boundary.push_back(line.at(k));
						}
						face_boundary.erase(face_boundary.begin() + j);
						break;
					}
					if (end == boundary.at(boundary.size() - 1))
					{
						for (int k = line.size() - 2; k >= 0; --k)
						{
							boundary.push_back(line.at(k));
						}
						face_boundary.erase(face_boundary.begin() + j);
						break;
					}
					if (start == boundary.at(0))
					{
						for (int k = 1; k < line.size(); ++k)
						{
							boundary.insert(boundary.begin(), line.at(k));
						}
						face_boundary.erase(face_boundary.begin() + j);
						break;
					}
					if (end == boundary.at(0))
					{
						for (int k = line.size() - 2; k >= 0; --k)
						{
							boundary.insert(boundary.begin(), line.at(k));
						}
						face_boundary.erase(face_boundary.begin() + j);
						break;
					}
				}
				int b = face_boundary.size();
				if (a == b)
				{
					size = 0;
				}
			}
			Point start = boundary.at(0);
			Point end = boundary.at(boundary.size() - 1);
			Point normal(start.hx() - end.hx(), start.hy() - end.hy(), start.hz() - end.hz());
			double distance = Distance_3D(start, end);
			if (distance > step)
			{
				double distance_normal = sqrt(normal.hx() *normal.hx() + normal.hy() *normal.hy() + normal.hz() *normal.hz());
				double x_step = normal.hx()*(step / distance_normal);
				double y_step = normal.hy()*(step / distance_normal);
				double z_step = normal.hz()*(step / distance_normal);
				while (distance > 2 * step)
				{
					Point p(end.hx() + x_step, end.hy() + y_step, end.hz() + z_step);
					boundary.push_back(p);
					m_points.push_back(p);
					edge_point_endid++;
					end = p;
					distance = Distance_3D(start, end);
				}
			}
		//	if (size - face_boundary.size() >= 2)
			{
				boundarys.push_back(boundary);
			}
		}
		if (boundarys.size() == 0)
		{
			cout << "error" << endl;
			getchar();
		}
		if (boundarys.size() == 1)
		{
			face.boundary = boundarys.at(0);
		}
		if (boundarys.size() > 1) {
			map<int, double> polygon_order;
			for (int j = 0; j < boundarys.size(); ++j)
			{
				vector<Point> boundary = boundarys.at(j);
			
				polygon_order.insert(make_pair(j, boundary.size()));
			}
			vector<PAIR> _polygon_order(polygon_order.begin(), polygon_order.end());
			sort(_polygon_order.begin(), _polygon_order.end(), CmpByValue());
			face.boundary = boundarys.at(_polygon_order.at(1).first);
			face.in_boundary = boundarys.at(_polygon_order.at(0).first);
		}
	
		m_faces.at(i) = face;
	}


}
void Apollonius3D::ComputeFace()
{
//	num_virture_site = 0;
	cout <<"num_virture_site"<< num_virture_site << endl;
	for (int i = 0; i < m_sites.size() - num_virture_site; ++i)
	{
		int id1 = m_sites.at(i).id;
		map<int, vector<int>>::iterator edge_it1;
		edge_it1 = cellboundary.find(m_sites.at(i).id);
		vector<int> EdgeGroup1;
		if (edge_it1 != cellboundary.end())
		{
			EdgeGroup1 = cellboundary.at(m_sites.at(i).id);
		}
		double num;
		if (open == true)
			num = num_virture_site;
		else num = 0;
		for (int j = i + 1; j < m_sites.size()-num; ++j)
		{
			int id2 = m_sites.at(j).id;
			map<int, vector<int>>::iterator edge_it2;
			edge_it2 = cellboundary.find(m_sites.at(j).id);
			vector<int> EdgeGroup2;
			if (edge_it2 != cellboundary.end())
			{
				EdgeGroup2 = cellboundary.at(m_sites.at(j).id);
			}
			vector<int> faceboundary = GetInsectOfVector(EdgeGroup1, EdgeGroup2);
		
			bool flag = false;
			if (faceboundary.size() >= 2)
			{
				Edge edge1 = GetEdge(faceboundary.at(0));
				for (int k = 1; k < faceboundary.size(); ++k)
				{
					Edge edge2 = GetEdge(faceboundary.at(k));
					if (edge1.start == edge2.start || edge1.start == edge2.end || edge1.end == edge2.start || edge1.end == edge2.end)
					{
						flag = true;
						break;
					}
				}
			}
		//	cout << "flag" << flag << endl;
			if (flag == true)
			{
				Face face;
				face.site1 = id1;
				face.site2 = id2;
				face.face_boundary = faceboundary;
				face.id = m_faces.size();

				Point s1 = GetSite(id1).site_position;
				Point s2 = GetSite(id2).site_position;
			//	cout << "success" << endl;
				double x = s2.x() - s1.x();
				double y = s2.y() - s1.y();
				double z = s2.z() - s1.z();
				double m = sqrt(x*x + y*y + z*z);
				Point normal(x / m, y / m, z / m);
				face.normal = normal;
				m_faces.push_back(face);
			}
			
		}
	}

}

void Apollonius3D::MeshTheFace()
{
	OrderTheFaceBoundary();
    #pragma omp parallel for
	for (int i = 0; i < m_faces.size(); ++i)
	{
		AddThePointToOneFace(i);
	}
	for (int i = 0; i < m_faces.size(); ++i)
	{
		MapPointToOneFace(i);		
	}
}

bool Apollonius3D::IfIntheFace(int faceid, Point p)
{
	Face face = m_faces.at(faceid);
	Point site1 = GetSite(face.site1).site_position;
	Point site2 = GetSite(face.site2).site_position;
	double w1 = GetSite(face.site1).weight;
	double w2 = GetSite(face.site2).weight;

	double distance1 = Distance_3D(site1, p);
	double distance2 = Distance_3D(site2, p);
	double a = distance1 - w1;
	double b = distance2 - w2;
	if (abs(1) - abs(2) <= 1e6 || abs(1) - abs(2) >= -1e6)
		return true;
	else
		return false;
}
void Apollonius3D::MapPointToOneFace(int faceid)
{
	Face face = m_faces.at(faceid);
	double a = face.normal.x();
	double b = face.normal.y();
	double c = face.normal.z();
	double x1 = GetSite(face.site1).site_position.hx();
	double y1 = GetSite(face.site1).site_position.hy();
	double z1 = GetSite(face.site1).site_position.hz();

	double x2 = GetSite(face.site2).site_position.hx();
	double y2 = GetSite(face.site2).site_position.hy();
	double z2 = GetSite(face.site2).site_position.hz();

	double w1 = GetSite(face.site1).weight;
	double w2 = GetSite(face.site2).weight;

	for (int i = 0; i < face.Points.size(); ++i)
	{
		Point p = face.Points.at(i);

		double x0 = p.x();
		double y0 = p.y();
		double z0 = p.z();

		double aer1 = 2 * (a*(x0 - x1) + b*(y0 - y1) + c*(z0 - z1));
		double beta1 = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1);
		double aer2 = 2 * (a*(x0 - x2) + b*(y0 - y2) + c*(z0 - z2));
		double beta2 = (x0 - x2)*(x0 - x2) + (y0 - y2)*(y0 - y2) + (z0 - z2)*(z0 - z2);

		double w = (w1 - w2)*(w1 - w2);
		double A = (aer1 - aer2)*(aer1 - aer2) - 4 * w;
		double B = 2 * (beta1 + beta2 - w)* (aer1 + aer2) - 4 * (aer1*beta2 + aer2*beta1);
		double C = (beta1 + beta2 - w)*(beta1 + beta2 - w) - 4 * beta1*beta2;

		double t1 = (sqrt(B*B - 4 * A*C) - B) / (2 * A);
		double t2 = (-sqrt(B*B - 4 * A*C) - B) / (2 * A);

		double t;
		if (abs(t1) > abs(t2))
			t = t2;
		else t = t1;
		if (w1 == w2)
		{
			Point p1(x0, y0, z0);
			face.Points.at(i) = p1;
			m_points.push_back(p1);
		}
		else
		{
			Point p1(x0 + t* a, y0 + t*b, z0 + t*c);
			if (IfIntheFace(faceid, p1))//点在去曲面上
			{
				face.Points.at(i) = p1;
				m_points.push_back(p1);
			}
		}
	}
	for (int i = 0; i < face.boundary.size(); ++i)
	{
		Point p = face.boundary.at(i);

		double x0 = p.x();
		double y0 = p.y();
		double z0 = p.z();

		double aer1 = 2 * (a*(x0 - x1) + b*(y0 - y1) + c*(z0 - z1));
		double beta1 = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1);
		double aer2 = 2 * (a*(x0 - x2) + b*(y0 - y2) + c*(z0 - z2));
		double beta2 = (x0 - x2)*(x0 - x2) + (y0 - y2)*(y0 - y2) + (z0 - z2)*(z0 - z2);

		double w = (w1 - w2)*(w1 - w2);
		double A = (aer1 - aer2)*(aer1 - aer2) - 4 * w;
		double B = 2 * (beta1 + beta2 - w)* (aer1 + aer2) - 4 * (aer1*beta2 + aer2*beta1);
		double C = (beta1 + beta2 - w)*(beta1 + beta2 - w) - 4 * beta1*beta2;

		double t1 = (sqrt(B*B - 4 * A*C) - B) / (2 * A);
		double t2 = (-sqrt(B*B - 4 * A*C) - B) / (2 * A);

		double t;
		if (abs(t1) > abs(t2))
			t = t2;
		else t = t1;

		if (w1 == w2)
		{
			Point p1(x0, y0, z0);
			face.boundary.at(i) = p1;
			int index = face.boundary_index.at(i);
			if (index > -1)
				m_points.at(index) = p1;
		}
		else
		{
			Point p1(x0 + t* a, y0 + t*b, z0 + t*c);
			face.boundary.at(i) = p1;
			int index = face.boundary_index.at(i);
			if (index > -1)
				m_points.at(index) = p1;
		}
	}
	m_faces.at(faceid) = face;
}
//void Apollonius3D::AddThePointToOneFace(int k)
//{
//	string num = to_string(k);
//	Face face = m_faces.at(k);
//	vector<Point> boundary_3d = face.boundary;
//	vector<double> normal;
//	vector<double> v;
//	vector<double>  theta;
//	normal.push_back(face.normal.x());
//	normal.push_back(face.normal.y());
//	normal.push_back(face.normal.z());
//
//	RotateXYZ(normal, theta, boundary_3d);
//	Polygon2 boundary_poly;
//	double move_distance = boundary_3d.at(0).z();
//	vector<Point_2> boundary_2d;
//	for (int j = 0; j < boundary_3d.size(); ++j)
//	{
//		Point_2 p(boundary_3d.at(j).x(), boundary_3d.at(j).y());
//		boundary_poly.push_back(p);
//		boundary_2d.push_back(p);
//		int cur_index = getPointId(face.boundary.at(j));
//		face.boundary_index.push_back(cur_index);
//	}
//
//	int index_start = edge_point_endid+1;
//	double x_min_bound = boundary_poly.bbox().xmin();
//	double x_max_bound = boundary_poly.bbox().xmax();
//	double y_min_bound = boundary_poly.bbox().ymin();
//	double y_max_bound = boundary_poly.bbox().ymax();
//	double x_lenth = x_max_bound - x_min_bound;
//	double y_lenth = y_max_bound - y_min_bound;
//
//	int numofPoint = abs(boundary_poly.area()) / step / step;
//
//	vector<Point_2> m_ranpoints;
//	m_ranpoints = RandomPointInPolygon(boundary_poly, numofPoint);
//
//	CVT cvt(m_ranpoints, boundary_poly);
//	CVT_LBFGS cvt_lbfgs(cvt);
//	cvt_lbfgs.run();
//	vector<Point_2> cvt_site = cvt.m_new_site;
//
//	Point_2 p1(x_min_bound - x_lenth, y_min_bound - y_lenth);
//	Point_2 p2(x_min_bound - x_lenth, y_max_bound + y_lenth);
//	Point_2 p3(x_max_bound + x_lenth, y_max_bound + y_lenth);
//	Point_2 p4(x_max_bound + x_lenth, y_min_bound - y_lenth);
//	CDT  cdt;
//	for (int i = 0; i < boundary_2d.size(); ++i)
//	{
//		cdt.insert_constraint(boundary_2d.at(i), boundary_2d.at((i + 1) % boundary_2d.size()));
//	}
//	for (int i = 0; i < cvt_site.size(); ++i)
//	{
//		cdt.insert(cvt_site.at(i));
//	}
//	cdt.insert(p1);
//	cdt.insert(p2);
//	cdt.insert(p3);
//	cdt.insert(p4);
//	for (auto fit = cdt.faces_begin(); fit != cdt.faces_end(); fit++)
//	{
//		vector<int> indexs;
//		vector<Point_2>  delaunay_point;
//		for (int i = 0; i < 3; i++)
//		{
//			Point_2 p = fit->vertex(i)->point();
//			delaunay_point.push_back(p);
//			//find the id if in the bound;
//			auto it1 = std::find(boundary_2d.begin(),boundary_2d.end(),p);
//			if (it1 != boundary_2d.end())
//			{
//				auto id = it1 - boundary_2d.begin();
//				indexs.push_back(face.boundary_index.at(id));
//			}
//			//find the id if in the face
//			auto it2 = std::find(cvt_site.begin(), cvt_site.end(), p);
//			if (it2 != cvt_site.end())
//			{
//				auto id = index_start + it2 - cvt_site.begin();
//				indexs.push_back(id);
//			}
//		}
//		Polygon2  delaunay;
//		delaunay.push_back(delaunay_point.at(0));
//		delaunay.push_back(delaunay_point.at(1));
//		delaunay.push_back(delaunay_point.at(2));
//		Point_2 center = cvt.compute2DPolygonCentroid(delaunay);
//		if ((CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(boundary_poly.vertices_begin(), boundary_poly.vertices_end(), center, K())) && indexs.size() == 3)
//		{
//			tuple<int, int, int> indextuple(indexs.at(0), indexs.at(1), indexs.at(2));
//			if (indexs.at(0) > -1 && indexs.at(1) > -1 && indexs.at(2) > -1)
//				face.connect.push_back(indextuple);
//		}
//
//	}
//	vector<Point> cvt_site_3d;
//	for (int i = 0; i < cvt_site.size(); ++i)
//	{
//		Point p_3(cvt_site.at(i).hx(), cvt_site.at(i).hy(), move_distance);
//		cvt_site_3d.push_back(p_3);
//	}
//	ReRotateXYZ(normal, theta, cvt_site_3d);
//	for (int i = 0; i < cvt_site_3d.size(); ++i)
//	{
//		face.Points.push_back(cvt_site_3d.at(i));
//	}
//	m_faces.at(k) = face;
//}

void Apollonius3D::AddThePointToOneFace(int k)
{
	string num = to_string(k);
	Face face = m_faces.at(k);
	vector<Point> boundary_3d = face.boundary;
	vector<Point> in_boundary_3d = face.in_boundary;
	vector<double> normal;
	vector<double> v;
	vector<double>  theta;
	normal.push_back(face.normal.x());
	normal.push_back(face.normal.y());
	normal.push_back(face.normal.z());
	int index_start = edge_point_endid + 1;
	RotateXYZ(normal, theta, boundary_3d);
	if(in_boundary_3d.size()>0)
	RotateXYZ(normal, theta, in_boundary_3d);
	Polygon2 boundary_poly;
	Polygon2 in_boundary_poly;
	double move_distance = boundary_3d.at(0).z();
	vector<Point_2> boundary_2d;
	vector<Point_2> in_boundary_2d;
	for (int j = 0; j < boundary_3d.size(); ++j)
	{
		Point_2 p(boundary_3d.at(j).x(), boundary_3d.at(j).y());
		boundary_poly.push_back(p);
		boundary_2d.push_back(p);
		int cur_index = getPointId(face.boundary.at(j));
		face.boundary_index.push_back(cur_index);
	}
	for (int j = 0; j < in_boundary_3d.size(); ++j)
	{
		Point_2 p(in_boundary_3d.at(j).x(), in_boundary_3d.at(j).y());
		in_boundary_poly.push_back(p);
		in_boundary_2d.push_back(p);
		int cur_index = getPointId(face.in_boundary.at(j));
		face.in_boundary_index.push_back(cur_index);
	}
	double x_min_bound = boundary_poly.bbox().xmin();
	double x_max_bound = boundary_poly.bbox().xmax();
	double y_min_bound = boundary_poly.bbox().ymin();
	double y_max_bound = boundary_poly.bbox().ymax();
	double x_lenth = x_max_bound - x_min_bound;
	double y_lenth = y_max_bound - y_min_bound;

	int numofPoint = abs(boundary_poly.area()) / step / step;

	vector<Point_2> m_ranpoints;
	m_ranpoints = RandomPointInPolygon(boundary_poly, numofPoint);

	CVT cvt(m_ranpoints, boundary_poly);
	CVT_LBFGS cvt_lbfgs(cvt);
	cvt_lbfgs.run();
//	cvt.draw_CVT_to_obj("E://code//Apollonius//output//cvt" + num + ".obj");
//	cvt.draw_site_to_obj("E://code//Apollonius//output//site" + num + ".obj");
	vector<Point_2> cvt_site = cvt.m_new_site;
	if (in_boundary_2d.size() != 0)
	{
		for (int j = 0; j < cvt_site.size(); ++j)
		{
			Point_2 p = cvt_site.at(j);
			if ((CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(in_boundary_poly.vertices_begin(), in_boundary_poly.vertices_end(), p, K())))
			{
				cvt_site.erase(cvt_site.begin()+j);
				j--;
			}
		}
	}
	
	Point_2 p1(x_min_bound - x_lenth, y_min_bound - y_lenth);
	Point_2 p2(x_min_bound - x_lenth, y_max_bound + y_lenth);
	Point_2 p3(x_max_bound + x_lenth, y_max_bound + y_lenth);
	Point_2 p4(x_max_bound + x_lenth, y_min_bound - y_lenth);
	CDT  cdt;
	for (int i = 0; i < boundary_2d.size(); ++i)
	{
		cdt.insert_constraint(boundary_2d.at(i), boundary_2d.at((i + 1) % boundary_2d.size()));
	}
	for (int i = 0; i < in_boundary_2d.size(); ++i)
	{
		cdt.insert_constraint(in_boundary_2d.at(i), in_boundary_2d.at((i + 1) % in_boundary_2d.size()));
	}
	for (int i = 0; i < cvt_site.size(); ++i)
	{
		cdt.insert(cvt_site.at(i));
	}
	cdt.insert(p1);
	cdt.insert(p2);
	cdt.insert(p3);
	cdt.insert(p4);
	for (auto fit = cdt.faces_begin(); fit != cdt.faces_end(); fit++)
	{
		vector<int> indexs;
		vector<Point_2>  delaunay_point;
		for (int i = 0; i < 3; i++)
		{
			Point_2 p = fit->vertex(i)->point();
			delaunay_point.push_back(p);
			//find the id if in the bound;
			auto it1 = std::find(boundary_2d.begin(), boundary_2d.end(), p);
			if (it1 != boundary_2d.end())
			{
				auto id = it1 - boundary_2d.begin();
				indexs.push_back(face.boundary_index.at(id));
			}
			auto it3 = std::find(in_boundary_2d.begin(), in_boundary_2d.end(), p);
			if (it3 != in_boundary_2d.end())
			{
				auto id = it3 - in_boundary_2d.begin();
				indexs.push_back(face.in_boundary_index.at(id));
			}
			//find the id if in the face
			auto it2 = std::find(cvt_site.begin(), cvt_site.end(), p);
			if (it2 != cvt_site.end())
			{
				auto id = index_start + it2 - cvt_site.begin();
				indexs.push_back(id);
			}
		}
		Polygon2  delaunay;
		delaunay.push_back(delaunay_point.at(0));
		delaunay.push_back(delaunay_point.at(1));
		delaunay.push_back(delaunay_point.at(2));
		Point_2 center = cvt.compute2DPolygonCentroid(delaunay);
		if ((CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(boundary_poly.vertices_begin(), boundary_poly.vertices_end(), center, K()))
			&& indexs.size() == 3
			&& (CGAL::ON_UNBOUNDED_SIDE == CGAL::bounded_side_2(in_boundary_poly.vertices_begin(), in_boundary_poly.vertices_end(), center, K())))
		{
			tuple<int, int, int> indextuple(indexs.at(0), indexs.at(1), indexs.at(2));
			if (indexs.at(0) > -1 && indexs.at(1) > -1 && indexs.at(2) > -1)
				face.connect.push_back(indextuple);
		}

	}
	vector<Point> cvt_site_3d;
	for (int i = 0; i < cvt_site.size(); ++i)
	{
		Point p_3(cvt_site.at(i).hx(), cvt_site.at(i).hy(), move_distance);
		cvt_site_3d.push_back(p_3);
	}
	ReRotateXYZ(normal, theta, cvt_site_3d);
	for (int i = 0; i < cvt_site_3d.size(); ++i)
	{
		face.Points.push_back(cvt_site_3d.at(i));
	}
	m_faces.at(k) = face;
}
int Apollonius3D::getPointId(Point p)
{
	int index;
	auto it = std::find(m_points.begin(),m_points.end(),p);
	if (it != m_points.end())
	{
		index = it - m_points.begin();
	}
	else
	{
		index = -1;
	}
	return index;
}
Face Apollonius3D::GetFace(int id)
{
	Face face;
	if (m_faces.at(id).id == id)
		return m_faces.at(id);
	else
		for (int i = 0; i < m_faces.size(); ++i)
		{
			if (m_faces.at(i).id == id)
				return m_faces.at(i);
		}
	face.id = -1;
	return face;
}
void Apollonius3D::ComputerCell()
{
	for (int i = 0; i < m_sites.size(); ++i)
	{
		Site site = m_sites.at(i);
		if (site.id >= 0)
		{
			Cell_3D cell;
			int id = site.id;
			vector<int> faceids;
			vector<int> siteids;
			for (int j = 0; j < m_faces.size(); ++j)
			{
				Face face = m_faces.at(j);
				int id1 = face.site1;
				int id2 = face.site2;
				if (id == id1)
				{
					faceids.push_back(face.id);
					siteids.push_back(id2);
				}
				if (id == id2)
				{
					faceids.push_back(face.id);
					siteids.push_back(id1);
				}
			}
			cell.site_id = id;
			cell.face_ids = faceids;
			cell.site_ids = siteids;
			m_cells.push_back(cell);
		}
	}
	for (int i = 0; i < m_cells.size(); ++i)
	{
		Cell_3D cell = m_cells.at(i);
		Point p = GetSite(cell.site_id).site_position;
		double area = 0;
		double x = 0;
		double y = 0;
		double z = 0;
		for (int j = 0; j < cell.face_ids.size(); ++j)
		{
			int  faceid = cell.face_ids.at(j);
			Face face = GetFace(faceid);
			if (face.id == -1)
				break;
			vector<Point> points;
			for (int k = 0; k < edge_point_endid + 1; ++k)
			{
				points.push_back(m_points.at(k));
			}
			for (int k = 0; k < face.Points.size(); ++k)
			{
				points.push_back(face.Points.at(k));
			}
			vector<tuple<int, int, int>> connect = face.connect;
			for (int k = 0; k < connect.size(); k++)
			{
				auto index0 = get<0>(connect.at(k));
				auto index1 = get<1>(connect.at(k));
				auto index2 = get<2>(connect.at(k));
				auto p0 = points.at(index0);
				auto p1 = points.at(index1);
				auto p2 = points.at(index2);
				Eigen::Matrix4d m;
				m << 1, 1, 1, 1,
					p.hx(), p0.hx(), p1.hx(), p2.hx(),
					p.hy(), p0.hy(), p1.hy(), p2.hy(),
					p.hz(), p0.hz(), p1.hz(), p2.hz();
				double are = abs(m.determinant()) / 6;
				x = x + (p0.hx() + p1.hx() + p2.hx() + p.hx())*are / 4;
				y = y + (p0.hy() + p1.hy() + p2.hy() + p.hy())*are / 4;
				z = z + (p0.hz() + p1.hz() + p2.hz() + p.hz())*are / 4;
				area = area + are;
			}
		}
		x = x / area;
		y = y / area;
		z = z / area;
		Point center(x,y,z);
		cell.center = center;
		cell.id = i;
		m_cells.at(i) = cell;
	}
}
void Apollonius3D::out_cell_one(string file_cell, int id)
{
	ofstream out(file_cell);
	for (int i = 0; i < m_cells.size(); ++i)
	{
		Cell_3D cell = m_cells.at(i);
		if (cell.id == id)
		{
			for (int j = 0; j < edge_point_endid + 1; ++j)
			{
				out << "v " << m_points.at(j) << endl;
			}
			int num = 0;
			for (int j = 0; j < cell.face_ids.size(); ++j)
			{
				int face_id = cell.face_ids.at(j);
				Face face = GetFace(face_id);
				vector<Point> points = face.Points;
				vector<tuple<int, int, int>> connect = face.connect;
				for (int k = 0; k < points.size(); ++k)
				{
					out << "v " << points.at(k) << endl;
				}
				for (int k = 0; k < connect.size(); k++)
				{
					auto index0 = get<0>(connect.at(k));
					auto index1 = get<1>(connect.at(k));
					auto index2 = get<2>(connect.at(k));
					if (index0 > edge_point_endid)
						index0 = index0 + num;
					if (index1 > edge_point_endid)
						index1 = index1 + num;
					if (index2 > edge_point_endid)
						index2 = index2 + num;
					out << "f " << index0 + 1 << " " << index1 + 1 << " " << index2 + 1 << endl;
				}
				num = num + points.size();
			}
		}
	}
	out.close();
}



