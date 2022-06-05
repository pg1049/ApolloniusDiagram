#include "Edge.h"
#include "Site.h"
#include "solve_equation.h"

typedef Edge_3DApollonius::Edge_3D Edge;
ostream& operator<<(ostream& out, Edge edge)
{
	out << "Id:" << edge.id << "     focus:    ";
	for (int i = 0; i < edge.focus.size(); ++i)
	{
		out << "(" << edge.focus.at(i).site_position << ") ";
	}
	out << endl;
	//	out << "(" << edge.start << ")    (" << edge.end << ") ";
	out << "(" << edge.start_id << ")    (" << edge.end_id << ") " << endl;
	out << edge.line.size() << endl;
	out << edge.normal << endl;
}
bool operator==( Edge edge1,  Edge edge2)
{
	if (edge1.id == edge2.id)
		return true;
	else return false;
}
void Edge_3DApollonius::insert_poinBetweenEndpoint(Edge& hyperbola1, double step)
{
	Point startpoint = hyperbola1.start;
	Point endpoint = hyperbola1.end;
	double distance = sqrt(((startpoint.hx() - endpoint.hx()) * (startpoint.hx() - endpoint.hx())) +
		((startpoint.hy() - endpoint.hy()) * (startpoint.hy() - endpoint.hy())) +
		((startpoint.hz() - endpoint.hz()) * (startpoint.hz() - endpoint.hz())));
	hyperbola1.line.push_back(startpoint);
	Point _nor = hyperbola1.normal;
	bool flag = true;
	while (distance > 1* step)
	{
	//	cout << distance<<" "<< step << endl;
		Point normal1 = get_Fun_face_normal(hyperbola1.focus.at(0).site_position, hyperbola1.focus.at(1).site_position, startpoint);
		Point normal2 = get_Fun_face_normal(hyperbola1.focus.at(0).site_position, hyperbola1.focus.at(2).site_position, startpoint);
		double x1 = normal1.hx();
		double y1 = normal1.hy();
		double z1 = normal1.hz();
		double x2 = normal2.hx();
		double y2 = normal2.hy();
		double z2 = normal2.hz();
		Point normal((y1*z2 - y2*z1), -(x1*z2 - x2*z1), (x1*y2 - x2*y1));
		double result = dot_product(normal, _nor);
		if (result < 0)
		{
			Point renormal(-normal.hx(), -normal.hy(), -normal.hz());
			normal = renormal;
		}
		_nor = normal;
		double distance_normal = sqrt(normal.hx() *normal.hx() + normal.hy() *normal.hy() + normal.hz() *normal.hz());
		double x = startpoint.hx() + normal.hx()*(step / distance_normal);
		double y = startpoint.hy() + normal.hy()*(step / distance_normal);
		double z = startpoint.hz() + normal.hz()*(step / distance_normal);
		Point p(x, y, z);
		vector<Site> m_site;
		m_site = hyperbola1.focus;

		vector<double> P_result = solve_equation(m_site, p);
	
		if (P_result.size() > 0)
		{
			Point p_result(P_result.at(0), P_result.at(1), P_result.at(2));
			hyperbola1.line.push_back(p_result);
			startpoint = p_result;
		}
		else
		{
			startpoint = p;
		}
		distance = sqrt(((startpoint.hx() - endpoint.hx()) * (startpoint.hx() - endpoint.hx())) +
			((startpoint.hy() - endpoint.hy()) * (startpoint.hy() - endpoint.hy())) +
			((startpoint.hz() - endpoint.hz()) * (startpoint.hz() - endpoint.hz())));
	double dis1 = Distance_3D(startpoint, hyperbola1.focus.at(0).site_position) - hyperbola1.focus.at(0).weight;
	double dis2 = Distance_3D(startpoint, hyperbola1.site.site_position) - hyperbola1.site.weight;
	if (dis1 > dis2)
	{
		flag = false;
		hyperbola1.line.clear();
		break;
	}
	}
	if (flag == false)
	{
		_nor = hyperbola1.normal;
		startpoint = hyperbola1.start;
		hyperbola1.line.push_back(startpoint);
		Point _nor2(-_nor.hx(), -_nor.hy(), -_nor.hz());
		int count = 0;
		while (distance > 1 * step)
		{
			cout << distance << endl;
			Point normal1 = get_Fun_face_normal(hyperbola1.focus.at(0).site_position, hyperbola1.focus.at(1).site_position, startpoint);
			Point normal2 = get_Fun_face_normal(hyperbola1.focus.at(0).site_position, hyperbola1.focus.at(2).site_position, startpoint);
			double x1 = normal1.hx();
			double y1 = normal1.hy();
			double z1 = normal1.hz();
			double x2 = normal2.hx();
			double y2 = normal2.hy();
			double z2 = normal2.hz();
			Point normal((y1*z2 - y2*z1), -(x1*z2 - x2*z1), (x1*y2 - x2*y1));
			double result = dot_product(normal, _nor2);
			if (result < 0)
			{
				Point renormal(-normal.hx(), -normal.hy(), -normal.hz());
				normal = renormal;
			}
			_nor2 = normal;
			double distance_normal = sqrt(normal.hx() *normal.hx() + normal.hy() *normal.hy() + normal.hz() *normal.hz());
			double x = startpoint.hx() + normal.hx()*(step / distance_normal);
			double y = startpoint.hy() + normal.hy()*(step / distance_normal);
			double z = startpoint.hz() + normal.hz()*(step / distance_normal);
			Point p(x, y, z);
			vector<Site> m_site;
			m_site = hyperbola1.focus;

			vector<double> P_result = solve_equation(m_site, p);

			if (P_result.size() > 0)
			{
				Point p_result(P_result.at(0), P_result.at(1), P_result.at(2));
				hyperbola1.line.push_back(p_result);
				startpoint = p_result;
			}
			else
			{
				startpoint = p;
			}
			distance = sqrt(((startpoint.hx() - endpoint.hx()) * (startpoint.hx() - endpoint.hx())) +
				((startpoint.hy() - endpoint.hy()) * (startpoint.hy() - endpoint.hy())) +
				((startpoint.hz() - endpoint.hz()) * (startpoint.hz() - endpoint.hz())));
			count++;
			if (count > 100)
				distance = 0;
		}
	}
	Point end;
	end = hyperbola1.line.at(hyperbola1.line.size() - 1);
	double dis = Distance_3D(end, endpoint);
	if (dis < 0.5*step&&hyperbola1.line.size()>1)
	{
		hyperbola1.line.erase(hyperbola1.line.begin() + hyperbola1.line.size() - 1);
	}
//	cout << "success" << endl;
	hyperbola1.line.push_back(endpoint);
}
Edge Edge_3DApollonius::creatEdge(vector<Site> sg, vector<Vertex> vg)
{
	Edge edge;
	Vertex vertex1 = vg.at(0);
	Vertex vertex2 = vg.at(1);
	edge.id = m_edges.size();
	edge.focus = sg;
	edge.start = vertex1.vertex_3d;
	edge.end = vertex2.vertex_3d;
	edge.start_id = vertex1.id;
	edge.end_id = vertex2.id;
	vector<Site>  all_site = vertex1.site_group;
	int id;
	for (int i = 0; i < all_site.size(); ++i)
	{
		int id1 = all_site.at(i).id;
		bool flag = true;
		for (int j = 0; j < sg.size(); ++j)
		{
			int id2 = sg.at(j).id;
			if (id1 == id2)
			{
				flag = false;
				break;
			}
		}
		if (flag == true)
		{
			id = i;
			break;
		}
	}
	Site S = all_site.at(id);
	Point normal1(edge.end.hx() - edge.start.hx(), edge.end.hy() - edge.start.hy(), edge.end.hz() - edge.start.hz());
	edge.normal = normal1;
	edge.site = S;
	return edge;

}
void Edge_3DApollonius::GetEdgeEndPoint()
{
	vector<VertexAndSite> m_vs1; //存放一组site生成两个vertex的情况
	vector<VertexAndSite> m_vs; //存放普通的情况
	vector<VertexAndSite> m_vs2;
	vector<Bound_vertex> bound_vertex;

	for (int i = 0; i < m_vertex.size(); ++i)
	{
		Vertex vertex1 = m_vertex.at(i);
		vector<Site> sitegroup1 = vertex1.site_group;
		vector<vector<Site>> S;
		bool flag = 0;
		for (int j = i + 1; j < m_vertex.size(); ++j)
		{
			Vertex vertex2 = m_vertex.at(j);
			vector<Site> sitegroup2 = vertex2.site_group;
			vector<Site>  sitegroup;
			GetIntersectOfSitegroup(sitegroup1, sitegroup2, sitegroup);
			bool flag = false;
			int count = 0;
			for (int k = 0; k < sitegroup.size(); ++k)
			{
				Site s = sitegroup.at(k);
				if (s.id < 0)
					count++;
			}
			if (count >= 3)
				flag = true;
			if (flag == true)
			{
				sitegroup.clear();
			}
			if (sitegroup.size() == sitegroup1.size() && sitegroup.size() == sitegroup2.size())
			{
				VertexAndSite vs;
				vs.sg = sitegroup;
				vs.vg.push_back(vertex1);
				vs.vg.push_back(vertex2);
				m_vs1.push_back(vs);
			}
			else
				if (sitegroup.size() >= 3)
				{
					VertexAndSite vs;
					vs.sg = sitegroup;
					vs.vg.push_back(vertex1);
					vs.vg.push_back(vertex2);
					m_vs.push_back(vs);
				}
		}
	}
	//for (int i = 0; i < m_vs.size(); ++i)
	//{
	//		vector<Site> sg1 = m_vs.at(i).sg;
	//		vector<Vertex> vg1 = m_vs.at(i).vg;
	//		for (int j = i + 1; j < m_vs.size(); ++j)
	//		{
	//				vector<Site> sg2 = m_vs.at(j).sg;
	//				vector<Vertex> vg2 = m_vs.at(j).vg;
	//				if (sg1 == sg2)
	//				{
	//					vg1.insert(vg1.end(), vg2.begin(), vg2.end());
	//					m_vs.erase(m_vs.begin() + j);
	//					--j;
	//				}
	//		}
	//		for (int j = 0; j < vg1.size(); ++j)
	//		{
	//			Vertex v1 = vg1.at(j);
	//			for (int k = j + 1; k < vg1.size(); ++k)
	//			{
	//				Vertex v2 = vg1.at(k);
	//				if (v1.id == v2.id)
	//				{
	//					vg1.erase(vg1.begin() + k);
	//					k--;
	//				}
	//			}
	//		}
	//		m_vs.at(i).vg = vg1;
	//}

	for (int i = 0; i < m_vs1.size(); ++i)
	{
		vector<Site_3D> sg = m_vs1.at(i).sg;
		vector<Vertex> vg = m_vs1.at(i).vg;
		double min_w = sg.at(0).weight;
		int id = 0;
		for (int j = 0; j < sg.size(); ++j)
		{
			if (sg.at(j).weight < min_w)
			{
				min_w = sg.at(j).weight;
				id = j;
			}
		}
		for(int j =0;j<sg.size()-2;++j)
			for (int m = j+1; m<sg.size()-1; ++m)
				for (int n = m+1; n < sg.size(); ++n)
				{
					vector<Site_3D> _sg = {sg.at(j),sg.at(m),sg.at(n)};
					if (j == id || m == id || n == id)
					{
						Edge edge = creatEdge(_sg, vg);
						m_edges.push_back(edge);
					}
				}
	}
	cout << "success" << endl;
	for (int i = 0; i < m_vs.size(); ++i)
	{
		vector<Site> sitegroup = m_vs.at(i).sg;
		vector<Vertex> vg = m_vs.at(i).vg;
		if (vg.size() == 2 && sitegroup.size() >= 3)//最普通的情况
		{
			Edge edge = creatEdge(sitegroup, vg);
			m_edges.push_back(edge);
		}	
		if (vg.size() > 2)
		{
			Point site1 = sitegroup.at(0).site_position;
			Point site2 = sitegroup.at(1).site_position;
			Point site3 = sitegroup.at(2).site_position;
			double x1 = site1.hx();
			double y1 = site1.hy();
			double z1 = site1.hz();
			double x2 = site2.hx();
			double y2 = site2.hy();
			double z2 = site2.hz();
			double x3 = site3.hx();
			double y3 = site3.hy();
			double z3 = site3.hz();
			Point center((x1 + x2 +x3)/3, (y1 + y2 +y3) / 3, (z1+z2+z3) / 3);
			double cx = center.hx();
			double cy = center.hy();
			double cz = center.hz();
			vector<double> v1 = {x1 - cx, y1 - cy, z1- cz};
			vector<double> v2 = { x2 - cx, y2 - cy, z2 - cz };
			vector<double> v = {v1.at(1)*v2.at(2)- v1.at(2)*v2.at(1),-(v1.at(0)*v2.at(2) - v1.at(2)*v2.at(0)),v1.at(0)*v2.at(1) - v1.at(1)*v2.at(0) };
			v = { v.at(0) + cx, v.at(1) + cy, v.at(2) + cz };
			map<int,double> vg_order;
			for (int j = 0; j < vg.size(); j++)
			{
				vector<double> _v = {vg.at(j).vertex_3d.hx() - cx, vg.at(j).vertex_3d.hy() - cy , vg.at(j).vertex_3d.hz() - cz };
				double lenth = Dot(v,_v);
				vg_order.insert(make_pair(j,lenth));
			}
			vector<PAIR> _vg_order(vg_order.begin(), vg_order.end());
			sort(_vg_order.begin(), _vg_order.end(), CmpByValue());
		   for(int j =0;j<_vg_order.size()-1;++j)
		   {
			   int id1 = _vg_order.at(j).first;
			   int id2 = _vg_order.at(j+1).first;
				Vertex v1 = vg.at(id1);
				Vertex v2 = vg.at(id2);
				if (!(v1.site_group == v2.site_group))
				{
					vector<Vertex> _vg = {v1 ,v2};
					Edge edge = creatEdge(sitegroup, _vg);
					m_edges.push_back(edge);
				}
			}
		}
	}
	cout << "success" << endl;
}

void Edge_3DApollonius::get_hyperbola(double step)
{
	clock_t starttime1, endtime1;
	starttime1 = clock();
	GetEdgeEndPoint();
	endtime1 = clock();
//	cout << "edge time:" << (double)(endtime1 - starttime1) / CLOCKS_PER_SEC << "s" << endl;
  #pragma omp parallel for
	for (int i = 0; i < m_edges.size(); ++i)
	{
		Edge edge = m_edges.at(i);
		insert_poinBetweenEndpoint(edge, step);
		m_edges.at(i) = edge;
	}
	
	for (int i = 0; i < m_edges.size(); ++i)
	{
		Edge edge = m_edges.at(i);
		if (edge.line.size() >= 2)
		{
			edge.id = i;
			m_edges.at(i) = edge;
		}
		else
		{
			m_edges.erase(m_edges.begin()+i);
			i--;
		}
	}
}
void Edge_3DApollonius::remove_pointofhyperbola(Edge&  edge, vector<double>& boundary)
{
	double x_start = boundary.at(0);
	double x_end = boundary.at(1);
	double y_start = boundary.at(2);
	double y_end = boundary.at(3);
	double z_start = boundary.at(4);
	double z_end = boundary.at(5);
	for (int j = 0; j < edge.line.size(); ++j)
	{
		double x = edge.line.at(j).hx();
		double y = edge.line.at(j).hy();
		double z = edge.line.at(j).hz();
		if (x<x_start || x>x_end || y<y_start || y>y_end || z<z_start || z>z_end)
		{
			edge.line.erase(edge.line.begin() + j);
			--j;
		}
	}
	if (edge.line.size() >= 2)
	{
		edge.start = edge.line.at(0);
		edge.end = edge.line.at(edge.line.size() - 1);
	}

}

map<int, vector<int>> Edge_3DApollonius::GetCellboundary()
{
	map<int, vector<int>> m_cell;
	for (int i = 0; i < m_sites.size(); ++i)
	{
		Site site = m_sites.at(i);
		int id = site.id;
		vector<int> cell_boundary;
		for (int j = 0; j < m_edges.size(); ++j)
		{
			Edge edge = m_edges.at(j);
			vector<Site> sites = edge.focus;
			for (int k = 0; k < sites.size(); ++k)
			{
				if (site == sites.at(k))
				{
					cell_boundary.push_back(edge.id);
					break;
				}
			}
		}
		if (cell_boundary.size() > 1)
		{
			m_cell[id] = cell_boundary;
		}
	}
	return m_cell;
}