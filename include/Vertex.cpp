#include "Vertex.h"
#include "solve_equation.h"
#include <time.h>
typedef Vertex_3DApollonius::range_vertex Range_vertex;
typedef Vertex_3DApollonius::Vertex_3D Vertex;
typedef Vertex_3DApollonius::VertexAndSite VertexAndSite;
ostream& operator<<(ostream& out, Vertex& m_vertex)
{
	out << "vertex: " << "(" << m_vertex.vertex_3d << ")";
	out << "      site group:";
	for (int j = 0; j < m_vertex.site_group.size(); ++j)
	{
		out << "(" << m_vertex.site_group.at(j).site_position << ")";
	}
	out << endl;
	return out;
}
bool operator==(const Vertex  V1, const Vertex V2)
{
	if (V1.site_group.size() != V2.site_group.size())
		return false;
	else
	{
		int flag = 0;
		int size = V1.site_group.size();
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				if (V1.site_group.at(i).id == V1.site_group.at(j).id)
				{
					flag++;
					break;
				}
			}
		if ((flag == size) && (V1.vertex_3d == V2.vertex_3d))
			return true;
		else
			return false;
	}
}
void Vertex_3DApollonius::vertex_3D_get()
{
	clock_t starttime, endtime;
	starttime = clock();
	Generate_Cube(m_cube_boundary, m_sites.size());
	first_update_cube();
	int  count = 0;
	clock_t starttime3, endtime3;
	starttime3 = clock();
	update_cube();
	endtime3 = clock();
	std::cout << "update_cube:" << (double)(endtime3 - starttime3) / CLOCKS_PER_SEC << "s" << endl;

	clock_t starttime4, endtime4;
	starttime4 = clock();
	remove_cube();
//	compute_vertex();
	endtime4 = clock();
	std::cout << "remove_cube:" << (double)(endtime4 - starttime4) / CLOCKS_PER_SEC << "s" << endl;
	double lenth_flag = DBL_MAX;

	clock_t starttime5, endtime5;
	starttime5= clock();
	while (lenth_flag > 0.0001)
	{
		refine_cube(lenth_flag);
		remove_cube();
	//	cout << m_cube.size() << endl;
	}
	endtime5 = clock();
	std::cout << "remove_cube:" << (double)(endtime5 - starttime5) / CLOCKS_PER_SEC << "s" << endl;

	for (int i = 0; i < m_cube.size(); ++i)
	{
		m_cube_result.push_back(m_cube.at(i));
	}
	m_cube.clear();
	confirm_vertex();
	endtime = clock();
	std::cout << "get vertex time:" << (double)(endtime - starttime) / CLOCKS_PER_SEC << "s" << endl;

}
Cube Vertex_3DApollonius::update_cube_one(Site site_one, Cube cube_one, bool& flag)
{
	for (int i = 0; i < cube_one.sites.size(); ++i)
	{
		auto site = cube_one.sites.at(i).first;
		if (site == site_one)
		{
			flag = false;
			return cube_one;
		}
	}
	double lenth = sqrt(3) / 2 * cube_one.lenthOfcube_sides;
	double distance = Distance_3D(site_one.site_position, cube_one.pointOfcube_cente) -site_one.weight;
	if (distance < cube_one.min_distance)
	{
		cube_one.min_distance = distance;
		vector<pair<Site_3D, double>> site_dis;
		site_dis.push_back(make_pair(site_one, distance));
		for (int i = 0; i < cube_one.sites.size(); ++i)
		{
			double dis = cube_one.sites.at(i).second;
			if (dis - cube_one.min_distance <= sqrt(3)*cube_one.lenthOfcube_sides)
			{
				site_dis.push_back(cube_one.sites.at(i));
			}
		}
		cube_one.sites.clear();
		cube_one.sites = site_dis;
		flag = true;
	}
	else
		if (distance > cube_one.min_distance)
		{
			if (distance - cube_one.min_distance <= sqrt(3)*cube_one.lenthOfcube_sides)
			{
				cube_one.sites.push_back(make_pair(site_one, distance));
				flag = true;
			}
		}
		else
			flag = false;
	return cube_one;
}
void Vertex_3DApollonius::first_update_cube()
{
	for (int i = 0; i < m_sites.size(); ++i)
	{
		Point p_position = get_SiteInCube(m_cube.at(0), m_sites.at(i).site_position);//判断点的坐标判断所在格子的位置
		Cube cube_one = get_CubeByIndex(m_cube, p_position);//根据格子的位置得到这个格子信息
		bool flag;
		m_cube.at(cube_one.id) = update_cube_one(m_sites.at(i), cube_one, flag);//更新一个格子；
		m_cube.at(cube_one.id).empty = 0;
	}
}
void Vertex_3DApollonius::update_cube()
{
	priority_queue<Event> m_prio_event;
#pragma omp parallel for
	for (int i = 0; i < m_cube.size(); ++i)
	{
		Event m_event;
		m_event.cube_one = m_cube.at(i);
		double range_min = m_cube.at(i).min_distance;
		m_event.distance = range_min;
		if (m_event.distance < DBL_MAX)
		{
            #pragma omp critical
			m_prio_event.push(m_event);
		}	
	}
	while (m_prio_event.size() > 0 )
	{
		Event m_event;
		m_event = m_prio_event.top();
		m_prio_event.pop();
		vector<Point> cube_neighbor = get_CubeNeighbor(m_event.cube_one);
#pragma omp parallel for
		for (int j = 0; j < cube_neighbor.size(); ++j)
		{
			Cube cube_one = get_CubeByIndex(m_cube, cube_neighbor.at(j));
			bool flag2 = false;
#pragma omp parallel for
			for (int k = 0; k < m_event.cube_one.sites.size(); ++k)
			{
				bool flag;
				cube_one = update_cube_one(m_event.cube_one.sites.at(k).first, cube_one, flag);
				if (flag == true)
					flag2 = true;
			}
			if (flag2 == true)
			{
				m_cube.at(cube_one.id) = cube_one;
				Event m_event1;
				m_event1.cube_one = cube_one;
				double range_min1 = cube_one.min_distance;
				m_event1.distance = range_min1;
               #pragma omp critical
				m_prio_event.push(m_event1);
			}
		}
		
	}
}

void Vertex_3DApollonius::compute_vertex()
{
	Cube box;
	vector<double> result;
	for (int i = 0; i < m_cube.size(); ++i)
	{	
	    box = m_cube.at(i);
		if (box.sites.size() == 4)
		{
			result = solve1(box);
			if (result.size() >= 3)
			{
				m_cube_result.push_back(box);
			}
		}
		else
			if (box.sites.size() >= 5)
		{
			auto sites = box.sites;
			vector<pair<Site_3D, double>> sites4;
			for(int j = 0; j<sites.size()-3;j++)
				for(int k = j+1;k<sites.size()-2;k++)
					for (int m = k + 1; m<sites.size()-1; m++)
						for (int n = m + 1; n < sites.size(); n++)
						{
							sites4.clear();
							sites4.push_back(sites.at(j));
							sites4.push_back(sites.at(k));
							sites4.push_back(sites.at(m));
							sites4.push_back(sites.at(n));
							Cube cube = box;
							cube.sites = sites4;
							result = solve1(cube);
							if (result.size() >= 3)
							{
								m_cube_result.push_back(cube);
							}
						}
		}
	}
}
void Vertex_3DApollonius::remove_cube()
{
	vector<Cube> m_box;
	std::cout << m_cube.size() << endl;
	//map<int, int> site_num;
	//ofstream out("cube.txt");
	//for (int i = 0; i < m_cube.size(); ++i)
	//{
	//	Cube cube = m_cube.at(i);
	//	int num = cube.sites.size();
	//	if (site_num.find(num) == site_num.end())
	//	{
	//		site_num[num] = 1;
	//	}
	//	else
	//	{
	//		int count = site_num[num];
	//		site_num[num] = count+1;
	//	}
	//	if (cube.sites.size() == 60)
	//	{
	//		out << "cube:" << endl;
	//		out << "center:"<<cube.pointOfcube_cente << endl;
	//		out << "lenth" << cube.lenthOfcube_sides << endl;
	//		for (int j = 0; j < cube.sites.size(); ++j)
	//		{
	//			out << "v " << cube.sites.at(j).first.site_position<<" "<<cube.sites.at(j).first.weight << endl;
	//		}
	//	}
	//}
	//out.close();
	//for (auto it = site_num.begin(); it != site_num.end(); it++)
	//{
	//	std::cout << it->first << " " << it->second << endl;
	//}
	//getchar();
#pragma omp parallel for
	for (int i = 0; i < m_cube.size(); i++)
	{
		Cube box = m_cube.at(i);
		if (box.sites.size() == 4)
		{
			vector<double> result = solve(box);
			if (result.size() >= 3)
			{
#pragma omp critical
				m_cube_result.push_back(box);
			}
		}
		if (box.sites.size() >= 5 && box.sites.size() <= 5)
		{
			auto sites = box.sites;
			vector<pair<Site_3D, double>> sites4;
			for (int j = 0; j<sites.size() - 3; j++)
				for (int k = j + 1; k<sites.size() - 2; k++)
					for (int m = k + 1; m<sites.size() - 1; m++)
						for (int n = m + 1; n < sites.size(); n++)
						{
							sites4.clear();
							sites4.push_back(sites.at(j));
							sites4.push_back(sites.at(k));
							sites4.push_back(sites.at(m));
							sites4.push_back(sites.at(n));
							Cube cube = box;
							cube.sites = sites4;
							vector<double> result = solve(cube);
							if (result.size() >= 3)
							{
#pragma omp critical
								m_cube_result.push_back(cube);
							}
						}
		}
		if (box.sites.size() > 5)
		{
#pragma omp critical
			m_box.push_back(box);
		}
	}
	m_cube.clear();
	m_cube = m_box;

}
vector<double> Vertex_3DApollonius::solve1(Cube &box)
{
	Solve_Equation solve_equation;
	vector<vector<double>>  site;
	vector<Site> sg;
#pragma omp parallel for
	for (int j = 0; j < box.sites.size(); ++j)
	{
		vector<double> site_one;
		site_one.push_back(box.sites.at(j).first.site_position.hx());
		site_one.push_back(box.sites.at(j).first.site_position.hy());
		site_one.push_back(box.sites.at(j).first.site_position.hz());
		site_one.push_back(box.sites.at(j).first.weight);
#pragma omp critical
		site.push_back(site_one);
		sg.push_back(box.sites.at(j).first);
	}


	double x_center = box.pointOfcube_cente.hx();
	double y_center = box.pointOfcube_cente.hy();
	double z_center = box.pointOfcube_cente.hz();
	Site s_center(box.pointOfcube_cente, 0);

	double box_lenth = box.lenthOfcube_sides;
	double x_start = x_center - box_lenth / 2;
	double x_end = x_center + box_lenth / 2;
	double y_start = y_center - box_lenth / 2;
	double y_end = y_center + box_lenth / 2;
	double z_start = z_center - box_lenth / 2;
	double z_end = z_center + box_lenth / 2;
	vector<double> range = { x_start,x_end,y_start,y_end, z_start,z_end };
	vector<double>  X = { x_center,y_center,z_center };
	vector<double> result = solve_equation.solve(3, X, site, range);
	return result;
}
vector<double> Vertex_3DApollonius::solve(Cube &box)
{
	vector<double> res;
	vector<Site> sites;
	for (int j = 0; j < box.sites.size(); ++j)
	{
		sites.push_back(box.sites.at(j).first);
	}
	vector<double> result = solve_equation(sites);
	double x_center = box.pointOfcube_cente.hx();
	double y_center = box.pointOfcube_cente.hy();
	double z_center = box.pointOfcube_cente.hz();
	double box_lenth = box.lenthOfcube_sides;

	double x_start = x_center - box_lenth / 2;
	double x_end = x_center + box_lenth / 2;
	double y_start = y_center - box_lenth / 2;
	double y_end = y_center + box_lenth / 2;
	double z_start = z_center - box_lenth / 2;
	double z_end = z_center + box_lenth / 2;
	double errro = 1e-5;
	Point center = box.pointOfcube_cente;
	if (result.size()==3)
	{
		double x = result.at(0);
		double y = result.at(1);
		double z = result.at(2);
		Point v(x,y,z);
		
		double distance = Distance_3D(v,center);
			if (x >= x_start - errro && x <= x_end + errro
		&&  y >= y_start - errro && y <= y_end + errro
		&&  z >= z_start - errro && z <= z_end + errro)
	//	if(distance - sqrt(3)*box.lenthOfcube_sides <= 0)
		{
			res = result;
			return res;
		}
	}
	if (result.size() == 6)
	{
		double x1 = result.at(0);
		double y1 = result.at(1);
		double z1 = result.at(2);
		double x2 = result.at(3);
		double y2 = result.at(4);
		double z2 = result.at(5);
		Point v1(x1, y1, z1);
		Point v2(x2, y2, z2);
		double distance1 = Distance_3D(v1, center);
		double distance2 = Distance_3D(v2, center);
		if ((x1 > (x_start - errro)) && (x1 < (x_end + errro))
		&& (y1 > (y_start - errro)) && (y1 <(y_end + errro))
		&& (z1 > (z_start - errro)) && (z1 <(z_end + errro)))
	//	if (distance1 - sqrt(3) * box.lenthOfcube_sides <= 0)
		{
			res.push_back(x1);
			res.push_back(y1);
			res.push_back(z1);
		}
		if ((x2 > (x_start - errro)) && (x2 < (x_end + errro))
			&& (y2 > (y_start - errro)) && (y2 < (y_end + errro))
			&& (z2 > (z_start - errro)) && (z2 < (z_end + errro)))
	//	if (distance2 - sqrt(3) * box.lenthOfcube_sides <= 0)
		{
			res.push_back(x2);
			res.push_back(y2);
			res.push_back(z2);
		}
		return res;
	}
	return res;
}
void Vertex_3DApollonius::refine_cube(double &lenth_flag)
{
	vector<Cube> m_box_new;
	int flag = 0;
	for (int i = 0; i < m_cube.size(); ++i)
	{
		//分割四点以上的cube
		vector<Cube> box_seg;
		if (m_cube.at(i).sites.size() >= 4)
		{
			box_seg = segment_box(m_cube.at(i));
			lenth_flag = m_cube.at(i).lenthOfcube_sides / 2;
			flag++;
			for (int j = 0; j < box_seg.size(); ++j)
			{
				m_box_new.push_back(box_seg.at(j));
			}
		}
	}
	if (flag == 0)
		lenth_flag = 0;

	m_cube.clear();
	for (int i = 0; i < m_box_new.size(); ++i)
	{
		Cube cube_one = m_box_new.at(i);
		set<Site> site;
		for (int j = 0; j < cube_one.sites.size(); ++j)
			site.insert(cube_one.sites.at(j).first);
		site = PKatCube(site, cube_one, false);
		vector<pair<Site_3D, double>> sites;
		for (auto it = site.begin(); it != site.end(); it++)
		{
			sites.push_back(make_pair(*it,0));
		}
		cube_one.sites = sites;
		if (site.size() >= 4)
		{
			m_cube.push_back(cube_one);
		}
	}
	m_box_new.clear();
}
void Vertex_3DApollonius::confirm_vertex()
{
	m_vertex.clear();


#pragma omp parallel for
	for (int i = 0; i < m_cube_result.size(); ++i)
	{
		vector<double> result;
		if(m_cube_result.at(i).sites.size()==4)
			result = solve(m_cube_result.at(i));
		else 
			result = solve1(m_cube_result.at(i));
		if (result.size() > 0)
		{
			Vertex vertex_one;
			for (int j = 0; j < m_cube_result.at(i).sites.size(); ++j)
			{
				vertex_one.site_group.push_back(m_cube_result.at(i).sites.at(j).first);
			}
			for (int j = 0; j < result.size() / 3; ++j)
			{
				Point p(result.at(3*j), result.at(3 * j+1), result.at(3 * j+2));
				vertex_one.vertex_3d = p;
				vertex_one.cube_id = m_cube_result.at(i).id;
				vertex_one.flag = 1;
				bool flag1 = true;
				int count = 0;
				double distance1 = Distance_3D(vertex_one.site_group.at(0).site_position, p) - vertex_one.site_group.at(0).weight;
				for (int j = 1; j < vertex_one.site_group.size(); ++j)
				{
					double distance2 = Distance_3D(vertex_one.site_group.at(j).site_position, p) - vertex_one.site_group.at(j).weight;
					if (abs(distance1 - distance2) > 10e-6)
					{
						flag1 = false;
						break;
					}
				}
				if (flag1 == true)
				{
#pragma omp critical			
					m_vertex.push_back(vertex_one);
				}
			}	
		}
	}	
	cout << m_vertex.size() << endl;

	for (int i = 0; i < m_vertex.size(); ++i)
	{
		vector<Site> sg1 = m_vertex.at(i).site_group;
		Point vertex1 = m_vertex.at(i).vertex_3d;
		double distance = Distance_3D(vertex1, sg1.at(0).site_position) - sg1.at(0).weight;
		for (int j = 0; j < m_sites.size(); ++j)
		{
			double distance1 = Distance_3D(vertex1, m_sites.at(j).site_position) - m_sites.at(j).weight;
			if (distance1 < distance - 1e-6)
			{
				m_vertex.erase(m_vertex.begin() + i);
				i--;
				break;
			}
		}
	}
	cout << m_vertex.size() << endl;

	for (int i = 0; i < m_vertex.size(); ++i)
	{
		int id1 = m_vertex.at(i).cube_id;
		vector<Site> sg1 = m_vertex.at(i).site_group;
		Point vertex1 = m_vertex.at(i).vertex_3d;
		bool flag = true;
		int count = 0;
		for (int j = 0; j < sg1.size(); ++j)
		{
			int id = sg1.at(j).id;
			if (id < 0)
				count++;
		}
		if (count > 3)
		//	flag = false;
		for (int j = i + 1; j < m_vertex.size(); ++j)
		{
			int id2 = m_vertex.at(j).cube_id;
			if (id1 != id2)
			{
				vector<Site> sg2 = m_vertex.at(j).site_group;
				Point vertex2 = m_vertex.at(j).vertex_3d;
				if ((m_vertex.at(i).site_group == m_vertex.at(j).site_group)
					&& Distance_3D(m_vertex.at(i).vertex_3d, m_vertex.at(j).vertex_3d) <10e-6 && (i != j) && flag == true)
				{
					m_vertex.erase(m_vertex.begin() + j);
					--j;
					if (j < 0) j = 0;
					m_vertex.at(j).flag = false;
				}
			}
		}
	}
	vector<Vertex> vertice;
	for (int i = 0; i < m_vertex.size(); ++i)
	{
		if (m_vertex.at(i).flag == true)
		{
			vertice.push_back(m_vertex.at(i));
		}
	}
	m_vertex.clear();
	m_vertex = vertice;
	cout << m_vertex.size() << endl;

	for (int i = 0; i < m_vertex.size(); ++i)
	{
		m_vertex.at(i).id = i;
	}
}
vector<Vertex> Vertex_3DApollonius::GetIntersectionOfVertex(vector<Vertex> vg1, vector<Vertex> vg2)
{
	vector<Vertex> vg;
	vg = vg1;
	for (int i = 0; i < vg2.size(); ++i)
	{
		Vertex v2 = vg2.at(i);

		bool flag = true;
		for (int j = 0; j < vg1.size(); ++j)
		{
			Vertex v1 = vg1.at(j);
			if (v1.vertex_3d == v2.vertex_3d)
			{
				flag = false;
				break;
			}
		}
		if (flag == true)
		{
			vg.push_back(v2);
		}
	}
	return vg;
}
void Vertex_3DApollonius::out_vertex( const char* filename3)
{
	ofstream out(filename3);
	for (int i = 0; i < m_vertex.size(); ++i)
	{
		out << "#" << i << "    vertex:" << m_vertex.at(i).vertex_3d << "    site_group: {";
		for (int j = 0; j < m_vertex.at(i).site_group.size(); ++j)
		{
			out << "(" << m_vertex.at(i).site_group.at(j).id << ")" << "   ";
			//out << "(" << m_vertex.at(i).site_group.at(j).site_position << ")" << "   ";
		}
		out << "}" << endl;
	}
	out.close();
}
void Vertex_3DApollonius::vertex_remove()
{
	for (int i = 0; i < m_vertex.size(); ++i)
	{
		Vertex v = m_vertex.at(i);
		vector<Site> sites = v.site_group;
		int count = 0;
		for (int j = 0; j < sites.size(); ++j)
		{
			Site s = sites.at(j);
			if (s.id == m_sites.at(m_sites.size() - 1).id || s.id == m_sites.at(m_sites.size() - 2).id || s.id == m_sites.at(m_sites.size() - 3).id)
			{
				count++;
			}
		}
		if (count > 0 && count < 4)
		{
			m_vertex.erase(m_vertex.begin()+i);
		}
	}
	for (int i = 0; i < m_vertex.size(); ++i)
	{
		m_vertex.at(i).id = i;
	}
}


