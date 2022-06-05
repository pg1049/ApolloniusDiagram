#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "Cube.h"
typedef Cube_3DApollonius::Cube Cube;
typedef Site_3DApollonius::Site_3D  Site;
using namespace std;

ostream& operator<<(ostream& out, Cube& cube)
{
	
	return out;
}
void Cube_3DApollonius::Generate_Cube(vector<double>& boundary, int sites_size)
{

	double  x_start = boundary.at(0);
	double  x_end = boundary.at(1);
	double y_start = boundary.at(2);
	double y_end = boundary.at(3);
	double z_start = boundary.at(4);
	double z_end = boundary.at(5);
	double x_lenth, y_lenth, z_lenth, cube_lenth, min_lenth,max_length;
	x_lenth = (x_end - x_start);
	y_lenth = (y_end - y_start);
	z_lenth = (z_end - z_start);
	min_lenth = x_lenth;
	max_length = x_lenth;
	if (y_lenth < min_lenth)
		min_lenth = y_lenth;
	if (z_lenth < min_lenth)
		min_lenth = z_lenth;

	if (y_lenth > max_length)
		max_length = y_lenth;
	if (z_lenth > max_length)
		max_length = z_lenth;
	int size = pow(sites_size, 1.0 / 3);
	//if (size > 10)
		cube_lenth = min_lenth / size;
	//cube_lenth = max_lenth / sites_size;
	//else
	//	cube_lenth = min_lenth / 10;
	
	cube_lenth = cube_lenth;
	int max_num = max_length / cube_lenth + 1;
	x_max_cube = max_num + 1;
	y_max_cube = max_num + 1;
	z_max_cube = max_num + 1;

	Point PP((x_start + x_end)*0.5 - x_max_cube*cube_lenth*0.5,
		(y_start + y_end)*0.5 - y_max_cube*cube_lenth*0.5,
		(z_start + z_end)*0.5 - z_max_cube*cube_lenth*0.5);
	source_point = PP;
	numOfcube = x_max_cube*y_max_cube*z_max_cube;


//	ofstream out("E://code//Apollonius_new//result//cube.obj");
	int count = 0;

	for (int i = 0; i < x_max_cube; ++i)
		for (int j = 0; j < y_max_cube; ++j)
			for (int k = 0; k < z_max_cube; ++k)
			{
				Point P(i, j, k);
				Site site_one;
				site_one.id = -INTPTR_MAX;
				Point p(DBL_MAX, DBL_MAX, DBL_MAX);
				site_one.site_position = p;
				site_one.weight = DBL_MAX;
			

				Cube cube;
				cube.id = count;
				cube.indexOfcube = P;
				cube.lenthOfcube_sides = cube_lenth;
				cube.min_distance = DBL_MAX;
				cube.empty = 1;
				double x_center = PP.hx() + (i + 0.5)*cube_lenth;
				double y_center = PP.hy() + (j + 0.5)*cube_lenth;
				double z_center = PP.hz() + (k + 0.5)*cube_lenth;
				Point center(x_center, y_center, z_center);
				cube.pointOfcube_cente = center;
				get_vertice_Cube(cube);
				count++;
				if (i == 0 || j == 0 || k == 0 || i == x_max_cube || j == y_max_cube || k == z_max_cube)
					cube.inOrbound = 0;
				else
					cube.inOrbound = 1;
				m_cube.push_back(cube);
			}
	m_cube_boundary.at(0) = PP.hx();
	m_cube_boundary.at(1) = m_boundary.at(0) + x_max_cube*cube_lenth;
	m_cube_boundary.at(2) = PP.hy();
	m_cube_boundary.at(3) = m_boundary.at(2) + y_max_cube*cube_lenth;
	m_cube_boundary.at(4) = PP.hz();
	m_cube_boundary.at(5) = m_boundary.at(4) + z_max_cube*cube_lenth;
	cout << numOfcube << endl;
	cout << cube_lenth << endl;
}

void Cube_3DApollonius::get_vertice_Cube(Cube &cube)
{
	int i, j, k;
	Point center = cube.pointOfcube_cente;
	double x = center.hx();
	double y = center.hy();
	double z = center.hz();
	double cube_lenth = cube.lenthOfcube_sides;
	Point p1(x - cube_lenth*0.5, y - cube_lenth*0.5, z - cube_lenth*0.5);
	Point p2(x + cube_lenth*0.5, y - cube_lenth*0.5, z - cube_lenth*0.5);
	Point p3(x + cube_lenth*0.5, y + cube_lenth*0.5, z - cube_lenth*0.5);
	Point p4(x - cube_lenth*0.5, y + cube_lenth*0.5, z - cube_lenth*0.5);
	Point p5(x - cube_lenth*0.5, y - cube_lenth*0.5, z + cube_lenth*0.5);
	Point p6(x + cube_lenth*0.5, y - cube_lenth*0.5, z + cube_lenth*0.5);
	Point p7(x + cube_lenth*0.5, y + cube_lenth*0.5, z + cube_lenth*0.5);
	Point p8(x - cube_lenth*0.5, y + cube_lenth*0.5, z + cube_lenth*0.5);

	cube.cube_vertice.clear();
	cube.cube_vertice.push_back(p1);
	cube.cube_vertice.push_back(p2);
	cube.cube_vertice.push_back(p3);
	cube.cube_vertice.push_back(p4);
	cube.cube_vertice.push_back(p5);
	cube.cube_vertice.push_back(p6);
	cube.cube_vertice.push_back(p7);
	cube.cube_vertice.push_back(p8);
}
double Cube_3DApollonius::Distance_point2point(Point P1, Point P2)
{
	return sqrt(((P1.hx() - P2.hx()) * (P1.hx() - P2.hx())) +
		((P1.hy() - P2.hy()) * (P1.hy() - P2.hy())) +
		((P1.hz() - P2.hz()) * (P1.hz() - P2.hz())));
}

Point Cube_3DApollonius::get_SiteInCube(Cube &cube_one, Point p_site)
{
	double box_lenth = cube_one.lenthOfcube_sides;
	double box_lenth_half = box_lenth / 2;
	double x_ref = source_point.hx();
	double y_ref = source_point.hy();
	double z_ref = source_point.hz();
	double x_site = p_site.hx();
	double y_site = p_site.hy();
	double z_site = p_site.hz();
	int x = (x_site - x_ref) / box_lenth;
	int y = (y_site - y_ref) / box_lenth;
	int z = (z_site - z_ref) / box_lenth;
	Point p(x, y, z);
	return p;
}

Cube Cube_3DApollonius::get_CubeByIndex(vector<Cube>& m_cube, Point index)
{
	Cube cube_one;
	int i = index.hx();
	int j = index.hy();
	int k = index.hz();
	int ind = i * y_max_cube*z_max_cube + j * z_max_cube + k;
	cube_one = m_cube.at(ind);
	return cube_one;

}

vector<Point> Cube_3DApollonius::get_CubeNeighbor(Cube &cube_one)
{
	vector<Point> CubeNeighbor;
	Point cube_position = cube_one.indexOfcube;
	int x, y, z;
	x = cube_position.hx();
	y = cube_position.hy();
	z = cube_position.hz();

	Point CubeNeighbor1(x - 1, y, z);
	if (x - 1 >= 0)
		CubeNeighbor.push_back(CubeNeighbor1);
	Point CubeNeighbor2(x + 1, y, z);
	if (x + 1 < x_max_cube)
		CubeNeighbor.push_back(CubeNeighbor2);
	Point CubeNeighbor3(x, y - 1, z);
	if (y - 1 >= 0)
		CubeNeighbor.push_back(CubeNeighbor3);
	Point CubeNeighbor4(x, y + 1, z);
	if (y + 1 < y_max_cube)
		CubeNeighbor.push_back(CubeNeighbor4);
	Point CubeNeighbor5(x, y, z - 1);
	if (z - 1 >= 0)
		CubeNeighbor.push_back(CubeNeighbor5);
	Point CubeNeighbor6(x, y, z + 1);
	if (z + 1 < z_max_cube)
		CubeNeighbor.push_back(CubeNeighbor6);
	return  CubeNeighbor;
}

vector<Cube> Cube_3DApollonius::segment_box(Cube &box)
{
	vector<Cube> box_segment;
	Cube box_seg;
	box_seg.lenthOfcube_sides = box.lenthOfcube_sides / 2;
	box_seg.indexOfcube = box.indexOfcube;
	box_seg.sites = box.sites;
	box_seg.min_distance = DBL_MAX;
	Point center = box.pointOfcube_cente;
	double cube_lenth = box.lenthOfcube_sides / 2;
	double x1, x2, y1, y2, z1, z2;
	x1 = center.hx() - cube_lenth / 2;
	x2 = center.hx() + cube_lenth / 2;
	y1 = center.hy() - cube_lenth / 2;
	y2 = center.hy() + cube_lenth / 2;
	z1 = center.hz() - cube_lenth / 2;
	z2 = center.hz() + cube_lenth / 2;

	Point center1(x1, y1, z1);
	Point center2(x1, y2, z1);
	Point center3(x2, y2, z1);
	Point center4(x2, y1, z1);
	Point center5(x1, y1, z2);
	Point center6(x1, y2, z2);
	Point center7(x2, y2, z2);
	Point center8(x2, y1, z2);

	vector<Point> centers;
	centers.push_back(center1);
	centers.push_back(center2);
	centers.push_back(center3);
	centers.push_back(center4);
	centers.push_back(center5);
	centers.push_back(center6);
	centers.push_back(center7);
	centers.push_back(center8);
	for (int i = 0; i < centers.size(); ++i)
	{
		box_seg.cube_vertice.clear();
		box_seg.pointOfcube_cente = centers.at(i);
		get_vertice_Cube(box_seg);
		box_segment.push_back(box_seg);
	}
	return box_segment;
}

void Cube_3DApollonius::Cube_Inorbound(Cube &cube)
{
	double x_start = m_boundary.at(0);
	double x_end = m_boundary.at(1);
	double y_start = m_boundary.at(2);
	double y_end = m_boundary.at(3);
	double z_start = m_boundary.at(4);
	double z_end = m_boundary.at(5);

	double x_min = cube.pointOfcube_cente.hx() - cube.lenthOfcube_sides / 2;
	double x_max = cube.pointOfcube_cente.hx() + cube.lenthOfcube_sides / 2;

	double y_min = cube.pointOfcube_cente.hy() - cube.lenthOfcube_sides / 2;
	double y_max = cube.pointOfcube_cente.hx() + cube.lenthOfcube_sides / 2;

	double z_min = cube.pointOfcube_cente.hz() - cube.lenthOfcube_sides / 2;
	double z_max = cube.pointOfcube_cente.hx() + cube.lenthOfcube_sides / 2;

	if ((x_start == x_min) || (x_end == x_max) || (y_start == y_min) || (y_end == y_max) || (z_start == z_min) || (z_end == z_max))
	{
		cube.inOrbound = 0;
		cout << cube.inOrbound << endl;
	}
	else
		cube.inOrbound = 1;
}