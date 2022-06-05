
#define _CRT_SECURE_NO_WARNINGS
#include "float.h"
#include "Site.h"

typedef Site_3DApollonius::Site_3D  Site;
void Site_3DApollonius::Receive_sites(const char* filename)
{
	ifstream in(filename);
	string str;
	int i = 0;
	while (getline(in, str))
	{
		istringstream line(str);
		double x, y, z, weight;
		line >> x >> y >> z >> weight;
		Site site;
		Point p(x, y, z);
		site.site_position = p;
		site.weight = weight;
		site.id = i;
		m_sites.push_back(site);
		++i;
	}
	in.close();
	int count = 0;
	//for (int i = 0; i < m_sites.size(); ++i)
	//{
	//	Site site1 = m_sites.at(i);
	//	for (int j = i+1; j < m_sites.size(); ++j)
	//	{
	//		Site site2 = m_sites.at(j);
	//		Point P1 = site1.site_position;
	//		Point P2 = site2.site_position;
	//		double distance = sqrt(((P1.hx() - P2.hx()) * (P1.hx() - P2.hx())) +
	//			((P1.hy() - P2.hy()) * (P1.hy() - P2.hy())) +
	//			((P1.hz() - P2.hz()) * (P1.hz() - P2.hz())));
	//		if (distance <= abs(site1.weight - site2.weight))
	//			count++;
	//	}
	//}
	//cout << "hidden site: "<<count << endl;
}

void Site_3DApollonius::add_virtual_site(double boundary)
{
	double xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
	xmin1 = m_sites.at(0).site_position.x();
	xmax1 = m_sites.at(0).site_position.x();
	ymin1 = m_sites.at(0).site_position.y();
	ymax1 = m_sites.at(0).site_position.y();
	zmin1 = m_sites.at(0).site_position.z();
	zmax1 = m_sites.at(0).site_position.z();
	weight_max = m_sites.at(0).weight;

#pragma omp parallel for
	for (int i = 1; i < m_sites.size(); ++i)
	{
		double x, y, z, weight;
		x = m_sites.at(i).site_position.x();
		y = m_sites.at(i).site_position.y();
		z = m_sites.at(i).site_position.z();
		weight = m_sites.at(i).weight;
		if (x < xmin1) xmin1 = x;
		if (x > xmax1) xmax1 = x;
		if (y < ymin1) ymin1 = y;
		if (y > ymax1) ymax1 = y;
		if (z < zmin1) zmin1 = z;
		if (z > zmax1) zmax1 = z;
		if (weight > weight_max) weight_max = weight;
	}
	double xmin, xmax, ymin, ymax, zmin, zmax, xlenth, ylenth, zlenth, boundary_lenth;
	xmin = xmin1;
	xmax = xmax1;
	ymin = ymin1;
	ymax = ymax1;
	zmin = zmin1;
	zmax = zmax1;

	xlenth = (xmax - xmin) * 1 + weight_max * 2;
	ylenth = (ymax - ymin) * 1 + weight_max * 2;
	zlenth = (zmax - zmin) * 1 + weight_max * 2;

	boundary_lenth = xlenth;
	if (ylenth > boundary_lenth)
		boundary_lenth = ylenth;
	if (zlenth > boundary_lenth)
		boundary_lenth = zlenth;
	double x_point = (xmin + xmax) / 2;
	double y_point = (ymin + ymax) / 2;
	double z_point = (zmin + zmax) / 2;

	boundary_lenth = boundary_lenth;

	xlenth = boundary_lenth;
	ylenth = boundary_lenth;
	zlenth = boundary_lenth;

	double x_min = x_point - xlenth * 0.5;
	double y_min = y_point - ylenth * 0.5;
	double z_min = z_point - zlenth * 0.5;
	m_boundary.push_back(x_point - xlenth*0.5 );
	m_boundary.push_back(x_point + xlenth*0.5);
	m_boundary.push_back(y_point - ylenth*0.5);
	m_boundary.push_back(y_point + ylenth*0.5);
	m_boundary.push_back(z_point - zlenth*0.5);
	m_boundary.push_back(z_point + zlenth*0.5);
	virturre_site_size = 3;
	if (virturre_site_size > 1)
	{
		xlenth = xlenth * 1;
		ylenth = ylenth * 1;
		zlenth = zlenth * 1;
		x_min = x_point - xlenth * 0.5;
		y_min = y_point - ylenth * 0.5;
	    z_min = z_point - zlenth * 0.5;
		double x_len = xlenth / (virturre_site_size - 1);
		double y_len = ylenth / (virturre_site_size - 1);
		double z_len = zlenth / (virturre_site_size - 1);
		int id = -1;

		for (int i = 0; i<virturre_site_size; i++)
			for (int j = 0; j<virturre_site_size; j++)
				for (int k = 0; k < virturre_site_size; k++)
				{
					if (i == 0 || j == 0 || k == 0 || i == virturre_site_size - 1 || j == virturre_site_size - 1 || k == virturre_site_size - 1)
					{
						Point p(x_min + i*x_len, y_min + j*y_len, z_min + k*z_len);
						Site site;
						site.id = id;
						site.site_position = p;
						//site.weight = weight_max;
						site.weight = 0;
					//	m_sites.push_back(site);
					//	id--;
					//	num_virture_site++;
					}

				}
	}
	
	//m_boundary ----  x_min,x_max,y_min,y_max,z_min,z_max;
	double x_lenth_bbox = xlenth;
	double y_lenth_bbox = ylenth;
	double z_lenth_bbox = zlenth;
	double para = 10;
	m_cube_boundary.push_back(x_point - x_lenth_bbox * para);
	m_cube_boundary.push_back(x_point + x_lenth_bbox * para);
	m_cube_boundary.push_back(y_point - y_lenth_bbox * para);
	m_cube_boundary.push_back(y_point + y_lenth_bbox * para);
	m_cube_boundary.push_back(z_point - z_lenth_bbox * para);
	m_cube_boundary.push_back(z_point + z_lenth_bbox * para);
}

//bool operator==(const Point  P1, const Point P2)
//{
//	double distance = sqrt((P1.hx() - P2.hx())*(P1.hx() - P2.hx())
//		+ (P1.hy() - P2.hy())*(P1.hy() - P2.hy())
//		+ (P1.hz() - P2.hz())*(P1.hz() - P2.hz()));
//	if (distance < 0.00001)
//		return true;
//	else return  false;
//}
//
//bool operator==(const vector<Site> site1, const vector<Site> site2)
//{
//	if (site1.size() != site2.size())
//		return false;
//	else {
//		int num = 0;
//		for (int i = 0; i < site1.size(); ++i)
//			for (int j = 0; j < site2.size(); ++j)
//			{
//				if (site1.at(i).id == site2.at(j).id)
//					num++;
//			}
//		if (num == site1.size())
//			return true;
//		else
//			return false;
//	}
//}
//
//bool operator==(const Site site1, const Site site2)
//{
//	if (site1.id == site2.id)
//		return true;
//	else
//		return false;
//}

void Site_3DApollonius::GetIntersectOfSitegroup(vector<Site> &sg1, vector<Site> &sg2, vector<Site> &sg)
{
	for (int i = 0; i < sg1.size(); ++i)
	{
		Site site1 = sg1.at(i);
		for (int j = 0; j < sg2.size(); ++j)
		{
			Site  site2 = sg2.at(j);
			if (site1.id == site2.id)
			{
				sg.push_back(site1);
				break;
			}
		}
	}
}
bool Site_3DApollonius::IfSubset(vector<Site> &sg1, vector<Site> &sg2)
{
	vector<Site> sg;
	GetIntersectOfSitegroup(sg1, sg2, sg);
	if (sg.size() == sg1.size() && sg.size() < sg2.size())
		return true;
	else
		return false;
}
bool Site_3DApollonius::IfEquality(vector<Site> &sg1, vector<Site> &sg2)
{
	vector<Site> sg;
	GetIntersectOfSitegroup(sg1, sg2, sg);
	if (sg.size() == sg1.size() && sg.size() == sg2.size())
		return true;
	else
		return false;
}
bool Site_3DApollonius::IfCoplanarity(vector<Site> sg)
{
	if (sg.size() <= 3)
		return true;
	else
	{
		vector<Site> sg1;
		Point p1 = sg.at(0).site_position;
		Point p2 = sg.at(1).site_position;
		Point p3;
		sg1.push_back(sg.at(0));
		sg1.push_back(sg.at(1));
		for (int i = 2; i < sg.size(); ++i)
		{
			sg1.push_back(sg.at(i));
			if (!IfCollineation(sg1))
			{
				p3 = sg.at(i).site_position;
				sg.erase(sg.begin() + i);
				break;
			}
			else
			{
				sg1.pop_back();
			}
		}
		bool flag = true;
		double a = ((p2.y() - p1.y())*(p3.z() - p1.z()) - (p2.z() - p1.z())*(p3.y() - p1.y()));
		double b = ((p2.z() - p1.z())*(p3.x() - p1.x()) - (p2.x() - p1.x())*(p3.z() - p1.z()));
		double c = ((p2.x() - p1.x())*(p3.y() - p1.y()) - (p2.y() - p1.y())*(p3.x() - p1.x()));
		double d = -(a*p1.x() + b*p1.y() + c*p1.z());
		for (int i = 2; i < sg.size(); ++i)
		{
			Point p = sg.at(i).site_position;
			double fx = a*p.hx() + b*p.hy() + c*p.hz() + d;
			if (fx != 0)
				flag = false;
		}
		return flag;
	}

}
bool Site_3DApollonius::IfCollineation(vector<Site> sg)
{
	bool flag = true;
	if (sg.size() > 2)
	{
		double x1 = sg.at(0).site_position.hx();
		double y1 = sg.at(0).site_position.hy();
		double x2 = sg.at(1).site_position.hx();
		double y2 = sg.at(1).site_position.hy();
		if (x1 == x2&&y1 != y2)
		{
			for (int i = 2; i < sg.size(); ++i)
			{
				double x = sg.at(i).site_position.hx();
				double y = sg.at(i).site_position.hy();
				if (x != x1)
				{
					flag = false;
					break;
				}
			}
		}
		if (y1 == y2&&x1 != x2)
		{
			for (int i = 2; i < sg.size(); ++i)
			{
				double x = sg.at(i).site_position.hx();
				double y = sg.at(i).site_position.hy();
				if (y != y1)
				{
					flag = false;
					break;
				}
			}
		}
		if (y1 != y2&&x1 != x2)
			for (int i = 2; i < sg.size(); ++i)
			{
				double x = sg.at(i).site_position.hx();
				double y = sg.at(i).site_position.hy();
				if ((x - x1) / (x2 - x1) != (y - y1) / (y2 - y1))
				{
					flag = false;
					break;
				}
			}
	}
	return flag;
}
