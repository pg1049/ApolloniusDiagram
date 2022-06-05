#include "tool.h"
#include "solve_equation.h"
Point_3 get_Fun_face_normal(Point_3 s1, Point_3 s2, Point_3 p)
{
	double f1 = sqrt((p.hx() - s1.hx()) * (p.hx() - s1.hx()) + (p.hy() - s1.hy()) * (p.hy() - s1.hy()) + (p.hz() - s1.hz()) * (p.hz() - s1.hz()));
	double f2 = sqrt((p.hx() - s2.hx()) * (p.hx() - s2.hx()) + (p.hy() - s2.hy()) * (p.hy() - s2.hy()) + (p.hz() - s2.hz()) * (p.hz() - s2.hz()));
	double fx = (p.hx() - s1.hx()) / f1 - (p.hx() - s2.hx()) / f2;
	double fy = (p.hy() - s1.hy()) / f1 - (p.hy() - s2.hy()) / f2;
	double fz = (p.hz() - s1.hz()) / f1 - (p.hz() - s2.hz()) / f2;
	Point normal(fx, fy, fz);
	return normal;
};
Point_3 get_normal(vector<Point_3> face_one)
{
	Point_3 P1 = face_one.at(0);
	Point_3 P2 = face_one.at(1);
	Point_3 P3 = face_one.at(2);

	Point_3 normal1(P1.hx() - P2.hx(), P1.hy() - P2.hy(), P1.hz() - P2.hz());
	Point_3 normal2(P2.hx() - P3.hx(), P2.hy() - P3.hy(), P2.hz() - P3.hz());

	double x1 = normal1.hx();
	double y1 = normal1.hy();
	double z1 = normal1.hz();
	double x2 = normal2.hx();
	double y2 = normal2.hy();
	double z2 = normal2.hz();
	Point normal((y1*z2 - y2*z1), -(x1*z2 - x2*z1), (x1*y2 - x2*y1));
	return normal;
}
double dot_product(Point_3 P1, Point_3 P2)
{
	double result;
	result = P1.hx()*P2.hx() + P1.hy()*P2.hy() + P1.hz()*P2.hz();
	return result;
}
vector<double> cross_product(vector<double> v1, vector<double> v2)
{
	double x1 = v1.at(0);
	double y1 = v1.at(1);
	double z1 = v1.at(2);
	double x2 = v2.at(0);
	double y2 = v2.at(1);
	double z2 = v2.at(2);
	double x = y1*z2 - y2*z1;
	double y = -(x1*z2 - x2*z1);
	double z = x1*y2 - x2*y1;
	vector<double> vec = { x,y,z };
	return vec;
}
vector<double> solve_equation(vector<Site> m_site, Point p)
{
	Solve_Equation solve_equation;
	vector<vector<double>>  site;
	for (int j = 0; j < m_site.size(); ++j)
	{
		vector<double> site_one;
		site_one.push_back(m_site.at(j).site_position.hx());
		site_one.push_back(m_site.at(j).site_position.hy());
		site_one.push_back(m_site.at(j).site_position.hz());
		site_one.push_back(m_site.at(j).weight);
		site.push_back(site_one);
	}
	vector<double>  X;
	X.push_back(p.hx());
	X.push_back(p.hy());
	X.push_back(p.hz());
	vector<double> result = solve_equation.solve(3, X, site);
	return result;
}
vector<double> solve_equation(vector<Site> m_site, Point p, vector<double> range)
{
	Solve_Equation solve_equation;
	vector<vector<double>>  site;
	for (int j = 0; j < m_site.size(); ++j)
	{
		vector<double> site_one;
		site_one.push_back(m_site.at(j).site_position.hx());
		site_one.push_back(m_site.at(j).site_position.hy());
		site_one.push_back(m_site.at(j).site_position.hz());
		site_one.push_back(m_site.at(j).weight);
		site.push_back(site_one);
	}
	vector<double>  X;
	X.push_back(p.hx());
	X.push_back(p.hy());
	X.push_back(p.hz());
	vector<double> result = solve_equation.solve(3, X, site, range);
	return result;
}
vector<int> GetInsectOfVector(vector<int> v1, vector<int> v2)
{
	vector<int> result;
	for (int i = 0; i < v1.size(); ++i)
	{
		for (int j = 0; j < v2.size(); j++)
		{
			if (v1.at(i) == v2.at(j))
				result.push_back(v1.at(i));
		}
	}
	return result;
}

Point_3  codeRotateByX(Point_3 p, double rx)
{
	double x1 = p.x();//将变量拷贝一次，保证&x == &outx这种情况下也能计算正确
	double y1 = p.y();
	double z1 = p.z();
	double outy = cos(rx) * y1 - sin(rx) * z1;
	double outz = cos(rx) * z1 + sin(rx) * y1;
	double outx = x1;
	Point out_p(outx, outy, outz);
	return out_p;
}
Point_3  codeRotateByY(Point_3 p, double ry)
{
	double x1 = p.x();//将变量拷贝一次，保证&x == &outx这种情况下也能计算正确
	double y1 = p.y();
	double z1 = p.z();
	double outx = cos(ry) * x1 + sin(ry) * z1;
	double outz = cos(ry) * z1 - sin(ry) * x1;
	double outy = y1;
	Point out_p(outx, outy, outz);
	return out_p;
}

Point_3  codeRotateByZ(Point_3 p, double rz)
{
	double x1 = p.x();//将变量拷贝一次，保证&x == &outx这种情况下也能计算正确
	double y1 = p.y();
	double z1 = p.z();
	double outx = cos(rz) * x1 - sin(rz) * y1;
	double outy = sin(rz) * x1 + cos(rz) * y1;
	double outz = z1;
	Point_3 out_p(outx, outy, outz);
	return out_p;
}

void  RotateXYZ(vector<double> v1, vector<double> &theta, vector<Point>& boundary_point)
{
	vector<double> vXOY = { v1.at(0), v1.at(1), 0 };
	vector<double> vXOZ = { v1.at(0), 0, v1.at(2) };
	vector<double> vYOZ = { 0,  v1.at(1), v1.at(2) };
	vector<double> vx = { 1,0,0 };
	vector<double> vy = { 0,1,0 };
	vector<double> vz = { 0,0,1 };
	vector<Point> points;
	points.push_back(boundary_point.at(0));
	points.push_back(boundary_point.at(boundary_point.size() / 3));
	points.push_back(boundary_point.at(boundary_point.size() / 3 * 2));
	Point normal1 = get_normal(points);
	double theta1, theta2, theta3;
	theta1 = Angle(vXOY, vx);
	if (v1.at(1) >= 0)
		theta1 = -theta1;
	theta2 = Angle(vXOZ, vz);
	if (v1.at(0) >= 0)
		theta2 = -theta2;
	theta3 = Angle(vYOZ, vz);
	if (v1.at(1) <= 0)
		theta3 = -theta3;

	if ((v1.at(0) == 0) && (v1.at(1) != 0))
	{
		for (int i = 0; i < boundary_point.size(); ++i)
		{
			boundary_point.at(i) = codeRotateByX(boundary_point.at(i), theta3);
		}
	}
	else
		if (v1.at(0) != 0 && v1.at(1) == 0)
		{
			for (int i = 0; i < boundary_point.size(); ++i)
			{
				boundary_point.at(i) = codeRotateByY(boundary_point.at(i), theta2);
			}
		}
		else
			if (v1.at(0) != 0 || v1.at(1) != 0)
			{

				if (theta3 > PI / 2)
					theta3 = -(PI - theta3);
				if (theta3 < -PI / 2)
					theta3 = PI + theta3;

				Point n(v1.at(0), v1.at(1), v1.at(2));
				Point normal = codeRotateByX(n, theta3);
				vXOZ.clear();
				vXOZ.push_back(normal.hx());
				vXOZ.push_back(normal.hy());
				vXOZ.push_back(normal.hz());
				theta2 = Angle(vXOZ, vz);
				if (v1.at(0) >= 0)
					theta2 = -theta2;

				if (theta2 > PI / 2)
					theta2 = -(PI - theta2);
				if (theta2 < -PI / 2)
					theta2 = PI + theta2;
				for (int i = 0; i < boundary_point.size(); ++i)
				{
					boundary_point.at(i) = codeRotateByX(boundary_point.at(i), theta3);
					boundary_point.at(i) = codeRotateByY(boundary_point.at(i), theta2);
				}
			}
	theta.push_back(theta1);
	theta.push_back(theta2);
	theta.push_back(theta3);
}

void  ReRotateXYZ(vector<double> v1, vector<double> theta, vector<Point>& boundary_point)
{
	double theta1 = -theta.at(0);
	double theta2 = -theta.at(1);
	double theta3 = -theta.at(2);
	//	cout <<"rerotate:"<< theta1 / PI * 180 << " " << theta2 / PI * 180 << "  " << theta3 / PI * 180 << endl;
	if ((v1.at(0) == 0) && (v1.at(1) != 0))
	{
		for (int i = 0; i < boundary_point.size(); ++i)
		{
			boundary_point.at(i) = codeRotateByX(boundary_point.at(i), theta3);
		}
	}
	else
		if (v1.at(0) != 0 && v1.at(1) == 0)
		{
			for (int i = 0; i < boundary_point.size(); ++i)
			{
				boundary_point.at(i) = codeRotateByY(boundary_point.at(i), theta2);
			}
		}
		else
			if (v1.at(0) != 0 || v1.at(1) != 0)
			{
				for (int i = 0; i < boundary_point.size(); ++i)
				{
					boundary_point.at(i) = codeRotateByY(boundary_point.at(i), theta2);
					boundary_point.at(i) = codeRotateByX(boundary_point.at(i), theta3);
				}
			}
}


vector<Point_2> RandomPointInPolygon(Polygon2 poly, int num)
{
	int count = 0;
	priority_queue<Event_poly> m_poly;
	vector<Point_2> m_point;
	Event_poly eve;
	eve.poly = poly;
	eve.area = abs(poly.area());
	m_poly.push(eve);
	do
	{
		Event_poly event = m_poly.top();
		m_poly.pop();
		Polygon2 poly1 = event.poly;
		int vertexCount = poly1.size();
		Point_2 center = compute2DPolygonCentroid(poly1);
	//	m_point.push_back(center);
		for (int i = 0; i < vertexCount; ++i)
		{
			Point_2 p1 = poly1[i%vertexCount];
			Point_2 p2 = poly1[(i + 1) % vertexCount];
			Polygon2 poly2;
			poly2.push_back(center);
			poly2.push_back(p1);
			poly2.push_back(p2);
			Event_poly eve1;
			eve1.poly = poly2;
			eve1.area = abs(poly2.area());
			m_poly.push(eve1);
		}
	} while (m_poly.size() < num+1);
	count = 0;
//	cout <<"m_poly.size()"<< m_poly.size() << endl;
	while (count<num)
	{
		Event_poly event = m_poly.top();
		m_poly.pop();
		Polygon2 poly1 = event.poly;
		Point_2 center = compute2DPolygonCentroid(poly1);
		m_point.push_back(center);
		count++;
	}
	//for (int i = 0; i < m_point.size(); ++i)
	//{
	//	cout << m_point.at(i) << endl;
	//	if (CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(poly.vertices_begin(), poly.vertices_end(), m_point.at(i), K()))
	//	{
	//		cout << "success" << endl;
	//	}
	//}
		
	//cout << m_point.size() << endl;
//	getchar();
	return m_point;
}

Point_2 compute2DPolygonCentroid(Polygon2 poly)
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

double Distance_3D(Point P1, Point P2)
{
	return sqrt(((P1.hx() - P2.hx()) * (P1.hx() - P2.hx())) +
		((P1.hy() - P2.hy()) * (P1.hy() - P2.hy())) +
		((P1.hz() - P2.hz()) * (P1.hz() - P2.hz())));
}

double Distance_2D(Point_2 P1, Point_2 P2)
{
	return sqrt(((P1.hx() - P2.hx()) * (P1.hx() - P2.hx())) +
		((P1.hy() - P2.hy()) * (P1.hy() - P2.hy())));
}

double MaxLenthOfTriangle(Point_2 p1, Point_2 p2, Point_2 p3)
{
	double len1 = Distance_2D(p1, p2);
	double len2 = Distance_2D(p1, p3);
	double len3 = Distance_2D(p2, p3);
	double max = len1;
	if (len2 > max)
		max = len2;
	if (len3 > max)
		max = len3;
	return max;
}
Point_3 RotateByVector(Point_3 old_point, Point_3  vn, double theta)
{
	double old_x = old_point.hx();
	double old_y = old_point.hy();
	double old_z = old_point.hz();
	double vx = vn.hx();
	double vy = vn.hy();
	double vz = vn.hz();
	double r = theta;
	double c = cos(r);
	double s = sin(r);
	double new_x = (vx*vx*(1 - c) + c) * old_x + (vx*vy*(1 - c) - vz*s) * old_y + (vx*vz*(1 - c) + vy*s) * old_z;
	double new_y = (vy*vx*(1 - c) + vz*s) * old_x + (vy*vy*(1 - c) + c) * old_y + (vy*vz*(1 - c) - vx*s) * old_z;
	double new_z = (vx*vz*(1 - c) - vy*s) * old_x + (vy*vz*(1 - c) + vx*s) * old_y + (vz*vz*(1 - c) + c) * old_z;
	return Point_3(new_x, new_y, new_z);
}

void delay(int seconds)
{
	time_t start, current;

	time(&start);

	do
	{
		time(&current);
	} while ((current - start) < seconds);
}
//vector<double> solve_equation(vector<Site> m_site)
//{
//	vector<double>  result;
//	Site site1 = m_site.at(0);
//	Site site2 = m_site.at(1);
//	Site site3 = m_site.at(2);
//	Site site4 = m_site.at(3);
//
//	double x4 = site4.site_position.hx();
//	double y4 = site4.site_position.hy();
//	double z4 = site4.site_position.hz();
//	double r4 = site4.weight;
//
//	double x1 = site1.site_position.hx() - x4;
//	double y1 = site1.site_position.hy() - y4;
//	double z1 = site1.site_position.hz() - z4;
//	double r1 = site1.weight - r4;
//
//	double x2 = site2.site_position.hx() - x4;
//	double y2 = site2.site_position.hy() - y4;
//	double z2 = site2.site_position.hz() - z4;
//	double r2 = site2.weight - r4;
//
//	double x3 = site3.site_position.hx() - x4;
//	double y3 = site3.site_position.hy() - y4;
//	double z3 = site3.site_position.hz() - z4;
//	double r3 = site3.weight - r4;
//	if ((x1 == 0 && x2 == 0 && x3 == 0) || (y1 == 0 && y2 == 0 && y3 == 0) || (z1 == 0 && z2 == 0 && z3 == 0))
//	{
//		x1 = x1 + 1;
//		y1 = y1 + 1;
//		z1 = z1 + 1;
//		x2 = x2 + 1;
//		y2 = y2 + 1;
//		z2 = z2 + 1;
//		x3 = x3 + 1;
//		y3 = y3 + 1;
//		z3 = z3 + 1;
//		x4 = x4 - 1;
//		y4 = y4 - 1;
//		z4 = z4 - 1;
//	
//	}
//	Eigen::Matrix3d m;
//	m(0, 0) = x1;
//	m(0, 1) = y1;
//	m(0, 2) = z1;
//
//	m(1, 0) = x2;
//	m(1, 1) = y2;
//	m(1, 2) = z2;
//
//	m(2, 0) = x3;
//	m(2, 1) = y3;
//	m(2, 2) = z3;
//
//	double A_det;
//	A_det = m.determinant();
//	double d1 = (x1*x1 + y1*y1 + z1*z1 - r1*r1) / 2;
//	double d2 = (x2*x2 + y2*y2 + z2*z2 - r2*r2) / 2;
//	double d3 = (x3*x3 + y3*y3 + z3*z3 - r3*r3) / 2;
////	cout << m << endl;
////	cout << A_det << endl;
//	if (abs(A_det)!=0)
//	{
//		Eigen::Matrix3d m1;
//		m1 = m;
//		m1(0, 0) = d1;
//		m1(1, 0) = d2;
//		m1(2, 0) = d3;
//		double B1_det = m1.determinant();
//		m1(0, 0) = r1;
//		m1(1, 0) = r2;
//		m1(2, 0) = r3;
//		double C1_det = m1.determinant();
//
//		m1 = m;
//		m1(0, 1) = d1;
//		m1(1, 1) = d2;
//		m1(2, 1) = d3;
//		double B2_det = m1.determinant();
//		m1(0, 1) = r1;
//		m1(1, 1) = r2;
//		m1(2, 1) = r3;
//		double C2_det = m1.determinant();
//
//		m1 = m;
//		m1(0, 2) = d1;
//		m1(1, 2) = d2;
//		m1(2, 2) = d3;
//		double B3_det = m1.determinant();
//		m1(0, 2) = r1;
//		m1(1, 2) = r2;
//		m1(2, 2) = r3;
//		double C3_det = m1.determinant();
//		double a = C1_det * C1_det + C2_det * C2_det + C3_det * C3_det - A_det *A_det;
//		double b = -2 * (B1_det * C1_det + B2_det * C2_det + B3_det * C3_det);
//		double c = B1_det * B1_det + B2_det * B2_det + B3_det * B3_det;
//
//		double delta = b*b - 4 * a*c;
//		if (abs(a) == 0)
//		{
//			double r = -c / b;
//			Eigen::Matrix<double, 3, 1> D;
//			D << d1 - r1*r, d2 - r2*r, d3 - r3*r;
//			auto X = m.fullPivLu().solve(D);
//			double x = X(0, 0) + x4;
//			double y = X(1, 0) + y4;
//			double z = X(2, 0) + z4;
//			result.push_back(x);
//			result.push_back(y);
//			result.push_back(z);
//			return result;
//		}
//		if (delta >= 0)
//		{
//			double r_1 = (-b + sqrt(delta)) / (2 * a);
//			double r_2 = (-b - sqrt(delta)) / (2 * a);
//			Eigen::Matrix<double, 3, 1> D1, D2;
//			D1 << d1 - r1*r_1, d2 - r2*r_1, d3 - r3*r_1;
//			D2 << d1 - r1*r_2, d2 - r2*r_2, d3 - r3*r_2;
//			auto X1 = m.fullPivLu().solve(D1);
//			auto X2 = m.fullPivLu().solve(D2);
//			double x1 = X1(0, 0) + x4;
//			double y1 = X1(1, 0) + y4;
//			double z1 = X1(2, 0) + z4;
//			double x2 = X2(0, 0) + x4;
//			double y2 = X2(1, 0) + y4;
//			double z2 = X2(2, 0) + z4;
//			result.push_back(x1);
//			result.push_back(y1);
//			result.push_back(z1);
//			if (!(X1 == X2))
//			{
//				result.push_back(x2);
//				result.push_back(y2);
//				result.push_back(z2);
//			}
//		}
//		return result;
//	}
//	return result;
//}

vector<double> solve_equation( vector<Site> m_site)
{
	vector<double>  result;
	Site site1 = m_site.at(0);
	Site site2 = m_site.at(1);
	Site site3 = m_site.at(2);
	Site site4 = m_site.at(3);

	double x4 = site4.site_position.hx();
	double y4 = site4.site_position.hy();
	double z4 = site4.site_position.hz();
	double r4 = site4.weight;

	double x1 = site1.site_position.hx() - x4;
	double y1 = site1.site_position.hy() - y4;
	double z1 = site1.site_position.hz() - z4;
	double r1 = site1.weight - r4;

	double x2 = site2.site_position.hx() - x4;
	double y2 = site2.site_position.hy() - y4;
	double z2 = site2.site_position.hz() - z4;
	double r2 = site2.weight - r4;

	double x3 = site3.site_position.hx() - x4;
	double y3 = site3.site_position.hy() - y4;
	double z3 = site3.site_position.hz() - z4;
	double r3 = site3.weight - r4;

	double w1 = (x1*x1 + y1*y1 + z1*z1 - r1*r1) / 2;
	double w2 = (x2*x2 + y2*y2 + z2*z2 - r2*r2) / 2;
	double w3 = (x3*x3 + y3*y3 + z3*z3 - r3*r3) / 2;
	bool flag = true;

	Eigen::Matrix3d m, m1, m2, m3;

	m << x1, y1, z1,
		x2, y2, z2,
		x3, y3, z3;
	m1 << y1, z1, r1,
		y2, z2, r2,
		y3, z3, r3;
	m2 << x1, z1, r1,
		x2, z2, r2,
		x3, z3, r3;
	m3 << x1, y1, r1,
		x2, y2, r2,
		x3, y3, r3;

	double A_det, A_det1, A_det2, A_det3;
	A_det = m.determinant();
	A_det1 = m1.determinant();
	A_det2 = m2.determinant();
	A_det3 = m3.determinant();
	Eigen::Matrix<double, 3, 1> w;
	w << w1, w2, w3;
//	cout << A_det << " " << A_det1 << " " << A_det2 << " " << A_det3 << endl;
	if (A_det != 0)
	{
		Eigen::Matrix<double, 3, 1> r;
		r << -r1, -r2, -r3;
		auto m_inverse = m.inverse();
		auto alpha = m_inverse * r;
		auto beta = m_inverse * w;
		double a, b, c;
		a = alpha(0, 0)*alpha(0, 0) + alpha(1, 0)*alpha(1, 0) + alpha(2, 0)*alpha(2, 0) - 1;
		b = 2 * (alpha(0, 0)*beta(0, 0) + alpha(1, 0)*beta(1, 0) + alpha(2, 0)*beta(2, 0));
		c = beta(0, 0)*beta(0, 0) + beta(1, 0)*beta(1, 0) + beta(2, 0)*beta(2, 0);
		double delta = b*b - 4 * a*c;
		if (abs(a) == 0)
		{
			double ro = -c / b;
			if (ro > 0)
			{
				Eigen::Matrix<double, 3, 1> D;
				D = w + ro*r;
				auto X = m.fullPivLu().solve(D);
				double x = X(0, 0) + x4;
				double y = X(1, 0) + y4;
				double z = X(2, 0) + z4;
				result.push_back(x);
				result.push_back(y);
				result.push_back(z);
			}
			return result;
		}
		if (delta >= 0)
		{
			double ro_1 = (-b + sqrt(delta)) / (2 * a);
			double ro_2 = (-b - sqrt(delta)) / (2 * a);
		//	cout << ro_1 << " " << ro_2 << " " << endl;
			Eigen::Matrix<double, 3, 1> D1, D2;
			D1 = w + ro_1*r;
			D2 = w + ro_2*r;
			auto X1 = m.fullPivLu().solve(D1);
			auto X2 = m.fullPivLu().solve(D2);
			double x1 = X1(0, 0) + x4;
			double y1 = X1(1, 0) + y4;
			double z1 = X1(2, 0) + z4;
			double x2 = X2(0, 0) + x4;
			double y2 = X2(1, 0) + y4;
			double z2 = X2(2, 0) + z4;
			if (ro_1 >0)
			{
				result.push_back(x1);
				result.push_back(y1);
				result.push_back(z1);
			}
			if (ro_1!= ro_2&&ro_2 >0)
			{
				result.push_back(x2);
				result.push_back(y2);
				result.push_back(z2);
			}

		}
	}
	else if (A_det1 != 0)
	{
		Eigen::Matrix<double, 3, 1> x;
		x << -x1, -x2, -x3;
		auto m_inverse = m1.inverse();
		auto alpha = m_inverse * x;
		auto beta = m_inverse * w;
		double a, b, c;
		a = alpha(0, 0)*alpha(0, 0) + alpha(1, 0)*alpha(1, 0) - alpha(2, 0)*alpha(2, 0) + 1;
		b = 2 * (alpha(0, 0)*beta(0, 0) + alpha(1, 0)*beta(1, 0) - alpha(2, 0)*beta(2, 0));
		c = beta(0, 0)*beta(0, 0) + beta(1, 0)*beta(1, 0) - beta(2, 0)*beta(2, 0);
		double delta = b*b - 4 * a*c;
		if (abs(a) == 0)
		{
			double ro = -c / b;
			Eigen::Matrix<double, 3, 1> D;
			D = w + ro*x;
			auto X = m1.fullPivLu().solve(D);
			double x = ro + x4;
			double y = X(0, 0) + y4;
			double z = X(1, 0) + z4;
			if (ro > 0)
			{
				result.push_back(x);
				result.push_back(y);
				result.push_back(z);
			}

			return result;
		}
		if (delta >= 0)
		{
			double ro_1 = (-b + sqrt(delta)) / (2 * a);
			double ro_2 = (-b - sqrt(delta)) / (2 * a);
			Eigen::Matrix<double, 3, 1> D1, D2;
			D1 = w + ro_1*x;
			D2 = w + ro_2*x;
			auto X1 = m1.fullPivLu().solve(D1);
			auto X2 = m1.fullPivLu().solve(D2);
			double x1 = ro_1 + x4;
			double y1 = X1(0, 0) + y4;
			double z1 = X1(1, 0) + z4;
			double x2 = ro_2 + x4;
			double y2 = X2(0, 0) + y4;
			double z2 = X2(1, 0) + z4;
			if (X1(2, 0) > 0)
			{
				result.push_back(x1);
				result.push_back(y1);
				result.push_back(z1);
			}
			if (!(ro_1 == ro_2) && (X2(2, 0) > 0))
			{
				result.push_back(x2);
				result.push_back(y2);
				result.push_back(z2);
			}

		}
	}
	else if (A_det2 != 0)
	{
		Eigen::Matrix<double, 3, 1> r;
		r << -y1, -y2, -y3;
		auto m_inverse = m2.inverse();
		auto alpha = m_inverse * r;
		auto beta = m_inverse * w;
		double a, b, c;
		a = alpha(0, 0)*alpha(0, 0) + alpha(1, 0)*alpha(1, 0) - alpha(2, 0)*alpha(2, 0) + 1;
		b = 2 * (alpha(0, 0)*beta(0, 0) + alpha(1, 0)*beta(1, 0) - alpha(2, 0)*beta(2, 0));
		c = beta(0, 0)*beta(0, 0) + beta(1, 0)*beta(1, 0) - beta(2, 0)*beta(2, 0);
		double delta = b*b - 4 * a*c;
		if (abs(a) == 0)
		{
			double ro = -c / b;
			Eigen::Matrix<double, 3, 1> D;
			D = w + ro*r;
			auto X = m2.fullPivLu().solve(D);
			double x = X(0, 0) + x4;
			double y = ro + y4;
			double z = X(1, 0) + z4;
			if (ro > 0)
			{
				result.push_back(x);
				result.push_back(y);
				result.push_back(z);
			}

			return result;
		}
		if (delta >= 0)
		{
			double ro_1 = (-b + sqrt(delta)) / (2 * a);
			double ro_2 = (-b - sqrt(delta)) / (2 * a);
			Eigen::Matrix<double, 3, 1> D1, D2;
			D1 = w + ro_1*r;
			D2 = w + ro_2*r;
			auto X1 = m2.fullPivLu().solve(D1);
			auto X2 = m2.fullPivLu().solve(D2);
			double x1 = X1(0, 0) + x4;
			double y1 = ro_1 + y4;
			double z1 = X1(1, 0) + z4;
			double x2 = X2(0, 0) + x4;
			double y2 = ro_2 + y4;
			double z2 = X2(1, 0) + z4;
			if (X1(2, 0) > 0)
			{
				result.push_back(x1);
				result.push_back(y1);
				result.push_back(z1);
			}
			if (!(ro_1 == ro_2) && (X2(2, 0) > 0))
			{
				result.push_back(x2);
				result.push_back(y2);
				result.push_back(z2);
			}

		}
	}
	else if (A_det3 != 0)
	{
		Eigen::Matrix<double, 3, 1> r;
		r << -z1, -z2, -z3;
		auto m_inverse = m3.inverse();
		auto alpha = m_inverse * r;
		auto beta = m_inverse * w;
		double a, b, c;
		a = alpha(0, 0)*alpha(0, 0) + alpha(1, 0)*alpha(1, 0) - alpha(2, 0)*alpha(2, 0) + 1;
		b = 2 * (alpha(0, 0)*beta(0, 0) + alpha(1, 0)*beta(1, 0) - alpha(2, 0)*beta(2, 0));
		c = beta(0, 0)*beta(0, 0) + beta(1, 0)*beta(1, 0) - beta(2, 0)*beta(2, 0);
		double delta = b*b - 4 * a*c;
		if (abs(a) == 0)
		{
			double ro = -c / b;
			Eigen::Matrix<double, 3, 1> D;
			D = w + ro*r;
			auto X = m3.fullPivLu().solve(D);
			double x = X(0, 0) + x4;
			double y = X(1, 0) + y4;
			double z = ro + z4;
			result.push_back(x);
			result.push_back(y);
			result.push_back(z);
			return result;
		}
		if (delta >= 0)
		{
			double ro_1 = (-b + sqrt(delta)) / (2 * a);
			double ro_2 = (-b - sqrt(delta)) / (2 * a);
			Eigen::Matrix<double, 3, 1> D1, D2;
			D1 = w + ro_1*r;
			D2 = w + ro_2*r;
			auto X1 = m3.fullPivLu().solve(D1);
			auto X2 = m3.fullPivLu().solve(D2);
			double x1 = X1(0, 0) + x4;
			double y1 = X1(1, 0) + y4;
			double z1 = ro_1 + z4;
			double x2 = X2(0, 0) + x4;
			double y2 = X2(1, 0) + y4;
			double z2 = ro_2 + z4;
			if (X1(2, 0) > 0)
			{
				result.push_back(x1);
				result.push_back(y1);
				result.push_back(z1);
			}
			if (!(ro_1 == ro_2) && (X2(2, 0) > 0))
			{
				result.push_back(x2);
				result.push_back(y2);
				result.push_back(z2);
			}

		}
	}
	return result;
}