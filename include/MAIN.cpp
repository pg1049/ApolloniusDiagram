#define _CRT_SECURE_NO_WARNINGS
#define ebsion 1e-5;
#include "Apollonius3D.h"
#include "VectorOperations.h"
#include "cvtlloyd.h"
typedef Apollonius3D::Face_3D Face;
typedef Apollonius3D::Cell_3D Cell;
using namespace std;
int main(int argc, char* argv[])
{
	Point p0(0, 0, 0);
	Point p1(1, 0, 0);
	Point p2(-1, 0, 0);
	Point p3(0, 1, 0);
	Point p4(0, -1, 0);
	Point p5(0, 0, 1);
	Point p6(0, 0, -1);

	Site site0(p0, 0.5);
	Site site1(p1, 0.3);
	Site site2(p2, 0.3);
	Site site3(p3, 0.3);
	Site site4(p4, 0.3);
	Site site5(p5, 0.3);
	Site site6(p6, 0.3);
	vector<Site> m_sites = { site0,site1,site2,site3,site4,site5,site6 };
	Apollonius3D apollonius(m_sites);
	apollonius.open = 0;
	apollonius.run();
	//	apollonius.out_SiteTxt("F://code//Apollonius//video//video3//result" + num + "//site.txt");
//	apollonius.out_VertexObj("F://code//Apollonius//video//video3//vertex.obj");
	apollonius.out_FaceObj("F://code//Apollonius//video//video3//result//face.obj");
	for (int i = 0; i < 20; i++)
	{
		string num = to_string(i);
		if (i < 10)
			num = "0" + num;
		string filename = "F://code//Apollonius//video//video3//result//result" + num + ".txt";
		ofstream out(filename);
		for (int j = 0; j < m_sites.size(); j++)
		{
			out << m_sites.at(j).site_position << " " << m_sites.at(j).weight - i * 0.01 << endl;
		}
		out.close();
	}
//	getchar();



}

