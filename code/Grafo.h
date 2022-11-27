//#pragma once
#include <vector>
using namespace std;
class Sphere_type
{
public:
	Sphere_type() {}
	Sphere_type(Vector3D ve, float ra, Vector3D ppSource, Vector3D ppTarget)
	{ // no requiero centro ya que siempre estara en el origen
		centro = ve; r = ra; Psource = ppSource; Ptarget = ppTarget;
	}


	Vector3D centro, Psource, Ptarget;
	float r;
};

class Edge
{
public:
	Edge()
	{
	}
	Edge(unsigned int v1, unsigned int v2)
	{
		vertex1 = v1; vertex2 = v2; collision = 0; cont_seguimiento = 0; num_collisions = 0; cont_fin = 0;
		ext_path = 0; end = 0; radio_path = 0.0; new_proc = 0;
	}
	unsigned int vertex1, vertex2, num_collisions, collisions[50],  //dir
		num_seguidores, cont_fin;
	float vectAngle, AngleRight, AngleLeft, originAngle, targetAngle, radio_path;
	int gammav1, gammav2, cont_seguimiento;
	bool collision, ext_path, end, new_proc;
	Sphere_type circumf;
	//vector<Sphere_type> sphere;
	Vector3D vect, origin;
	Vector3D seguimiento[60], stop[60]; //podrias requerir mas elementos
};

class Vertex
{
public:
	Vertex()
	{
	}
	Vertex(string iden, int etiq, int nivel)
	{
		name = iden;
		id = etiq;
		level = nivel;
		mark = -1;
		value = 0;
		color[0] = 0; color[1] = 0; color[2] = 0;
	}
	//
	string name;
	int id, level, mark;
	GLfloat centro[2];
	GLfloat color[3];
	vector<int> liga;
	bool value;
};

class Analiza_Par
{
public:
	Analiza_Par()
	{
		pivot = FALSE;
		edge_type = FALSE;
		rol = FALSE;
		follow= FALSE;
	}
	bool pivot, edge_type, rol, follow;
};
