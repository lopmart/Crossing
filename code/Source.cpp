#include <glew.h>
#include <freeglut.h>
#include "stdlib.h"
#include "stdio.h"
#include <string.h> //for fonts
#include <fstream> //to save into a file
#include <iostream>
#include <math.h>
#include <sstream>
#include <iomanip>
//#include <vector>
#include <algorithm>
#include <conio.h>
//----------------------------------Creados


#define PI 3.1416
//#define PI 3.1415926535897932384626433832795
//----------------------------------Creados
#include <D:/Visual2019/glm/glm.hpp>
using glm::mat4;
using glm::vec3;
#include <D:/Visual2019/glm/gtc/matrix_transform.hpp>
//mouse
#include "Vector3D.h"
#include "Grafo.h"
float mouseX, mouseY;

#include "Camera.h"
#include "MousePicker.h"
MousePicker picker;
CCamera Camera;
//-----------------------------
int ColisionNodo, NumNodes = 5000, NumEdges = 5300, NumEdgesSG[100], NumNodosSG[100], NumSG = -1;
int ScreenWidth = 512, ScreenHeight = 512;
//int ScreenWidth = 1000, ScreenHeight = 1000;
int Num[100]; //el indice de vSG en donde se empezara a desplegar en pantalla
Vertex v[5000], vSG[100][1000]; //100 SubGrafos de a 5000 nodos c/u
Edge edg[5300], edgSG[100][1000]; //100 SubGrafos de a 15000 aristas c/u
unsigned int MaxNivel = 0, NodosPorNivel[1000];
vector<int> PorNivel[1000]; //1000 niveles en G maximo
unsigned int NumLayers = 12;
short int capa = -1, MaxNumSG = 1;
unsigned int NivelLeader, CapaActual;
bool display_graph, display_nodo_label = 0, display_gamma = 0;
//int rootNode;
unsigned int NivelBegin[100], NivelEnd[100]; //100 SubGrafos maximo
//-----------------------------------------Luminosidad
//Add ambient light
GLfloat ambientColor[] = { 0.2f, 0.2f, 0.2f, 1.0f };
GLfloat specular[] = { 1.0f, 1.0f, 1.0f , 1.0f };
// Create light components
GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat position[] = { -1.5f, 1.0f, 2.0f, 1.0f };

//Add positioned light
GLfloat lightColor0[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat lightPos0[] = { 4.0f, 0.0f, 8.0f, 1.0f };
//Add directed light
GLfloat lightColor1[] = { 0.5f, 0.2f, 0.2f, 1.0f };
//Coming from the direction (-1, 0.5, 0.5)
GLfloat lightPos1[] = { -1.0f, 0.5f, 0.5f, 0.0f };
//MATERIAL
float mcolor[] = { 1.0f, 0.0f, 0.0f, 1.0f };
float colorBlue[] = { 0.0f, 0.0f, 1.0f, 1.0f };
float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
//--------------Camera
float transZ = -25.0; //-174; //-20; 
//float transZ = 1.4f;
//glm::vec3 CameraPos = glm::vec3(0, 0, 0);
glm::vec3 CameraPos = glm::vec3(10, 0, 0);
//----------------------escalamiento anillos
unsigned int sizeAnilloInterno = 3;
float MaxZoom;
//variables para ColisionRoutines.h
struct LigasIntersec { int nodo1, nodo2; } EdgeIntersec[100];
Edge edgTemp[100]; //para <e> o <d>
unsigned int EdgeInter[200], EdgeInterse[100];
unsigned int EdgeIntersecNum = 0, NumIntersec = 0; //ligas entre capas
unsigned int EdgeIntersecNumAro = 0, NumIntersecAro = 0; //ligas en el anillo
unsigned int NumColision00, NumColision01, NumColision10, NumColision11;
unsigned int ciclo00, ciclo01, ciclo10, ciclo11;
//====2-packing process
unsigned int sumaTotal = 0;
//====retomando colisiones el 21 de marzo
unsigned int Giro = 0;
unsigned int GiroMin = 0;
//unsigned int MaxNumColGamma = 100;
unsigned int MinNumColisions = 100;
bool despliega_par = 0, despliega_curva = 0;
unsigned int NumEdgeInterse, aristasBGcollision;
Vector3D LineSegment0, LineSegment1, LineSegment2;
unsigned int CircleEdge[20];
bool checaCircles = 0;
bool ext_path = 0; //arista ext_path que falto concluir su camino cuando su pivote finalizo
bool new_proc = 0; //arista pendiente
bool ext_path_new_proc = 0; //arista ext_path que falto concluir su camino cuando su pivote pending finalizo

bool pivot_running = 0; //comienza proceso
bool arranca = 1; //para los parametros de Intersec en Detecta_Colision_con_Aristas_Gama(...)
 //==== Analysis.h
//Analiza_Par s[20];
vector<Vector3D> pointCircle; //para despliegue
bool display_points_in_circle = 0;
//+Vector3D circle_intersection_point[50]; //solo guardare el primer punto,a. El segundo es vertex2
vector<float> Angles;
struct subject {
	float angle;
	int edge;
	bool label;
	int end;
	bool red;
	subject(float angle, int edge, bool label, int end, bool red) {
		this->angle = angle;
		this->edge = edge;
		this->label = label;
		this->end = end;
		this->red = red;
	}
};
vector<subject> edge_angles;
vector<Analiza_Par> edg_angles;
struct Interval {
	int start, end;
};
struct {
	unsigned int NumEdg, count;
}match[10];
unsigned int Num_match = 0;
//regresando al programa despues del analisis
struct {
	unsigned int num_collisions, cont_seguimiento, rotations, tipo, num_collisions_red;
	unsigned int inner_circles;
	float radio_path, t_original, t_phase1, t_phase2;
	bool ext_path;
}report[15];
bool sentido_rot = 0;
unsigned int Num_aristas_pendientes = 0;
float dist[15];
unsigned int crossing_black_black = 0;
//Tabla para recorridos
int aristas_a_procesar[20];
bool sentido_de_recorridos[2048][20];
int Num_Aristas_Procesar, recorido_a_realizar;
unsigned int BetaCirc = 0, BetaCircAdj = 0;
//
vector<int> randomRow[50];
bool imprime = 0;
//
#include "RutinasGML.h"
#include "ColisionRoutines.h"
#include "Analysis.h"

//=================================================
void reshape(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat)height / (GLfloat)width, 1.0, 440.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawstring(void* font, float x, float y, char* str)
{
	char* ch;
	glRasterPos3f(x, y, 0.0);
	for (ch = str; *ch; ch++)
	{
		glutBitmapCharacter(font, *ch);
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawEdgeBeta(int vertex1, int vertex2, int colour, int thick, unsigned int NumAr)
{
	glLineWidth(thick);
	glColor3f(0.0f, 0.0f, colour);

	glBegin(GL_LINES);
	glVertex2f(v[vertex1].centro[0], v[vertex1].centro[1]);
	glVertex2f(v[vertex2].centro[0], v[vertex2].centro[1]);
	glEnd();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawEdgeAlgo(float colour, int thick, unsigned int NumAr)
{
	Vector3D vect, vect1;
	float r1 = (MaxNivel - CapaActual - 1) * 5.0 - 1; //radio intermedio
	vect = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1]; vect.unitize();
	vect1 = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1]; vect1.unitize();

	glLineWidth(thick);
	glColor3f(0, colour, 0.0f);

	if (edg[NumAr].cont_seguimiento > 0) {
		glBegin(GL_LINES);
		//for (int i = 0; i < (edg[NumAr].cont_fin - 1); i++) {

		glVertex2f(0, 0);
		//glVertex2f(edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y);
		glVertex2f(r1 * vect.x, r1 * vect.y);
		glColor3f(0.5, 0.0f, 0.5);
		glVertex2f(0, 0);
		//glVertex2f(edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y);
		glVertex2f(r1 * vect1.x, r1 * vect1.y);

		glColor3f(1, 0.0f, 0.0f);
		glVertex2f(r1 * vect.x, r1 * vect.y);
		glVertex2f(r1 * vect1.x, r1 * vect1.y);
		//}
		//glVertex2f(edg[NumAr].stop[edg[NumAr].cont_fin - 1].x, edg[NumAr].stop[edg[NumAr].cont_fin - 1].y);
		//glVertex2f(v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1]);
		glEnd();
	}
}//edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1]
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawEdgeGamma(int vertex1, int vertex2, int colour, int thick, unsigned int NumAr)
{
	glLineWidth(thick);
	glColor3f(colour, 0.0f, 0.0f);

	if (edg[NumAr].cont_fin > 1) {
		glBegin(GL_LINES);
		for (int i = 0; i < (edg[NumAr].cont_fin - 1); i++) {
			glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
			glVertex2f(edg[NumAr].stop[i + 1].x, edg[NumAr].stop[i + 1].y);
		}
		glVertex2f(edg[NumAr].stop[edg[NumAr].cont_fin - 1].x, edg[NumAr].stop[edg[NumAr].cont_fin - 1].y);
		glVertex2f(v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1]);
		glEnd();
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawEdge(int vertex1, int vertex2, int colour, int thick, unsigned int NumAr)
{
	glLineWidth(thick);
	glColor3f(colour, 0.0f, 0.0f);
	if (!edg[NumAr].end) { //aun no finaliza
		glBegin(GL_LINES);
		glVertex2f(v[vertex1].centro[0], v[vertex1].centro[1]);
		glVertex2f(v[vertex2].centro[0], v[vertex2].centro[1]);
		glEnd();
	}
	/*if (edg[NumAr].cont_seguimiento > 0 ) {
		glColor3f(0, 0.0f, 0.0f);
		glBegin(GL_LINES);
		glVertex2f(edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y);
		glVertex2f(edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y);
		for (int i = 0; i < (edg[NumAr].cont_seguimiento - 1); i++) {
			glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
			glVertex2f(edg[NumAr].stop[i + 1].x, edg[NumAr].stop[i + 1].y);
		}
		glEnd();
	}*/
	/*else if (edg[NumAr].cont_fin > 1) {
		glColor3f(1, 0.0f, 0.0f);
		glBegin(GL_LINES);
		for (int i = 0; i < (edg[NumAr].cont_fin - 1); i++) {
			glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
			glVertex2f(edg[NumAr].stop[i + 1].x, edg[NumAr].stop[i + 1].y);
		}
		glVertex2f(edg[NumAr].stop[edg[NumAr].cont_fin-1].x, edg[NumAr].stop[edg[NumAr].cont_fin-1].y);
		glVertex2f(v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1]);
		glEnd();
	}*/
	glLineWidth(1);

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawEdgeParametric(int vertex1, int vertex2, int colour, int thick)
{
	Vector3D vect;
	vect.x = v[vertex2].centro[0] - v[vertex1].centro[0];
	vect.y = v[vertex2].centro[1] - v[vertex1].centro[1];
	vect.z = 0;

	glLineWidth(thick);
	glColor3f(0.0f, 0.0f, colour);
	glBegin(GL_LINES);
	glVertex2f(v[vertex1].centro[0], v[vertex1].centro[1]);
	glVertex2f(v[vertex1].centro[0] + 1.0 * vect.x, v[vertex1].centro[1] + 1.0 * vect.y);
	glEnd();
	glLineWidth(1);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void particiona_circle_v1(void) {
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vect, vect1;
	float xc, yc, i1, incr, xcR, ycR;
	float r1, temp;
	//double temp;
	unsigned int cuentaGrados = 0, NumParticiones = 0, puntoSegmento = 0;

	r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda

	double angulo, result;
	unsigned int NumAr = 0;
	i1 = 0.0;
	incr = 2.0 / 360;
	//rpath = 1.0;
	r1 = r1 + 0.5;
	//NumParticiones = 360 / (10 + CapaActual);
	temp = 360.0 / (10.0 + CapaActual);
	NumParticiones = temp;

	if (CapaActual > 42) { NumParticiones = 13; temp = 27.692308; }
	else if (CapaActual > 39) { NumParticiones = 15; temp = 24.0; }
	else if (CapaActual > 35) { NumParticiones = 20; temp = 18.0; }
	else if (CapaActual > 31) { NumParticiones = 25; temp = 14.4; }
	else if (CapaActual > 25) { NumParticiones = 30; temp = 12.0; }
	else if (CapaActual > 19) { NumParticiones = 35; temp = 10.28; }
	else if (CapaActual > 14) { NumParticiones = 40; temp = 9.0; }
	else if (CapaActual > 7) { NumParticiones = 45; temp = 8.0; }
	else { NumParticiones = 50; temp = 7.2; }
	//
	edg[NumAr].cont_seguimiento = 0;

	double param, resultSin, resultCos;
	param = 0;
	puntoSegmento = 0;
	glLineWidth(2);
	glColor3f(0.86, 0.46f, 0.2f);
	glBegin(GL_LINES);
	for (int i = 0; i < NumParticiones; i++)
	{
		resultSin = sin(param * PI / 180);
		resultCos = cos(param * PI / 180);
		xc = r1 * resultCos;
		yc = r1 * resultSin;
		glVertex2f(0, 0);
		glVertex2f(xc, yc);

		if (i > 0) {
			glVertex2f(xcR, ycR);
			glVertex2f(xc, yc);
		}
		xcR = xc; ycR = yc;
		puntoSegmento++;
		param = param + temp; //grados
	}
	//edg[NumAr].cont_seguimiento = puntoSegmento;
	glEnd();
	glLineWidth(1);
	//
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void particiona_circle() {
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vect, vect1;
	float xc, yc, i1, incr, xcR, ycR;
	float r1, rpath;
	unsigned int cuentaGrados = 0, NumParticiones = 0;

	r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
	//el vector de la arista correspondiente
	//vect = edg[NumAr].vect;

	double angulo, result;

	i1 = 0.0;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI, si tuviera mas nodos habria que ajustar
	//i1 = edg[NumAr].AngleLeft * incr;
	rpath = 1.0;
	r1 = r1 + rpath;
	NumParticiones = 360 / (10 + CapaActual);
	//NumParticiones = 10;
	glLineWidth(2);
	//glColor3f(0.1, 0.3, 0.8f);
	glColor3f(0.86, 0.46f, 0.2f);
	glBegin(GL_LINES);
	while (cuentaGrados < NumParticiones) {
		xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
		yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
		glVertex2f(0, 0);
		glVertex2f(r1 * xc, r1 * yc);
		i1 += (10 + CapaActual) * incr;
		xcR = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
		ycR = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
		glVertex2f(r1 * xc, r1 * yc);
		glVertex2f(r1 * xcR, r1 * ycR);
		cuentaGrados++;
	}
	glEnd();
	glLineWidth(1);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void draw_particiones_creadas() {
	unsigned int NumAr = 0; // EdgeInterse[NumEdgeInterse];

	glBegin(GL_LINES);
	for (int i = 0; i < (edg[NumAr].cont_seguimiento - 1); i++) {
		glVertex2f(edg[NumAr].seguimiento[i].x, edg[NumAr].seguimiento[i].y);
		glVertex2f(edg[NumAr].seguimiento[i + 1].x, edg[NumAr].seguimiento[i + 1].y);
	}
	//cierras el ciclo
	glVertex2f(edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y);
	glVertex2f(edg[NumAr].seguimiento[0].x, edg[NumAr].seguimiento[0].y);
	glEnd();

}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void draw_particiones_creadas_todas(unsigned int NumAr) {

	glBegin(GL_LINES);
	if (edg[NumAr].cont_fin > 0) {
		glColor3f(0.0, 0.0, 1.0f);
		for (int i = 1; i < edg[NumAr].cont_fin; i++) {
			glVertex2f(edg[NumAr].stop[i - 1].x, edg[NumAr].stop[i - 1].y);
			glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
		}
		//union del utimo segmento con el primero
		glVertex2f(edg[NumAr].stop[edg[NumAr].cont_fin - 1].x, edg[NumAr].stop[edg[NumAr].cont_fin - 1].y);
		glVertex2f(edg[NumAr].stop[0].x, edg[NumAr].stop[0].y);
	}//if cont_fin > 0
	glEnd();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawCartesianPlane(unsigned int NumAr) {
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vect, vect1;
	float xc, yc, i1, incr, xcR, ycR;
	float r1;
	//Ejes XZ

	//el vector de la arista correspondiente
	vect = edg[NumAr].vect;
	glColor3f(0.0, 1.0f, 1.0f);
	glLineWidth(1);
	//Continúo con el producto punto
	double angulo, result;

	//verificando vector tangente Circunf
	glLineWidth(3);
	i1 = 0.0;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI, si tuviera mas nodos habria que ajustar
	i1 = edg[NumAr].AngleLeft * incr;
	xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
	yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
	//Right
	i1 = edg[NumAr].AngleRight * incr;
	xcR = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
	ycR = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
	//if (display_nodo_label) { //son las bases verdes de las aristas, tangentes al círculo mayor del par en cuestión
		//glColor3f(0, 0.533f, 0.0f);
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2f(v[edg[NumAr].vertex1].centro[0], v[edg[NumAr].vertex1].centro[1]);
	glVertex2f(v[edg[NumAr].vertex1].centro[0] + xc, v[edg[NumAr].vertex1].centro[1] + yc);
	glVertex2f(v[edg[NumAr].vertex1].centro[0], v[edg[NumAr].vertex1].centro[1]);
	glVertex2f(v[edg[NumAr].vertex1].centro[0] + xcR, v[edg[NumAr].vertex1].centro[1] + ycR);
	glEnd();
	//}
	glLineWidth(1);

	/*if (checaCircles) {
		printf("drawCartesianPlane( )\n ");
		printf("EdgeInterse[%d]: ", NumAr);
		printf("end=%d, cont_seguimiento=%d, cont_fin=%d \n", edg[NumAr].end, edg[NumAr].cont_seguimiento, edg[NumAr].cont_fin);
	}*/
	//glColor3f(0.1, 0.3, 0.8f);
	glBegin(GL_LINES);
	if (edg[NumAr].collision) {
		glColor3f(0, 0, 1.0f); glLineWidth(2);
		if (edg[NumAr].end) { glColor3f(0.1, 0.0, 0.8f); glLineWidth(1); } //finalizo
	//si no has finalizado despliega el vector de origen
		/*if (!edg[NumAr].end) {
			glVertex2f(0, 0);
			glVertex2f(edg[NumAr].origin.x, edg[NumAr].origin.y);
		}*/
		/*if (edg[NumAr].cont_seguimiento > 0) {
			glColor3f(0.0, 1.0, 0.0f);
			glVertex2f(edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y);
			glVertex2f(edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x, edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y);
			for (int i = 0; i < (edg[NumAr].cont_seguimiento - 1); i++) {
				glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
				glVertex2f(edg[NumAr].stop[i + 1].x, edg[NumAr].stop[i + 1].y);
			}
		}
		 else {*/
		if (edg[NumAr].cont_fin > 0) {
			glColor3f(0.0, 0.0, 0.0f); //para que los paths sean de color negro en lugar de azules
			//glVertex2f(edg[NumAr].stop[edg[NumAr].cont_fin - 1].x, edg[NumAr].stop[edg[NumAr].cont_fin - 1].y);
			//glVertex2f(edg[NumAr].seguimiento[edg[NumAr].cont_fin - 1].x, edg[NumAr].seguimiento[edg[NumAr].cont_fin - 1].y);
			for (int i = 1; i < (edg[NumAr].cont_fin); i++) {
				glVertex2f(edg[NumAr].stop[i - 1].x, edg[NumAr].stop[i - 1].y);
				glVertex2f(edg[NumAr].stop[i].x, edg[NumAr].stop[i].y);
			}
		}//if cont_fin > 0
		//}//else 

	}//if collision

	glEnd();

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void drawCircle(float radio, unsigned int NumAr)
{
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	float xc, yc, i1, incr;
	i1 = 0.0;
	incr = 2.0 / 80; //80 nodos para la vuelta entera 2 PI
	glColor3f(1, 0.0f, 0.0f);
	glBegin(GL_LINES);
	for (unsigned int i = 0; i < 80; i++)
	{
		xc = plane1[0] * radio * cos(i1 * PI) + (double)plane2[0] * radio * sin(i1 * PI);
		yc = plane1[1] * radio * cos(i1 * PI) + (double)plane2[1] * radio * sin(i1 * PI);
		glVertex3f(xc, yc, 0.0);
		i1 += incr;
	}
	//source
	glVertex3f(v[edg[NumAr].vertex1].centro[0], v[edg[NumAr].vertex1].centro[1], 0.0);
	//glVertex3f(v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1], 0.0);
	glVertex3f(edg[NumAr].circumf.Psource.x, edg[NumAr].circumf.Psource.y, 0.0);
	//target
	glVertex3f(v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1], 0.0);
	glVertex3f(edg[NumAr].circumf.Ptarget.x, edg[NumAr].circumf.Ptarget.y, 0.0);
	glEnd();
	glColor3f(0, 0.0f, 0.0f);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void display(void)
{
	float xc, yc; //posicion numero de nodos
	float xcLevel = 0, ycLevel = 0; //posicion numero de capas
	char collisionStr[20];
	char result1[30];
	char* text1 = new char[30];
	unsigned int Inicio_Capa, Fin_Capa;
	Vector3D temp;

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslatef(0, 0, transZ);
	//:::::::::::::::::::::::::::::::::::::::::::::::: 
	glColor3f(0.8, 0.8, 0.8);
	//glColor3f(0, 0, 0);
	glPointSize(3.0);
	if (display_graph)
	{ //de la capa 0 a la CapaActual+1 
		//Inicio_Capa = CapaActual; 
		Inicio_Capa = 0; Fin_Capa = MaxNivel - 1;
		if (despliega_par)  Fin_Capa = CapaActual + 1; //comentas
	  //:::::::::::::::::::::::::::::::::::::::::::::::: Nodes
		for (unsigned int k = Inicio_Capa; k < Fin_Capa + 1; k++)
		{
			for (std::vector<int>::iterator it = PorNivel[k].begin(); it != PorNivel[k].end(); ++it)
			{
				glColor3fv(v[*it].color);
				xc = v[*it].centro[0]; yc = v[*it].centro[1];
				glBegin(GL_POINTS);
				glVertex3f(xc, yc, 0.0);
				glEnd();

				/*if (display_nodo_label)
				{
					sprintf_s(collisionStr, 20, "%d", v[*it].id);
					strcpy_s(result1, _countof(result1), collisionStr);
					for (int i = 0; i < 5; i++)
						text1[i] = result1[i];
					text1[4] = '\0';
					drawstring(GLUT_BITMAP_HELVETICA_12, xc, yc, text1);
				}*/
			}
			sprintf_s(collisionStr, 20, " Layers: %d, %d", CapaActual, CapaActual + 1);
			strcpy_s(result1, _countof(result1), collisionStr);
			for (int i = 0; i < 20; i++)
				text1[i] = result1[i];
			text1[25] = '\0';
			glColor3f(1.0, 0.0, 0.0); //despliega en rojo el numero de nivel
			xcLevel += 0.25;
		}
		//drawstring(GLUT_BITMAP_HELVETICA_12, 0, transZ * 0.52, text1);

		//:::::::::::::::::::::::::::::::::::::::::::::::: Edges	
		glColor3f(0, 0, 0);

		if (despliega_par) { //<a>
			for (unsigned int i = 0; i < NumEdges; i++)
			{
				if (v[edg[i].vertex1].level > CapaActual) {
					if (edg[i].cont_fin > 0) drawCartesianPlane(i);
					else drawEdge(edg[i].vertex1, edg[i].vertex2, edg[i].collision, 1, i);
				}
			}//for i
		}
		else {
			//::::::::::::::::::::::::::::::::::::::::::::::
			for (unsigned int i = 0; i < EdgeIntersecNum; i++) //gamma (red and black)
			{
				if (EdgeInterse[i] == EdgeInterse[NumEdgeInterse]) {
					if (edg[EdgeInterse[NumEdgeInterse]].cont_fin < 1) {//aun sin recorrido
						glLineStipple(1, 0x3F07);
						glEnable(GL_LINE_STIPPLE);
						drawEdgeParametric(edg[EdgeInterse[NumEdgeInterse]].vertex1, edg[EdgeInterse[NumEdgeInterse]].vertex2, edg[EdgeInterse[NumEdgeInterse]].collision, 2);
						glDisable(GL_LINE_STIPPLE);
					}// if cont_fin < 1
					/*glLineWidth(2);
					glBegin(GL_LINES);
					   glVertex2f(0, 0);
					   glVertex2f(edg[EdgeInterse[NumEdgeInterse]].origin.x, edg[EdgeInterse[NumEdgeInterse]].origin.y);
					glEnd();
					glLineWidth(1);*/

					//drawEdgeAlgo(0.85, 2, EdgeInterse[NumEdgeInterse]);
				}
				else {
					if (edg[EdgeInterse[i]].cont_fin < 1) //aun sin recorrido
						drawEdgeBeta(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInterse[i]].collision, 1, EdgeInterse[i]);
				}
				drawCartesianPlane(EdgeInterse[i]);
				//draw_particiones_creadas_todas(EdgeInterse[i]);
				//etiqueta vertices
				if (display_nodo_label) {
					xc = v[edg[EdgeInterse[i]].vertex1].centro[0]; yc = v[edg[EdgeInterse[i]].vertex1].centro[1];
					sprintf_s(collisionStr, 20, "%d", v[edg[EdgeInterse[i]].vertex1].id);
					strcpy_s(result1, _countof(result1), collisionStr);
					for (int i = 0; i < 5; i++)
						text1[i] = result1[i];
					text1[4] = '\0';
					drawstring(GLUT_BITMAP_HELVETICA_12, xc, yc, text1);
					//etiqueta vertices
					xc = v[edg[EdgeInterse[i]].vertex2].centro[0]; yc = v[edg[EdgeInterse[i]].vertex2].centro[1];
					sprintf_s(collisionStr, 20, "%d", v[edg[EdgeInterse[i]].vertex2].id);
					strcpy_s(result1, _countof(result1), collisionStr);
					for (int i = 0; i < 5; i++)
						text1[i] = result1[i];
					text1[4] = '\0';
					drawstring(GLUT_BITMAP_HELVETICA_12, xc, yc, text1);
				}
			}

			/*for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) //Beta (aristas que forman el circulo)
			{
				drawEdgeBeta(edg[EdgeInter[i]].vertex1, edg[EdgeInter[i]].vertex2, 0, 1, EdgeInter[i]);
			}*/
		}//else  despliega_par
		//...if (checaCircles) particiona_circle(EdgeInterse[NumEdgeInterse]);
		//if (checaCircles) particiona_circle_v1();

	} //if (display_graph)
	//checaCircles = 0;
	//::::::::::::::::::::::::::::::::::::::::::::::
	glMatrixMode(GL_MODELVIEW);
	glColor3f(0, 0, 0);
	//glFlush();
	glutSwapBuffers();
	glutPostRedisplay();
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: comparativo para ordenar
bool compareInterval(Interval i1, Interval i2)
{
	return (i1.start < i2.start);
}
bool compareStruct(subject i1, subject i2)
{
	return (i1.angle > i2.angle);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: init()
void init()
{

	//:::::::::::::::::::::::::::::::::::::::::::::::: Luminosidad
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);

	glShadeModel(GL_SMOOTH);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	//
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);

	glMateriali(GL_FRONT, GL_AMBIENT, 0.1);
	glMateriali(GL_FRONT, GL_DIFFUSE, 0.396);
	glMateriali(GL_FRONT, GL_SPECULAR, 0.297254);
	glMateriali(GL_FRONT, GL_SHININESS, 12.8);

	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//----------------
	for (unsigned int i = 0; i < 100; i++) { NodosPorNivel[i] = 0; Num[i] = 0; }
	for (unsigned int i = 0; i < 40; i++)
	{
		NivelBegin[i] = 5 * i + 1; NivelEnd[i] = 5 * (i + 1);
	}

	CameraPos.z = transZ;
	picker = MousePicker(CameraPos);
	picker.update(CameraPos);
	srand((unsigned)time(0));
	//----------------
	printf("Cargando archivo GML\n");
	//GMLread("D:/semestres/2022/sem1-2022/GraphDrawing (recordando)/gml files/50Outerplanar/Outerplanar50_10.gml");
	GMLread("D:/semestres/2022/sem1-2022/GraphDrawing (recordando)/gml files/20Outerplanar/Outerplanar20_4.gml");
	//
	//GMLread("D:/OpenGL Programs/sem1-2021/GraphDrawingV3/outerPlanar/20Outerplanar2.gml");
	//GMLread("D:/OpenGL Programs/sem1-2021/GraphDrawingV3/outerPlanar/25Outerplanar1.gml");
	AdjustOuterGraph(0, MaxNivel - 1);
	NivelLeader = 0;
	//Display_Level_End = MaxNivel;
	display_graph = 1;
	//CapaActual = MaxNivel - 2;
	CapaActual = 18; //48; //23; //18; //33
	MaxZoom = -440; //-250.0;
	printf("\n CapaMasProfunda=%d, CapaMenosProfunda=%d, MaxNivel=%d\n", CapaActual + 1, CapaActual, MaxNivel);
	printf("Creacion de aristas ---------------------------------------------------------------------------\n");
	//creacion de aristas tipo gamma
	GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
	//creacion de aristas tipo beta
	GeneraLineasEnCapa(CapaActual + 1);
	printf("Deteccion de colisiones ---------------------------------------------------------------------------\n");
	ColisionEntreLineasDeDiferenteCapa();
	//ColisionEntreLineasDeDiferenteCapaVsLineasCapaMasProfunda();
	detecta_colision_con_circulos_interiores();
	printf("[%d]. NumIntersec (%d) +  NumIntersecAro (%d) = %d \n", Giro, NumIntersec, NumIntersecAro, NumIntersec + NumIntersecAro);
	MinNumColisions = NumIntersec + NumIntersecAro;
	NumEdgeInterse = 0; aristasBGcollision = 0;
	//for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		//if (edg[EdgeInterse[i]].collision) Calcula_vector_origin(EdgeInterse[i]);
	//}
	//colision con anillos internos
	printf("Estoy en Capas: %d, %d \n", CapaActual, CapaActual + 1);
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		printf("%d. edg[%d]. ", i, EdgeInterse[i]);
		Dist_Edge_Origin(EdgeInterse[i], i);
	}
	
	randomRow[0].push_back(29);   randomRow[0].push_back(10);   randomRow[0].push_back(5);   randomRow[0].push_back(2);   randomRow[0].push_back(45);   randomRow[0].push_back(59);
	randomRow[1].push_back(3);   randomRow[1].push_back(6);   randomRow[1].push_back(4);
	randomRow[2].push_back(228);   randomRow[2].push_back(244);   randomRow[2].push_back(126);   randomRow[2].push_back(34);   randomRow[2].push_back(28);   randomRow[2].push_back(109);   randomRow[2].push_back(144);   randomRow[2].push_back(95);
	randomRow[3].push_back(27);   randomRow[3].push_back(8);   randomRow[3].push_back(9);   randomRow[3].push_back(17);   randomRow[3].push_back(21);
	randomRow[4].push_back(26);   randomRow[4].push_back(11);   randomRow[4].push_back(18);   randomRow[4].push_back(30);   randomRow[4].push_back(5);
	randomRow[5].push_back(6);   randomRow[5].push_back(12);   randomRow[5].push_back(13);   randomRow[5].push_back(15);
	randomRow[6].push_back(7);   randomRow[6].push_back(20);   randomRow[6].push_back(16);   randomRow[6].push_back(5);   randomRow[6].push_back(28);
	randomRow[7].push_back(18);   randomRow[7].push_back(2);   randomRow[7].push_back(22);   randomRow[7].push_back(16);   randomRow[7].push_back(13);
	randomRow[8].push_back(15);   randomRow[8].push_back(15);   randomRow[8].push_back(4);   randomRow[8].push_back(8);
	randomRow[9].push_back(2);   randomRow[9].push_back(3);   randomRow[9].push_back(5);
	randomRow[10].push_back(29);   randomRow[10].push_back(23);   randomRow[10].push_back(7);   randomRow[10].push_back(12);   randomRow[10].push_back(17);
	randomRow[11].push_back(62);   randomRow[11].push_back(39);   randomRow[11].push_back(63);   randomRow[11].push_back(44);   randomRow[11].push_back(58);   randomRow[11].push_back(61);
	randomRow[12].push_back(59);   randomRow[12].push_back(29);   randomRow[12].push_back(48);   randomRow[12].push_back(61);   randomRow[12].push_back(53);   randomRow[12].push_back(37);
	randomRow[13].push_back(0);   randomRow[13].push_back(25);   randomRow[13].push_back(17);   randomRow[13].push_back(15);   randomRow[13].push_back(18);
	randomRow[14].push_back(7);   randomRow[14].push_back(50);   randomRow[14].push_back(54);   randomRow[14].push_back(11);   randomRow[14].push_back(25);   randomRow[14].push_back(60);
	randomRow[15].push_back(5);   randomRow[15].push_back(9);   randomRow[15].push_back(3);   randomRow[15].push_back(1);
	randomRow[16].push_back(56);   randomRow[16].push_back(6);   randomRow[16].push_back(42);   randomRow[16].push_back(7);   randomRow[16].push_back(36);   randomRow[16].push_back(55);
	randomRow[17].push_back(1);   randomRow[17].push_back(3);
	randomRow[18].push_back(7);   randomRow[18].push_back(13);   randomRow[18].push_back(4);   randomRow[18].push_back(5);   randomRow[18].push_back(11);
}
//=================================================
void keyboard(unsigned char key, int x, int y)
{
	unsigned int count = 0, index;
	float tmpx[2], tmpy[2];
	bool inter = 0;
	unsigned int NumColisionOld, NumColisionNew;
	float distance, r1;
	switch (key) {
	case '+':
		checaCircles = !checaCircles;
		break;
	case '-':
		sizeAnilloInterno--;
		ScaleRingLayer(CapaActual + 1, 0.8);
		//for(int i=CapaActual; i<MaxNivel; i++)
			//ScaleRingLayer(i, 0.8);
		break;
	case 'z':
		printf("Estoy en Capas: %d, %d \n", CapaActual, CapaActual + 1);
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
			printf("%d. edg[%d]. ", i, EdgeInterse[i]);
			Dist_Edge_Origin(EdgeInterse[i], i);
		}

		count = CapaActual + 2; r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4;
		if (CapaActual + 1 < MaxNivel - 1)
			while (count < MaxNivel) {
				r1 = r1 - 5;
				printf("Calculando Interseccion de circulo Capa=%d con r=%2.2f \n", count, r1);
				for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
					printf("%d. edg[%d]. ", i, EdgeInterse[i]);
					//Circle_Edge_Cross(EdgeInterse[i], r1);
					if (dist[i] < r1) { printf(" cross "); report[i].inner_circles = report[i].inner_circles + 2; }
					printf("\n");
				}
				count++; printf("\n");
			}//while
		break;
	case 'x':
		if (despliega_par) { //<a>
			printf("CapaActual=%d\n", CapaActual);
			for (unsigned int i = 0; i < NumEdges; i++)
			{
				if (v[edg[i].vertex1].level <= CapaActual + 2 || v[edg[i].vertex2].level <= CapaActual + 2) {
					//if (edg[i].cont_fin > 0) {}
					//else {
					printf("edg[%d], vertex1=%d(%d), vertex2=%d(%d)\n", i, edg[i].vertex1, v[edg[i].vertex1].level, edg[i].vertex2, v[edg[i].vertex2].level);
					//}
				}
			}//for i
		}
		break;
	case 'i': //--------------------------------------------------------------------------------------
		//count = pow(2, Num_Aristas_Procesar);
		if (recorido_a_realizar < pow(2, Num_Aristas_Procesar)) {
			printf("==============================recorrido =%d.   \n", recorido_a_realizar);
			for (int j = 0; j < Num_Aristas_Procesar; j++) {
				printf("  >  >  >  >  >  >  >  >  >%d(%d), \n", aristas_a_procesar[j], sentido_de_recorridos[recorido_a_realizar][j]);
				sentido_rot = sentido_de_recorridos[recorido_a_realizar][j];
				NumEdgeInterse = aristas_a_procesar[j];
				Calcula_vector_origin(EdgeInterse[NumEdgeInterse]);

				Start_Running();
			}
			printf("   \n");

			recorido_a_realizar++;
		}//if
		break;
	case 'g': //--------------------------------------------------------------------------------------
		Num_Aristas_Procesar = 0; r1 = 0.0;
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
			printf("%d. edg[%d]. ", i, EdgeInterse[i]);
			if (edg[EdgeInterse[i]].collision) //si es arista gamma roja 
			{
				edg[EdgeInterse[i]].radio_path = r1;
				aristas_a_procesar[Num_Aristas_Procesar] = i;
				printf(" aristas_a_procesar[%d]=%d  ", Num_Aristas_Procesar, aristas_a_procesar[Num_Aristas_Procesar]);
				Num_Aristas_Procesar++;
				r1 = r1 + 0.5;
			}//if
			printf("   \n");
		}//for

		count = pow(2, Num_Aristas_Procesar);
		printf("Tendremos %d recorridos\n", count);
		//inicializas el primer renglon con ceros
		printf(" i=%d.    ", 0);
		for (int j = 0; j < Num_Aristas_Procesar; j++) {
			sentido_de_recorridos[0][j] = 0;
			printf("  %d,   ", sentido_de_recorridos[0][j]);
		}
		printf("   \n");
		//llenas tabla
		for (int i = 1; i < count; i++) {
			printf(" i=%d.    ", i);
			for (int j = 0; j < Num_Aristas_Procesar; j++) {
				index = pow(2, j);
				if (i % index == 0) sentido_de_recorridos[i][j] = !sentido_de_recorridos[i - 1][j];
				else  sentido_de_recorridos[i][j] = sentido_de_recorridos[i - 1][j];
				printf(" (%d)   ", sentido_de_recorridos[i][j]);
			}
			printf("   \n");
		}
		recorido_a_realizar = 0; //para contabilizar recorridos, primer renglon para todas las aristas
		//checaCircles = 1;
		break;
	case 'm':
		aristasBGcollision = 0; checaCircles = 0; NumEdgeInterse = 0;
		crossing_black_black = 0;
		//EdgeIntersecNum = 0; EdgeIntersecNumAro = 0;
		printf("Incrementa Capa Actual \n ");
		if (CapaActual < (MaxNivel - 2)) CapaActual++;
		printf("CapaActual=%d \n", CapaActual);
		GiroMin = 0; MinNumColisions = 100;
		GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
		GeneraLineasEnCapa(CapaActual + 1);
		break;
	case 'n':
		pointCircle.clear(); edg_angles.clear();
		aristasBGcollision = 0; checaCircles = 0; NumEdgeInterse = 0;
		crossing_black_black = 0;
		//EdgeIntersecNum = 0; EdgeIntersecNumAro = 0;
		printf("Decrementa Capa Actual \n ");
		if (CapaActual > 0) CapaActual--;
		printf("CapaActual=%d \n", CapaActual);
		GiroMin = 0; MinNumColisions = 100;
		GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
		GeneraLineasEnCapa(CapaActual + 1);
		break;
	case 'h':
		display_nodo_label = !display_nodo_label;
		break;
	case 'a': //Incrementa CapaActual
		despliega_par = !despliega_par;
		CapaActual = 0;
		printf("CapaActual=%d\n", CapaActual);
		break;
	case 'b':  //gamma 
		//GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
		display_gamma = !display_gamma;
		break;
	case 'c':
		printf("\n CapaMasProfunda=%d, CapaMenosProfunda=%d, MaxNivel=%d\n", CapaActual + 1, CapaActual, MaxNivel);
		printf("Creacion de aristas ---------------------------------------------------------------------------\n");
		//creacion de aristas tipo gamma
		GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
		//creacion de aristas tipo beta
		GeneraLineasEnCapa(CapaActual + 1);
		break;
	case 'd': //gamma vs gamma
		ColisionEntreLineasDeDiferenteCapa();
		break;
	case 'e':
		printf("************************<e>\n");
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) edg[EdgeInterse[i]].collision = 0;
		for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) edg[EdgeInter[i]].collision = 0;
		ColisionEntreLineasDeDiferenteCapa();
		detecta_colision_con_circulos_interiores();
		printf("[%d]. NumIntersec (%d) +  NumIntersecAro (%d) = %d \n", Giro, NumIntersec, NumIntersecAro, NumIntersec + NumIntersecAro);
		//_getch();
		break;
	case 'f': //gamma vs beta
		printf("\n"); distance = 0.0;
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
			r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
			printf("edge[%d] vs circle(%f) \n", EdgeInterse[i], r1);
			distance = Circle_Edge_Cross_v1(EdgeInterse[i], r1);
			if (edg[EdgeInterse[i]].collision) {
				printf("Seguir checando con circulos de capas: \n");
				for (unsigned int j = CapaActual + 2; j < MaxNivel; j++) {
					r1 = (MaxNivel - j) * 5.0 - 4; //radio menor, capa mas profunda
					printf("     <C%i,%2.2f>", j, r1);
					//Circle_Edge_Cross_v1(EdgeInterse[i], r1);
					if (distance < r1) printf("*");
					printf(", ");
				}
				printf("\n");
			}//arista marcada
		}
		break;
	case '1':
		//printf("************************<1>   IN \n");
		proceso_Todo();
		//save_in_file();
		break;
	case '.':
		printf("Giro= ");
		for (Giro = 0; Giro < GiroMin + 1; Giro++)
		{
			printf("%d,", Giro);
			GiraHorario(CapaActual);
		}
		break;
	case 32:  //--------------------------------------------------------------------------------------
		inicia_antes_de_recorridos();
		//Deteccion de colisiones
		ColisionEntreLineasDeDiferenteCapa();
		//ColisionEntreLineasDeDiferenteCapaVsLineasCapaMasProfunda();
		detecta_colision_con_circulos_interiores();
		printf("[%d]. NumIntersec (%d) +  NumIntersecAro (%d)\n", Giro, NumIntersec, NumIntersecAro);
		break; //
	case '2':
		proceso_calcula_crossings_Todo();
		break;
	case 27:
		pointCircle.clear(); edg_angles.clear();
		exit(0);
		break;
	}
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: ArrowKey
void ArrowKey(int key, int x, int y)
{
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vArista;
	unsigned int collisionOld;
	Vector3D vect;
	switch (key) {
	case GLUT_KEY_RIGHT:
		sentido_rot = 0; //0 counter clockwise
		break;
	case GLUT_KEY_LEFT:
		sentido_rot = 1; //1 clockwise
		break;
	case GLUT_KEY_UP: //RIGHT
		if ((abs(edg[EdgeInterse[NumEdgeInterse]].targetAngle - edg[EdgeInterse[NumEdgeInterse]].originAngle)) > 1.5) {
			edg[EdgeInterse[NumEdgeInterse]].originAngle = edg[EdgeInterse[NumEdgeInterse]].originAngle + 1;
			printf("(GLUT_KEY_UP) originAngle=%f\n", edg[EdgeInterse[NumEdgeInterse]].originAngle);
			r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
			i1 = 0.0;
			incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
			i1 = edg[EdgeInterse[NumEdgeInterse]].originAngle * incr;
			xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
			yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
			edg[EdgeInterse[NumEdgeInterse]].origin.x = r1 * xc;
			edg[EdgeInterse[NumEdgeInterse]].origin.y = r1 * yc;
			edg[EdgeInterse[NumEdgeInterse]].seguimiento[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1] = edg[EdgeInterse[NumEdgeInterse]].origin;
			//calcula direccion del vector de movimiento
			vect.x = edg[EdgeInterse[NumEdgeInterse]].origin.x - vect.x;
			vect.y = edg[EdgeInterse[NumEdgeInterse]].origin.y - vect.y;
			vect.unitize();
			//deteccion de colision con otras aristas
			/*collisionOld = edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1];
			Detecta_Colision_con_Aristas_Gama(EdgeInterse[NumEdgeInterse], vect);
			if (edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1] < collisionOld && edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1>1)
			{
				edg[EdgeInterse[NumEdgeInterse]].num_collisions++;
				printf("Entro al IF. collisionOld=%d, ", collisionOld);
				printf("edg[%d].collisions[%d]=%d, ", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1, edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1]);
				printf("edg[%d].num_collisions=%d \n", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].num_collisions);
			}
			//ajuste originAngle
			if (edg[EdgeInterse[NumEdgeInterse]].originAngle > 360) edg[EdgeInterse[NumEdgeInterse]].originAngle = 0;
			else if (edg[EdgeInterse[NumEdgeInterse]].originAngle < 0) edg[EdgeInterse[NumEdgeInterse]].originAngle = 360;
			if (Calcula_vector_secante(EdgeInterse[NumEdgeInterse]))
			{
				Crea_mas_vectores_seguimiento(EdgeInterse[NumEdgeInterse], r1);
			}
			printf("cont_seguimiento=%d, cont_fin=%d \n", edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento, edg[EdgeInterse[NumEdgeInterse]].cont_fin);
			*/
		}
		break;
	case GLUT_KEY_DOWN: //LEFT
		if ((abs(edg[EdgeInterse[NumEdgeInterse]].targetAngle - edg[EdgeInterse[NumEdgeInterse]].originAngle)) > 1.5) {
			edg[EdgeInterse[NumEdgeInterse]].originAngle = edg[EdgeInterse[NumEdgeInterse]].originAngle - 1;
			printf("(GLUT_KEY_DOWN) originAngle=%f\n", edg[EdgeInterse[NumEdgeInterse]].originAngle);
			r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
			i1 = 0.0;
			incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
			i1 = edg[EdgeInterse[NumEdgeInterse]].originAngle * incr;
			xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
			yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
			edg[EdgeInterse[NumEdgeInterse]].origin.x = r1 * xc;
			edg[EdgeInterse[NumEdgeInterse]].origin.y = r1 * yc;
			edg[EdgeInterse[NumEdgeInterse]].seguimiento[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1] = edg[EdgeInterse[NumEdgeInterse]].origin;
			//calcula direccion del vector de movimiento
			vect.x = edg[EdgeInterse[NumEdgeInterse]].origin.x - vect.x;
			vect.y = edg[EdgeInterse[NumEdgeInterse]].origin.y - vect.y;
			vect.unitize();
			//deteccion de colision con otras aristas
			collisionOld = edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1];
			Detecta_Colision_con_Aristas_Gama(EdgeInterse[NumEdgeInterse], vect);
			if (edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1] < collisionOld && edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1>1)
			{
				edg[EdgeInterse[NumEdgeInterse]].num_collisions++;
				printf("Entro al IF. collisionOld=%d, ", collisionOld);
				printf("edg[%d].collisions[%d]=%d, ", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1, edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1]);
				printf("edg[%d].num_collisions=%d \n", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].num_collisions);
			}
			//ajuste originAngle
			if (edg[EdgeInterse[NumEdgeInterse]].originAngle > 360) edg[EdgeInterse[NumEdgeInterse]].originAngle = 0;
			else if (edg[EdgeInterse[NumEdgeInterse]].originAngle < 0) edg[EdgeInterse[NumEdgeInterse]].originAngle = 360;
			if (Calcula_vector_secante(EdgeInterse[NumEdgeInterse]))
			{
				Crea_mas_vectores_seguimiento(EdgeInterse[NumEdgeInterse], r1);
			}
			printf("cont_seguimiento=%d, cont_fin=%d \n", edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento, edg[EdgeInterse[NumEdgeInterse]].cont_fin);
		}
		break;
	case GLUT_KEY_PAGE_DOWN:
		break;
	case GLUT_KEY_PAGE_UP:
		break;
	}
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: glutMouse
void glutMouse(int button, int state, int x, int y)
{

}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: glutMotion
void glutMotion(int x, int y)
{

}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: mouseWheel
void mouseWheel(int wheel, int direction, int x, int y)
{
	if (direction > 0)
	{
		//printf(" direction > 0 \n");
		//if (transZ > 0.9)
		if (transZ < -1.2)
		{
			printf(" Zoom In, transZ=%2.2f ", transZ);//Acercate
			transZ = transZ + 1.0;
			printf(" y sale con transZ=%2.2f  \n", transZ);
			CameraPos.z = transZ;
			picker.update(CameraPos);
		}

	}
	else
	{
		//if (transZ > -2)
		if ((transZ - 0.05) > (MaxZoom)) //ajuste por los decimales
		{
			//printf(" transZ=%f  > MaxZoom=%f\n", transZ, MaxZoom);
			printf(" Zoom Out, transZ=%2.2f  ", transZ);//Alejate
			transZ = transZ - 1.0;
			printf(" y sale con transZ=%2.2f  \n", transZ);
			CameraPos.z = transZ;
			picker.update(CameraPos);
		}
	}
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::: main
int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	/*GLenum err = glewInit();
	if (GLEW_OK != err) {
		cout << "Hello World! ";
		return 1;
	}*/
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(1000, 20);
	glutCreateWindow("Graphs 20 Outer");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	init();
	glutSpecialFunc(ArrowKey);
	glutMouseFunc(glutMouse);
	glutMotionFunc(glutMotion);
	glutMouseWheelFunc(mouseWheel);
	glutMainLoop();
	return 0;
}