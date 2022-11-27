#pragma once
//#include <math.h> 
int max(float a, float b)
{
	if (a > b) return a;
	return b;
}
int min(int a, int b)
{
	if (a < b) return a;
	return b;
}
bool onSegment(float p1, float p2, float q1, float q2, float r1, float r2)
{
	if (q1 <= max(p1, r1) && q1 >= min(p1, r1) &&
		q2 <= max(p2, r2) && q2 >= min(p2, r2))
		return true;

	return false;
}
int orientation(float p1x, float p1y, float q1x, float q1y, float r1x, float r1y)
{

	float val = (q1y - p1y) * (r1x - q1x) - (q1x - p1x) * (r1y - q1y);
	//if (val == 0) return 0;  // colinear
	//printf("orientation() p1x=%f, p1y=%f, q1x=%f, q1y=%f, r1x=%f, r1y=%f, val=%f\n", p1x, p1y, q1x, q1y, r1x, r1y, val);
	if (abs(val) < 0.00001) return 0;  // colinear

									   //return (abs(val) > 0.0001) ? 1 : 2; // clock or counterclock wise
	if (val < 0) return 1;
	else return 2;
}

bool doIntersect(float p1x, float p1y, float q1x, float q1y, float p2x, float p2y, float q2x, float q2y)
{
	int o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y);
	int o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y);
	int o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y);
	int o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y);

	//printf("doIntersect() o1=%d, o2=%d, o3=%d, o4=%d \n", o1, o2, o3, o4);
	// General case 
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases 
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1 
	if (o1 == 0 && onSegment(p1x, p1y, p2x, p2y, q1x, q1y)) return true;

	// p1, q1 and q2 are colinear and q2 lies on segment p1q1 
	if (o2 == 0 && onSegment(p1x, p1y, q2x, q2y, q1x, q1y)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2 
	if (o3 == 0 && onSegment(p2x, p2y, p1x, p1y, q2x, q2y)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2 
	if (o4 == 0 && onSegment(p2x, p2y, q1x, q1y, q2x, q2y)) return true;

	return false; // Doesn't fall in any of the above cases 
}

bool Intersect(int linea1nodo1, int linea1nodo2, int linea2nodo1, int linea2nodo2)
{
	bool inter = 0;
	if (linea1nodo1 == linea2nodo1) inter = 0;
	else if (linea1nodo1 == linea2nodo2) inter = 0;
	else if (linea1nodo2 == linea2nodo1) inter = 0;
	else if (linea1nodo2 == linea2nodo2) inter = 0;
	else //si no hay vertice en comun, entonces llamo a la funcion
	{
		//printf("Llama a doIntersec \n");
		if (doIntersect(v[linea1nodo1].centro[0], v[linea1nodo1].centro[1],
			v[linea1nodo2].centro[0], v[linea1nodo2].centro[1],
			v[linea2nodo1].centro[0], v[linea2nodo1].centro[1],
			v[linea2nodo2].centro[0], v[linea2nodo2].centro[1])) inter = 1;
	}
	if (inter) return true;
	return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
//bool Intersect_one_vertex(int linea1nodo1, int linea1nodo2, int linea2nodo1, float xc, float yc)
bool Intersect_one_vertex(int linea1nodo1, int linea1nodo2, float xBegin, float yBegin, float xc, float yc)
{
	bool inter = 0;
	if (doIntersect(v[linea1nodo1].centro[0], v[linea1nodo1].centro[1],
		v[linea1nodo2].centro[0], v[linea1nodo2].centro[1],
		xBegin, yBegin,
		xc, yc)) inter = 1;

	if (inter) return true;
	return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool Intersect_two_vertex(float xBegin1, float yBegin1, float xc1, float yc1, float xBegin, float yBegin, float xc, float yc)
{
	bool inter = 0;
	if (doIntersect(xBegin1, yBegin1, 	xc1, yc1, xBegin, yBegin, xc, yc)) inter = 1;

	if (inter) return true;
	return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GeneraAngulos_para_Tangente_Circumf(unsigned NumAr)
{
	float angulo1, angulo2, angulo;
	Vector3D vect;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	float xc, yc, i1, incr, xcR, ycR;

	//Primero calculo el vector al punto  edg[NumAr].vertex1
	vect.x = v[edg[NumAr].vertex1].centro[0];
	vect.y = v[edg[NumAr].vertex1].centro[1];
	vect.z = 0;
	vect.unitize();
	angulo = vect.dotproduct(Vector3D(1, 0, 0));
	angulo = acos(angulo) * 180 / 3.1415;
	//ajuste para cuadradntes II, IV
	if (vect.x >= 0) {
		if (vect.y >= 0) {
			//printf("Cuadrante I, AngleRight= %f grados \n", angulo);
		}
		else {
			//printf("Cuadrante IV, AngleRight= %f grados  \n", 360 - angulo);
			angulo = 360 - angulo;
		}
	}
	else {
		//printf("Cuadrantes II o III - ");
		if (vect.y >= 0) {
			//printf("Cuadrante II, AngleRight= %f grados  \n", angulo);
		}
		else {
			//printf("Cuadrante III, AngleRight= %f grados  \n", 360 - angulo);
			angulo = 360 - angulo;
		}
	}
	//edg[NumAr].AngleRight = angulo;
	//Segundo, calculo el vector normal a la circunferencia, es decir, al punto edg[NumAr].vertex1
	angulo1 = angulo - 90;
	if (angulo1 > 0) edg[NumAr].AngleRight = angulo1;
	else edg[NumAr].AngleRight = angulo1 + 360;

	angulo2 = angulo + 90;
	if (angulo2 > 360) edg[NumAr].AngleLeft = angulo2 - 360;
	else edg[NumAr].AngleLeft = angulo2;
	//printf("edg[%d].AngleLeft=%f, edg[%d].AngleRight=%f \n\n", NumAr, edg[NumAr].AngleLeft, NumAr, edg[NumAr].AngleRight);
	
	//contador que permite dirigir el vector de la arista
	i1 = 0.0;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumAr].AngleLeft*incr;
	xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
	yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);
	//Right
	i1 = edg[NumAr].AngleRight*incr;
	xcR = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
	ycR = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);

	//dot product de la arista con el vectores Left, Right
	angulo = edg[NumAr].vect.dotproduct(Vector3D(xc, yc, 0));
	//edg[NumAr].dir = trunc( acos(angulo) * 180 / 3.1415 );
	//printf("vect.dotproduct(Left)= %f, equivalente a %f grados, dir=%d\n", angulo, acos(angulo) * 180 / 3.1415, edg[NumAr].dir);
	angulo = edg[NumAr].vect.dotproduct(Vector3D(xcR, ycR, 0));
	//printf("vect.dotproduct(Right)= %f, equivalente a %f grados\n", angulo, acos(angulo) * 180 / 3.1415);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Calcula_vector_Arista(unsigned NumAr)
{
	//printf("= = = Calcula_vector_Arista( ) \n");
	//printf(" edg[%d]=<v[%d]=(%2.2f, %2.2f), v[%d]=(%2.2f, %2.2f)>    \n ", NumAr, edg[NumAr].vertex1, v[edg[NumAr].vertex1].centro[0], v[edg[NumAr].vertex1].centro[1], edg[NumAr].vertex2, v[edg[NumAr].vertex2].centro[0], v[edg[NumAr].vertex2].centro[1]);
	float angulo;
	edg[NumAr].vect.x = v[edg[NumAr].vertex2].centro[0] - v[edg[NumAr].vertex1].centro[0];
	edg[NumAr].vect.y = v[edg[NumAr].vertex2].centro[1] - v[edg[NumAr].vertex1].centro[1];
	edg[NumAr].vect.z = 0;
	edg[NumAr].vect.unitize();
	angulo = edg[NumAr].vect.dotproduct(Vector3D(1, 0, 0));
	edg[NumAr].vectAngle = acos(angulo) * 180 / 3.1415;
	//ajuste para cuadradntes II, IV
	if (edg[NumAr].vect.x >= 0) {
		if (edg[NumAr].vect.y >= 0) {
			//printf("Cuadrante I, angulo= %f grados \n", edg[NumAr].vectAngle);
		}
		else {
			//printf("Cuadrante IV, angulo= %f grados  \n", 360 - edg[NumAr].vectAngle);
			edg[NumAr].vectAngle = 360 - edg[NumAr].vectAngle;
		}
	}
	else {
		//printf("Cuadrantes II o III - ");
		if (edg[NumAr].vect.y >= 0) {
			//printf("Cuadrante II, angulo= %f grados  \n", edg[NumAr].vectAngle);
		}
		else {
			//printf("Cuadrante III, angulo= %f grados  \n", 360 - edg[NumAr].vectAngle);
			edg[NumAr].vectAngle = 360 - edg[NumAr].vectAngle;
		}
	}
	GeneraAngulos_para_Tangente_Circumf(NumAr);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Crea_mas_vectores_seguimiento_extended_path(unsigned NumAr, float radioMenor)
{
	float xc, yc, i1, incr, r2;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D originExt;
	printf("........................Creando vectores edg[%d].seguimiento[ ] ........................\n", NumAr);
	//printf("r1=%f \n", radioMenor);
	i1 = 0.0; r2 = radioMenor + 5;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumAr].originAngle*incr;
	xc = plane1[0] * 1 * cos(i1*PI) + (double) plane2[0] * 1 * sin(i1*PI);
	yc = plane1[1] * 1 * cos(i1*PI) + (double) plane2[1] * 1 * sin(i1*PI);
	//extension de origin
	originExt.x = ((radioMenor + 1) + edg[NumAr].radio_path)*xc;
	originExt.y = ((radioMenor + 1) + edg[NumAr].radio_path)*yc;
	//nuevo, sale de stop[]
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento] = edg[NumAr].origin;
	edg[NumAr].stop[edg[NumAr].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
	//anterior
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
	//colision global
	edg[NumAr].num_collisions = edg[NumAr].num_collisions + edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
	//informa de las colisiones en el segmento
	printf("Crea_vectores--------  edg[%d].collisions[%d]=%d, ", NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
	printf("edg[%d].num_collisions=%d \n", NumAr, edg[NumAr].num_collisions);
	//incr num de vectores de seguimiento
	edg[NumAr].cont_seguimiento++;
	//inicializa
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;

	//verifica las otras aristas que te estan siguiendo
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		if (edg[EdgeInterse[i]].ext_path && EdgeInterse[i] != NumAr) {
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento] = edg[NumAr].origin;
			printf("edg[%d].cont_seguimiento= %d, edg[%d].radio_path=%f, suma radio=%f  \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento, EdgeInterse[i], edg[EdgeInterse[i]].radio_path, (radioMenor + 1) + edg[EdgeInterse[i]].radio_path);
			
			//incrementas un poco mas para el camino de la siguiente arista que te sigue
			originExt.x = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*xc;
			originExt.y = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*yc;
			//
			edg[EdgeInterse[i]].stop[edg[EdgeInterse[i]].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].stop[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//anterior
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].seguimiento[%d], edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento-1, EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//incr num de vectores de seguimiento
			edg[EdgeInterse[i]].cont_seguimiento++;
		}//if ext_path
	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Crea_mas_vectores_seguimiento_new_proc(unsigned NumAr, float radioMenor) {
	float xc, yc, i1, incr, r2;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D originExt;

	printf("Creando vectores de seguimiento[ ] (new_proc)........................\n");
	printf("r1=%f \n", radioMenor);
	i1 = 0.0; r2 = radioMenor + 5;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumAr].originAngle*incr;
	xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
	yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);
	//extension de origin
	originExt.x = ((radioMenor + 1) + edg[NumAr].radio_path)*xc;
	originExt.y = ((radioMenor + 1) + edg[NumAr].radio_path)*yc;
	//nuevo, sale de stop[]
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento] = edg[NumAr].origin;
	edg[NumAr].stop[edg[NumAr].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
	//anterior
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
	//colision global
	edg[NumAr].num_collisions = edg[NumAr].num_collisions + edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
	//informa de las colisiones en el segmento
	//printf("Crea_vectores--------  edg[%d].collisions[%d]=%d, ", NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
	//printf("edg[%d].num_collisions=%d \n", NumAr, edg[NumAr].num_collisions);
	//incr num de vectores de seguimiento
	edg[NumAr].cont_seguimiento++;
	//inicializa
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;

	//verifica las otras aristas que te estan siguiendo
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		if (edg[EdgeInterse[i]].ext_path && edg[EdgeInterse[i]].new_proc && EdgeInterse[i] != NumAr && !edg[EdgeInterse[i]].end) {
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento] = edg[NumAr].origin;
			printf("edg[%d].cont_seguimiento= %d, edg[%d].radio_path=%f, suma radio=%f  \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento, EdgeInterse[i], edg[EdgeInterse[i]].radio_path, (radioMenor + 1) + edg[EdgeInterse[i]].radio_path);

			//incrementas un poco mas para el camino de la siguiente arista que te sigue
			originExt.x = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*xc;
			originExt.y = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*yc;
			//
			edg[EdgeInterse[i]].stop[edg[EdgeInterse[i]].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].stop[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//anterior
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].seguimiento[%d], edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento - 1, EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//incr num de vectores de seguimiento
			edg[EdgeInterse[i]].cont_seguimiento++;
		}//if ext_path
	}//for
}
// ::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Crea_mas_vectores_seguimiento_extended_path_new_proc(unsigned NumAr, float radioMenor) {
	float xc, yc, i1, incr, r2;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D originExt;

	printf("Creando vectores de seguimiento[ ] (new_proc)........................\n");
	printf("r1=%f \n", radioMenor);
	i1 = 0.0; r2 = radioMenor + 5;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumAr].originAngle*incr;
	xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
	yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);
	//extension de origin
	originExt.x = ((radioMenor + 1) + edg[NumAr].radio_path)*xc;
	originExt.y = ((radioMenor + 1) + edg[NumAr].radio_path)*yc;
	//nuevo, sale de stop[]
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento] = edg[NumAr].origin;
	edg[NumAr].stop[edg[NumAr].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
	//anterior
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
	//colision global
	edg[NumAr].num_collisions = edg[NumAr].num_collisions + edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
	//informa de las colisiones en el segmento
	//printf("Crea_vectores--------  edg[%d].collisions[%d]=%d, ", NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
	//printf("edg[%d].num_collisions=%d \n", NumAr, edg[NumAr].num_collisions);
	//incr num de vectores de seguimiento
	edg[NumAr].cont_seguimiento++;
	//inicializa
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;

	//verifica las otras aristas que te estan siguiendo
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		if (edg[EdgeInterse[i]].ext_path && edg[EdgeInterse[i]].new_proc && EdgeInterse[i] != NumAr && !edg[EdgeInterse[i]].end) {
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento] = edg[NumAr].origin;
			printf("edg[%d].cont_seguimiento= %d, edg[%d].radio_path=%f, suma radio=%f  \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento, EdgeInterse[i], edg[EdgeInterse[i]].radio_path, (radioMenor + 1) + edg[EdgeInterse[i]].radio_path);

			//incrementas un poco mas para el camino de la siguiente arista que te sigue
			originExt.x = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*xc;
			originExt.y = ((radioMenor + 1) + edg[EdgeInterse[i]].radio_path)*yc;
			//
			edg[EdgeInterse[i]].stop[edg[EdgeInterse[i]].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].stop[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//anterior
			edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
			printf("edg[%d].seguimiento[%d], edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento - 1, EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
			//incr num de vectores de seguimiento
			edg[EdgeInterse[i]].cont_seguimiento++;
		}//if ext_path
	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Crea_mas_vectores_seguimiento(unsigned NumAr, float radioMenor)
{
	float xc, yc, i1, incr, r2;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D originExt;
	unsigned GrabaCollisions = 0;

	printf("Creando vectores de seguimiento[ ] ........................\n");
	printf("r1=%f \n", radioMenor);
	i1 = 0.0; r2 = radioMenor + 5;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumAr].originAngle*incr;
	xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1.0 * sin(i1*PI);
	yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1.0 * sin(i1*PI);
	//extension de origin
	originExt.x = (radioMenor+1+edg[NumAr].radio_path)*xc;
	originExt.y = (radioMenor+1+ edg[NumAr].radio_path)*yc;
	//nuevo, sale de stop[]
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento] = edg[NumAr].origin;
	edg[NumAr].stop[edg[NumAr].cont_seguimiento] = Vector3D(originExt.x, originExt.y, 0);
	//anterior
	edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1] = Vector3D(originExt.x, originExt.y, 0);
	//colision global
	edg[NumAr].num_collisions = edg[NumAr].num_collisions + edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
	GrabaCollisions = edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
	//informa de las colisiones en el segmento
	printf("Crea_vectores--------  edg[%d].collisions[%d]=%d, ", NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
	printf("edg[%d].num_collisions=%d, edg[%d].radio_path=%f \n", NumAr, edg[NumAr].num_collisions, NumAr, edg[NumAr].radio_path);
	//incr num de vectores de seguimiento
	edg[NumAr].cont_seguimiento++;
	//inicializa
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool Calcula_vector_secante(unsigned NumAr)
{
	//printf("- - - bool Calcula_vector_secante( . . . ) - - - \n");
	Vector3D vect1; //origin, pero al revés por que temina en origen
	Vector3D vect2; //desde source hasta origin
	double angulo;

	vect1.x = edg[NumAr].origin.x;
	vect1.y = edg[NumAr].origin.y;
	vect1.z = 0; vect1.unitize();
	//vect2.x = v[edg[NumAr].vertex1].centro[0]- edg[NumAr].origin.x;
	//vect2.y = v[edg[NumAr].vertex1].centro[1]- edg[NumAr].origin.y;
	vect2.x = edg[NumAr].stop[edg[NumAr].cont_seguimiento-1].x - edg[NumAr].origin.x;
	vect2.y = edg[NumAr].stop[edg[NumAr].cont_seguimiento-1].y - edg[NumAr].origin.y;
	vect2.z = 0; vect2.unitize();
	angulo = vect1.dotproduct(vect2);
	//printf("angulo= %f ,   ", angulo);
	angulo = acos(angulo) * 180.0 / 3.1415;
	//printf(" %f    ", angulo);
	if (angulo > 90) return true;
	else return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Calcula_vector_origin_extended_path(unsigned NumAr)
{
	Vector3D vect;
	float angulo, r1;
	printf("void Calcula_vector_origin(%d) Extended Path. . . . . . . . \n", NumAr);
	vect.z = 0;
	if (edg[NumAr].cont_seguimiento <= 0) edg[NumAr].cont_seguimiento=1;
	printf("edg[%d].seguimiento[%d].x=%2.2f\n", NumAr, edg[NumAr].cont_seguimiento - 1);
	printf("edg[%d].seguimiento[%d].y=%2.2f\n", NumAr, edg[NumAr].cont_seguimiento - 1);
	vect.x = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x;
	vect.y = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y;
	printf("edg[%d].seguimiento[%d].x=%2.2f\n", NumAr, edg[NumAr].cont_seguimiento - 1, vect.x);
	printf("edg[%d].seguimiento[%d].y=%2.2f\n", NumAr, edg[NumAr].cont_seguimiento - 1, vect.y);
	vect.unitize();
	r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
	edg[NumAr].origin.x = vect.x;
	edg[NumAr].origin.y = vect.y;
	edg[NumAr].origin.z = 0;
	//calculo de originAngle  (unsigned long long) 
	angulo = edg[NumAr].origin.dotproduct(Vector3D(1, 0, 0));
	edg[NumAr].originAngle =acos(angulo) * 180 / 3.1415;
	//ajuste para cuadrantes II, IV
	if (edg[NumAr].origin.x >= 0) {
		if (edg[NumAr].origin.y >= 0) {}
		else { edg[NumAr].originAngle = 360 - edg[NumAr].originAngle; }
	}
	else {
		if (edg[NumAr].origin.y >= 0) {}
		else { edg[NumAr].originAngle = 360 - edg[NumAr].originAngle; }
	}
	edg[NumAr].origin.x = r1*vect.x;
	edg[NumAr].origin.y = r1*vect.y;
	printf("originAngle= %f grados, ", edg[NumAr].originAngle);
	//calculo angulo del origen a target
	vect.z = 0;
	vect.x = v[edg[NumAr].vertex2].centro[0];
	vect.y = v[edg[NumAr].vertex2].centro[1];
	vect.unitize();
	angulo = vect.dotproduct(Vector3D(1, 0, 0));
	angulo = acos(angulo) * 180 / 3.1415;

	if (vect.x >= 0) {
		if (vect.y >= 0) {}
		else { angulo = 360 - angulo; }
	}
	else {
		if (vect.y >= 0) {}
		else { angulo = 360 - angulo; }
	}
	printf("Target Angle= %f grados. dif=%f  \n", angulo, abs(angulo - edg[NumAr].originAngle));
	edg[NumAr].targetAngle = angulo;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Calcula_vector_origin(unsigned NumAr)
{
	Vector3D vect;
	float angulo, r1;
	//printf("void Calcula_vector_origin(%d). . . . . . . . \n", NumAr);
	vect.z = 0;
	//vect.x = (v[edg[NumAr].gammav1].centro[0] + v[edg[NumAr].gammav2].centro[0]) / 2;
	//vect.y = (v[edg[NumAr].gammav1].centro[1] + v[edg[NumAr].gammav2].centro[1]) / 2;
	vect.x = v[edg[NumAr].vertex1].centro[0];
	vect.y = v[edg[NumAr].vertex1].centro[1];
	vect.unitize();
	r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
	edg[NumAr].origin.x = vect.x;
	edg[NumAr].origin.y = vect.y;
	edg[NumAr].origin.z = 0;
	//calculo de originAngle
	angulo = edg[NumAr].origin.dotproduct(Vector3D(1, 0, 0));
	edg[NumAr].originAngle = (double) acos(angulo) * 180 / 3.1415;
	//ajuste para cuadrantes II, IV
	if (edg[NumAr].origin.x >= 0) {
		if (edg[NumAr].origin.y >= 0) {	}
		else { edg[NumAr].originAngle = 360 - edg[NumAr].originAngle; }
	}
	else { 
		if (edg[NumAr].origin.y >= 0) { }
		else { edg[NumAr].originAngle = 360 - edg[NumAr].originAngle; }
	}
	//printf("edg[%d].origin.x = %f*%f\n",NumAr, r1, vect.x);
	//printf("edg[%d].origin.y = %f*%f\n", NumAr, r1, vect.y);
	edg[NumAr].origin.x = r1*vect.x;
	edg[NumAr].origin.y = r1*vect.y;
	edg[NumAr].seguimiento[0] = edg[NumAr].origin;
	edg[NumAr].cont_seguimiento = 1;
	edg[NumAr].num_seguidores = 1;
	edg[NumAr].stop[0] = Vector3D(v[edg[NumAr].vertex1].centro[0], v[edg[NumAr].vertex1].centro[1], 0);

	//printf("originAngle= %f grados, ", edg[NumAr].originAngle);
	//calculo angulo del origen a target
	vect.z = 0;
	vect.x = v[edg[NumAr].vertex2].centro[0];
	vect.y = v[edg[NumAr].vertex2].centro[1];
	vect.unitize();
	r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
	angulo = vect.dotproduct(Vector3D(1, 0, 0));
	angulo = acos(angulo) * 180.0 / 3.1415;
	//ajuste para cuadrantes II, IV
	if (vect.x >= 0) {
		if (vect.y >= 0) {}
		else { angulo = 360 - angulo; }
	}
	else {
		if (vect.y >= 0) {}
		else { angulo = 360 - angulo; }
	}
	//printf("Target Angle= %f grados. dif=%f  \n", angulo, abs(angulo- edg[NumAr].originAngle));
	edg[NumAr].targetAngle = angulo; edg[NumAr].collisions[0] = 0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Checa_Capa1(int CapaMenosProfunda) {
	int nodo;
	printf("*   *   *   * CapaMenosProfunda=%d  (size=%d)\n", CapaMenosProfunda, PorNivel[CapaMenosProfunda].size());
	//EdgeIntersecNum = 0;
	for (int j = 0; j < PorNivel[CapaMenosProfunda].size(); j++) {
		nodo = v[PorNivel[CapaMenosProfunda][j]].id;
		printf("v[%d]=>", v[nodo].id);
		//if (CapaMenosProfunda < 2)
		for (int i = 0; i < v[nodo].liga.size(); i++) {
			printf("%d (%d,%d),", i, edg[v[nodo].liga[i]].vertex1, edg[v[nodo].liga[i]].vertex2);
			/*if (v[edg[v[nodo].liga[i]].vertex1].level != v[edg[v[nodo].liga[i]].vertex2].level)
				if (v[edg[v[nodo].liga[i]].vertex1].level == CapaMasProfunda + 1 || v[edg[v[nodo].liga[i]].vertex2].level == CapaMasProfunda + 1) {
					EdgeInterse[EdgeIntersecNum] = v[nodo].liga[i];
					Calcula_vector_Arista(v[nodo].liga[i]);
					EdgeIntersecNum++;
				}*/

		}
		printf("\n");
	}
	_getch();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GeneraLineasEntreCapas(int CapaMasProfunda, int CapaMenosProfunda)
{
	double angulo, result;
	int nodo;
	//printf("CapaMasProfunda=%d, CapaMenosProfunda=%d  (size=%d)\n", CapaMasProfunda, CapaMenosProfunda, PorNivel[CapaMenosProfunda].size());
	EdgeIntersecNum = 0;

	for (std::vector<int>::iterator it = PorNivel[CapaMenosProfunda].begin(); it != PorNivel[CapaMenosProfunda].end(); ++it)
		for (int i : v[*it].liga) {
			//if (CapaMenosProfunda == 4) 	printf("%d (%d,%d),", i, edg[i].vertex1, edg[i].vertex2);


			if (v[edg[i].vertex1].level != v[edg[i].vertex2].level)
			{
				if (v[edg[i].vertex1].level == CapaMasProfunda + 1 || v[edg[i].vertex2].level == CapaMasProfunda + 1)
				{
					EdgeInterse[EdgeIntersecNum] = i;
					//printf("EdgeInterse[%d]=%d,   <%d, %d>\n", EdgeIntersecNum, EdgeInterse[EdgeIntersecNum], edg[EdgeInterse[EdgeIntersecNum]].vertex1, edg[EdgeInterse[EdgeIntersecNum]].vertex2);
					//value aun no lo uso, asi que puedo eliminar estas 3 lineas
					v[edg[i].vertex1].value = v[edg[i].vertex1].level - CapaActual - 1;
					v[edg[i].vertex2].value = v[edg[i].vertex2].level - CapaActual - 1;
					//printf("v[%d].value=%d, v[%d].value=%d\n", edg[i].vertex1, v[edg[i].vertex1].value, edg[i].vertex2, v[edg[i].vertex2].value);
					Calcula_vector_Arista(i);
					EdgeIntersecNum++;
					
				}//if
				
			}//if
			//if (CapaMenosProfunda == 4) 	printf("*%d (%d,%d),", i, edg[i].vertex1, edg[i].vertex2);
		}//for
	printf("\n Fueron %d Lineas entre las capas %d y %d\n\n", EdgeIntersecNum, CapaMenosProfunda, CapaMasProfunda);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool InEdgeInter(int newEdge) //busca repetidos
{
	bool flag = 0;
	for(unsigned int i=0; i<EdgeIntersecNumAro; i++)
		if (EdgeInter[i] == newEdge) { flag = 1; break; }

	return flag;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GeneraLineasEnCapa(unsigned int CapaMasProfunda)
{
	int count = 0;
	EdgeIntersecNumAro = 0;
	printf("CapaMasProfunda=%d, CapaMenosProfunda=%d  (size=%d)\n", CapaMasProfunda, CapaMasProfunda-1, PorNivel[CapaMasProfunda].size());
	for (std::vector<int>::iterator it = PorNivel[CapaMasProfunda].begin(); it != PorNivel[CapaMasProfunda].end(); ++it)
		for (int i : v[*it].liga) {
			//if (CapaMasProfunda == 4)  printf("%d (%d,%d),", i, edg[i].vertex1, edg[i].vertex2);
			
			if (v[edg[i].vertex1].level == v[edg[i].vertex2].level)
				//if (v[edg[i].vertex1].level == CapaMasProfunda+1 && v[edg[i].vertex2].level == CapaMasProfunda+1)
			{
				//printf("(%d)edg[%d]=<%d,%d> ", *it, i, edg[i].vertex1, edg[i].vertex2);
				/*if (count > 0)
				{
					if(EdgeInter[count] < i) EdgeInter[count] = i;
				}
				else
				 EdgeInter[count] = i;*/
				if (!InEdgeInter(i)) {
					//printf("!"); 
					EdgeInter[count] = i;
					//printf("EdgeInter[%d]=%d<%d, %d>, collision=%d, ", count, EdgeInter[count], edg[EdgeInter[count]].vertex1, edg[EdgeInter[count]].vertex2, edg[EdgeInter[count]].collision);
					count++; EdgeIntersecNumAro = count;
				}

			}
		}
	//*printf("Aristas en aro interno (Beta)=====\n");
	//*printf("Fueron %d Lineas en el anillo de la capa %d (mas profunda)\n\n", count, CapaMasProfunda);
	EdgeIntersecNumAro = count;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Detecta_Colision_con_Aristas_Gama_Extended_Path(unsigned int NumAr)
{
	Vector3D vect, vect1;
	float r1;

	//printf("Detecta_Colision_con_Aristas_Gama (arista %d) \n", NumAr);
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;
	r1 = (MaxNivel - CapaActual - 1)*5.0; //radio intermedio
	for (unsigned int i = 0; i < EdgeIntersecNum; i++)
	{
		vect.x = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x; 
		vect.y = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y;
		vect.z = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].z;
		vect.unitize();
		vect1.x = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x; 
		vect1.y = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y;
		vect1.z = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].z;
		vect1.unitize();
		//printf("edge %d, ", EdgeInterse[i]);
		if (EdgeInterse[i] != NumAr) 
			if (edg[EdgeInterse[i]].collision) {
				//printf("Red ");
				if (edg[EdgeInterse[i]].ext_path && edg[EdgeInterse[i]].cont_fin==0) { //seguidora existente con ext_path
					
					report[i].rotations++;
					printf("Seguidora Ext Path, cont_seguimiento= %d, rotations=%d", edg[EdgeInterse[i]].cont_seguimiento, report[i].rotations);
					edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1].x = edg[NumAr].origin.x; //solo actualizas
					edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1].y = edg[NumAr].origin.y;
					edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1].z = edg[NumAr].origin.z;

					//printf("edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento - 1);
					printf("follower edg[%d], cont_seguimiento=%d, angulos %f-%f, Follower edg[%d] edg[EdgeInterse[i]].radio_path=%f \n", EdgeInterse[i], edg[NumAr].cont_seguimiento, edg[EdgeInterse[i]].targetAngle, edg[NumAr].originAngle, EdgeInterse[i], edg[EdgeInterse[i]].radio_path);
					if ((abs(edg[EdgeInterse[i]].targetAngle - edg[NumAr].originAngle)) < 1.5) {
						edg[EdgeInterse[i]].cont_fin = edg[EdgeInterse[i]].cont_seguimiento;
						printf("\n FINALIZA seguidora de ext_path edg[%d], con edg[%d].cont_fin=%d ..............................\n", EdgeInterse[i], EdgeInterse[i], edg[EdgeInterse[i]].cont_fin);
						edg[EdgeInterse[i]].cont_seguimiento = 0;
						edg[EdgeInterse[i]].end = 1;
						report[i].tipo = 3; report[i].cont_seguimiento = edg[EdgeInterse[i]].cont_fin;
						//report[i].rotations = report[NumEdgeInterse].rotations - report[i].rotations;
						printf(" ext_path tiene  %d rotaciones y su follower recien terminada tiene %d rotaciones \n ", report[NumEdgeInterse].rotations, report[i].rotations);
					}
				}
				else if (edg[EdgeInterse[i]].cont_seguimiento == 0 && edg[EdgeInterse[i]].cont_fin==0) { //Nueva seguidora
						if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
							r1* vect.x, r1* vect.y, r1* vect1.x, r1* vect1.y)) {
						
						edg[EdgeInterse[i]].seguimiento[0] = edg[NumAr].origin;
						edg[EdgeInterse[i]].cont_seguimiento = 1;
						edg[EdgeInterse[i]].stop[0] = Vector3D(v[edg[EdgeInterse[i]].vertex1].centro[0], v[edg[EdgeInterse[i]].vertex1].centro[1], 0);
						edg[NumAr].num_seguidores++; //secuencia de aristas seguidoras para determinar el radio de sus caminos
						edg[EdgeInterse[i]].radio_path = edg[NumAr].num_seguidores;
						edg[EdgeInterse[i]].ext_path = 1;
						printf("Red collision, Seguidora Nueva edg[%d].radio_path=%f \n ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path);
						printf("Pivote tiene  %d rotaciones \n ", report[NumEdgeInterse].rotations);
						report[i].rotations = report[NumEdgeInterse].rotations;
					}//if Intersec
				}

			}//Red edge
			/*else if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
				r1* vect.x, r1* vect.y, r1* vect1.x, r1* vect1.y)) {
					printf("Black collision");
					edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]++;
			}*///Black Edge

		//printf("\n");
	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Detecta_Colision_con_Aristas_Gama(unsigned int NumAr, Vector3D dir)
{
	Vector3D vect, vect1, vect2;
	float r1, nodo1X, nodo1Y, nodo2X, nodo2Y;
	double dot;
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;
	r1 = (MaxNivel - CapaActual - 1)*5.0-1; //radio intermedio
	for (unsigned int i = 0; i < EdgeIntersecNum; i++)
	{ 
		//para determinar el vector extendido del Intersec
		vect.x = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x; 
		vect.y = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y;
		vect.z = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].z;
		//cout << "i=" << i << ", edg[" << NumAr << "].cont_seguimiento=" << edg[NumAr].cont_seguimiento << "  x=" << vect.x << " y=" << vect.y << " z=" << vect.z << "\n";
		vect1.x = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x;
		vect1.y = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y;
		vect1.z = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].z;

		vect2.x = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].x - edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x;
		vect2.y = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1].y - edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y;
		vect.unitize(); vect1.unitize(); vect2.unitize();
		nodo1X = r1* vect.x;
		nodo1Y = r1* vect.y;
		nodo2X = r1* vect1.x;
		nodo2Y = r1* vect1.y;
		if (arranca) {
			dot = vect2.dotproduct(dir);
			if (dot < 0.08) {
				nodo1X = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].x;
				nodo1Y = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1].y;
				nodo2X = edg[NumAr].origin.x;
				nodo2Y = edg[NumAr].origin.y;
			}
			else arranca = 0;
		}//arranca if
		else {
		}//arranca else
		//
		if (EdgeInterse[i] == NumAr) printf("");
		else {
			if (edg[EdgeInterse[i]].cont_seguimiento > 0) //follower
			{//creado, hubo colision en su momento
				edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento-1] = edg[NumAr].origin; //solo actualizas
				//printf("follower edg[%d], cont_seguimiento=%d, angulos %f-%f \n", EdgeInterse[i], edg[NumAr].cont_seguimiento, edg[EdgeInterse[i]].targetAngle, edg[NumAr].originAngle);

				if ((abs(edg[EdgeInterse[i]].targetAngle - edg[NumAr].originAngle)) < 1.5) {
					printf("\n FINALIZA seguidora edg[%d] ......................cuando pivote tenia  %d, %d colisiones, segm %d \n", EdgeInterse[i], edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1], edg[NumAr].num_collisions, edg[NumAr].cont_seguimiento - 1);
					//
					report[i].tipo = 3;
					report[i].num_collisions = edg[EdgeInterse[i]].num_collisions;
					report[i].cont_seguimiento = edg[EdgeInterse[i]].cont_seguimiento;
					report[i].radio_path = edg[EdgeInterse[i]].radio_path;
					report[i].rotations = report[NumEdgeInterse].rotations - report[i].rotations;
					//
					edg[EdgeInterse[i]].cont_fin = edg[EdgeInterse[i]].cont_seguimiento; 
					edg[EdgeInterse[i]].cont_seguimiento = 0;
					edg[EdgeInterse[i]].end = 1;
					printf("edg[%d], ", EdgeInterse[i]);
					printf("num_collisions=%d, ", edg[EdgeInterse[i]].num_collisions);
					printf("cont_seguimiento=%d, ", edg[EdgeInterse[i]].cont_seguimiento);
					printf("cont_fin=%d, end=%d, tipo=%d", edg[EdgeInterse[i]].cont_fin, edg[EdgeInterse[i]].end, report[i].tipo);
					printf("\n");
					if (edg[EdgeInterse[i]].cont_seguimiento - 1 == 0) //es el primer segmento
						edg[EdgeInterse[i]].num_collisions = edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] - edg[EdgeInterse[i]].num_collisions;
					else
						edg[EdgeInterse[i]].num_collisions = edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] + edg[EdgeInterse[i]].num_collisions;

					printf("edg[%d].num_collisions=%d \n", EdgeInterse[i], edg[EdgeInterse[i]].num_collisions);
					printf("Follower de Pivote tiene  %d rotaciones  y %d segmentos\n ", report[NumEdgeInterse].rotations, edg[EdgeInterse[i]].cont_fin);
				}
			}

			else if (edg[EdgeInterse[i]].collision) //aristas rojas 
			{
				if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
					nodo1X, nodo1Y, nodo2X, nodo2Y))
				{
					//checa producto punto del segmento pivote con la arista roja
					vect2.unitize();
					dot = vect2.dotproduct(Vector3D(1, 0, 0)); //segmento de la pivote
					double angulo = acos(dot) * 180 / 3.1415;

					//ajuste para cuadrantes II, IV
					if (vect2.x >= 0) {
						if (vect2.y >= 0) {}
						else { angulo = 360 - angulo; }
					}
					else {
						if (vect2.y >= 0) {}
						else { angulo = 360 - angulo; }
					}
					//el vector de la otra arista
					vect.x = v[edg[EdgeInterse[i]].vertex2].centro[0] - v[edg[EdgeInterse[i]].vertex1].centro[0];
					vect.y = v[edg[EdgeInterse[i]].vertex2].centro[1] - v[edg[EdgeInterse[i]].vertex1].centro[1];
					vect.z = 0;
					vect.unitize();
					dot = vect.dotproduct(Vector3D(1, 0, 0));
					double anguloSeg = acos(dot) * 180 / 3.1415;
					//ajuste para cuadrantes II, IV
					if (vect.x >= 0) {
						if (vect.y >= 0) {}
						else { anguloSeg = 360 - anguloSeg; }
					}
					else {
						if (vect.y >= 0) {}
						else { anguloSeg = 360 - anguloSeg; }
					}
					//
					printf("CHOCA  con arista. angulo segm %f, y angulo otra arista=%f\n", angulo, anguloSeg);
					//producto punto del segm con la otra arista
					dot = vect.dotproduct(vect2);
					anguloSeg = acos(dot) * 180 / 3.1415;
					angulo = anguloSeg;
					//ajuste para cuadrantes II, IV
					if (vect.x >= 0) {
						if (vect.y >= 0) {}
						else { anguloSeg = 360 - anguloSeg; }
					}
					else {
						if (vect.y >= 0) {}
						else { anguloSeg = 360 - anguloSeg; }
					}
					printf("Producto punto %f, y angulo =%f, %f\n", dot, angulo, anguloSeg);
					if (dot > 0.0) {
						//if (edg[EdgeInterse[i]].collision) { 
						printf(" - - - - (RED) crea edg[%d].seguimiento[0] cuando pivote tenia  %d, %d colisiones, segm %d\n ",
							EdgeInterse[i], edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1], edg[NumAr].num_collisions, edg[NumAr].cont_seguimiento - 1);
						printf("Pivote tiene  %d rotaciones \n ", report[NumEdgeInterse].rotations);
						report[i].rotations = report[NumEdgeInterse].rotations;
						edg[EdgeInterse[i]].seguimiento[0] = edg[NumAr].origin;
						edg[EdgeInterse[i]].cont_seguimiento = 1;
						edg[EdgeInterse[i]].stop[0] = Vector3D(v[edg[EdgeInterse[i]].vertex1].centro[0], v[edg[EdgeInterse[i]].vertex1].centro[1], 0);
						edg[EdgeInterse[i]].radio_path = edg[NumAr].radio_path + 0.5 * edg[NumAr].num_seguidores;
						edg[NumAr].num_seguidores++;
						edg[EdgeInterse[i]].num_collisions = edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1];
						printf("Red collision, Seguidora Nueva edg[%d].radio_path=%f \n ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path);
						//} //if collision
						/*else {
							edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]++; //aristas negras
							printf("colision con arista negra, pivote tiene %d seguidoras y  edg[%d].collisions[%d]=%d\n", edg[NumAr].num_seguidores, NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
						}*/
					} //dot > 0
				}//else if Intersect
			}//rojas
			else { //negras
			/*if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
				nodo1X, nodo1Y, nodo2X, nodo2Y)) {
				edg[NumAr].num_collisions++;
				printf("Colision con arista negra i=%d, EdgeInterse[i]=%d, num_collisions[%d]=%d \n", i, EdgeInterse[i], NumAr, edg[NumAr].num_collisions);
			}*/
			}//else del else if Intersec
		}//else
			
	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Detecta_Colision_con_Aristas_Gama_New_Proc(unsigned int NumAr) {
	Vector3D vect, vect1;
	float r1;
	printf("Detecta_Colision_con_Aristas_Gama (arista %d) \n", NumAr);
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;
	r1 = (MaxNivel - CapaActual - 1)*5.0; //radio intermedio
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		vect = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1]; vect.unitize();
		vect1 = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1]; vect1.unitize();

		if (EdgeInterse[i] != NumAr)
			if (edg[EdgeInterse[i]].collision) {
				if (edg[EdgeInterse[i]].new_proc && edg[EdgeInterse[i]].cont_fin == 0 && edg[EdgeInterse[i]].ext_path) { //seguidora existente con new_proc
					printf("Seguidora Ext Path, edg[%d].cont_seguimiento= %d, ", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento);
					edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1] = edg[NumAr].origin; //solo actualizas
					printf("edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento - 1);
					if ((abs(edg[EdgeInterse[i]].targetAngle - edg[NumAr].originAngle)) < 1.5) {
						edg[EdgeInterse[i]].cont_fin = edg[EdgeInterse[i]].cont_seguimiento;
						printf("\n FINALIZA seguidora de New_Proc edg[%d], con edg[%d].cont_fin=%d ..............................\n", EdgeInterse[i], EdgeInterse[i], edg[EdgeInterse[i]].cont_fin);
						
						report[i].cont_seguimiento = edg[EdgeInterse[i]].cont_seguimiento;
						report[i].radio_path = edg[EdgeInterse[i]].radio_path;
						report[i].tipo = 6;
						report[i].rotations = report[NumEdgeInterse].rotations - report[i].rotations;
						edg[EdgeInterse[i]].cont_seguimiento = 0;
						edg[EdgeInterse[i]].end = 1;
						Num_aristas_pendientes--; //decrementa pending edges
					}
				}
				else if (edg[EdgeInterse[i]].cont_seguimiento == 0 && edg[EdgeInterse[i]].cont_fin == 0 && edg[EdgeInterse[i]].new_proc) { //Nueva seguidora
					if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
						r1* vect.x, r1* vect.y, r1* vect1.x, r1* vect1.y)) {
						
						//
						edg[EdgeInterse[i]].seguimiento[0] = edg[NumAr].origin;
						edg[EdgeInterse[i]].cont_seguimiento = 1;
						edg[EdgeInterse[i]].stop[0] = Vector3D(v[edg[EdgeInterse[i]].vertex1].centro[0], v[edg[EdgeInterse[i]].vertex1].centro[1], 0);

						//edg[EdgeInterse[i]].radio_path = edg[NumAr].radio_path + 0.5;
						edg[EdgeInterse[i]].radio_path = edg[NumAr].radio_path + 0.5*edg[NumAr].num_seguidores;
						edg[EdgeInterse[i]].ext_path = 1;
						edg[NumAr].num_seguidores++;
						printf("Red collision, Seguidora Nueva edg[%d].radio_path=%f \n ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path);
						report[i].rotations = report[NumEdgeInterse].rotations;
					}//if Intersec
				}

			}//Red edge
			/*else if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
				r1* vect.x, r1* vect.y, r1* vect1.x, r1* vect1.y)) {
				printf("Black collision");
				edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]++;
				printf("edg[%d].collisions[%d]=%d, ", NumAr, edg[NumAr].cont_seguimiento - 1, edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1]);
				printf("edg[%d].num_collisions=%d \n", NumAr, edg[NumAr].num_collisions);
			}*///Black Edge

			//printf("\n");

	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Detecta_Colision_con_Aristas_Gama_Extended_Path_New_Proc(unsigned int NumAr) {
	Vector3D vect, vect1;
	float r1;
	//printf("Detecta_Colision_con_Aristas_Gama (arista %d) \n", NumAr);
	edg[NumAr].collisions[edg[NumAr].cont_seguimiento - 1] = 0;
	r1 = (MaxNivel - CapaActual - 1)*5.0; //radio intermedio
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		vect = edg[NumAr].stop[edg[NumAr].cont_seguimiento - 1]; vect.unitize();
		vect1 = edg[NumAr].seguimiento[edg[NumAr].cont_seguimiento - 1]; vect1.unitize();
		//printf("edge %d, ", EdgeInterse[i]);
		if (EdgeInterse[i] != NumAr)
			if (edg[EdgeInterse[i]].collision) {
				if (edg[EdgeInterse[i]].new_proc && edg[EdgeInterse[i]].cont_fin == 0 && edg[EdgeInterse[i]].ext_path) { //seguidora existente con new_proc
					printf("Seguidora Ext Path, cont_seguimiento= %d, ", edg[EdgeInterse[i]].cont_seguimiento);
					edg[EdgeInterse[i]].seguimiento[edg[EdgeInterse[i]].cont_seguimiento - 1] = edg[NumAr].origin; //solo actualizas
					printf("edg[%d].seguimiento[%d] \n", EdgeInterse[i], edg[EdgeInterse[i]].cont_seguimiento - 1);
					report[i].rotations++;
					if ((abs(edg[EdgeInterse[i]].targetAngle - edg[NumAr].originAngle)) < 1.5) {
						edg[EdgeInterse[i]].cont_fin = edg[EdgeInterse[i]].cont_seguimiento;
						printf("\n FINALIZA seguidora de New_Proc Ext_Path edg[%d], con edg[%d].cont_fin=%d ..............................\n", EdgeInterse[i], EdgeInterse[i], edg[EdgeInterse[i]].cont_fin);

						report[i].cont_seguimiento = edg[EdgeInterse[i]].cont_seguimiento;
						report[i].radio_path = edg[EdgeInterse[i]].radio_path;
						report[i].tipo = 6; report[i].cont_seguimiento = edg[EdgeInterse[i]].cont_fin;
						//report[i].rotations = report[NumEdgeInterse].rotations - report[i].rotations;
						edg[EdgeInterse[i]].cont_seguimiento = 0;
						edg[EdgeInterse[i]].end = 1;
						Num_aristas_pendientes--; //decrementa pending edges
					}
				}
				else if (edg[EdgeInterse[i]].cont_seguimiento == 0 && edg[EdgeInterse[i]].cont_fin == 0 && edg[EdgeInterse[i]].new_proc) { //Nueva seguidora
					if (Intersect_one_vertex(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2,
						r1* vect.x, r1* vect.y, r1* vect1.x, r1* vect1.y)) {

						edg[EdgeInterse[i]].seguimiento[0] = edg[NumAr].origin;
						edg[EdgeInterse[i]].cont_seguimiento = 1;
						edg[EdgeInterse[i]].stop[0] = Vector3D(v[edg[EdgeInterse[i]].vertex1].centro[0], v[edg[EdgeInterse[i]].vertex1].centro[1], 0);

						//edg[EdgeInterse[i]].radio_path = edg[NumAr].radio_path + 0.5;
						edg[EdgeInterse[i]].radio_path = edg[NumAr].radio_path + 0.5*edg[NumAr].num_seguidores;
						edg[EdgeInterse[i]].ext_path = 1;
						edg[NumAr].num_seguidores++;
						printf("Red collision, New_Proc Seguidora Nueva edg[%d].radio_path=%f \n ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path);
						report[i].rotations = report[NumEdgeInterse].rotations;
					}//if Intersec
				}

			}//Red edge

	}//for
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ColisionEntreLineasDeDiferenteCapa()
{
	NumIntersec = 0; //a ver cuanto pares de lineas se intersectan
	if (EdgeIntersecNum < 2)
	{
		printf("Una o menos ligas \n");
	}//if
	else if (EdgeIntersecNum < 3) //es solo una pareja
	{
		//*printf("Es una sola pareja, EdgeIntersecNum=%d\n", EdgeIntersecNum);
		//if (Intersect(EdgeIntersec[0].nodo1, EdgeIntersec[0].nodo2, EdgeIntersec[1].nodo1, EdgeIntersec[1].nodo1))
		if (Intersect(edg[EdgeInterse[0]].vertex1, edg[EdgeInterse[0]].vertex2, edg[EdgeInterse[1]].vertex1, edg[EdgeInterse[1]].vertex2))
		{
			//edg[EdgeInterse[0]].collision = 1; edg[EdgeInterse[1]].collision = 1;
			NumIntersec = 1;
		}
	} //else if
	else //mas de 3 ligas, se forman parejas
	{
		//printf("Las parejas formadas son: \n");  //EdgeIntersecNum = EdgeIntersecNum - 1;
		for (unsigned int i = 0; i < EdgeIntersecNum - 1; i++)
		{
			//printf("i=%d\n", i);
			for (unsigned int j = i + 1; j <= EdgeIntersecNum - 1; j++)
			{
				//printf("j=%d, ", j);
				if (Intersect(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInterse[j]].vertex1, edg[EdgeInterse[j]].vertex2))
				{
					//--edg[EdgeInterse[i]].collision = 1; edg[EdgeInterse[j]].collision = 1;
					//printf("[%d] nodos <%d, %d> vs [%d] nodos <%d, %d>, collision=%d, %d\n", i, edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, j, edg[EdgeInterse[j]].vertex1, edg[EdgeInterse[j]].vertex2, edg[EdgeInterse[i]].collision, edg[EdgeInterse[j]].collision);
					NumIntersec++;
					report[i].num_collisions_red++; report[j].num_collisions_red++;
					//printf("NumIntersec=%d\n", NumIntersec);
					//si ambas son gamma negra, entonces contabilizo aparte para el tercer reporte
					//if (!edg[EdgeInterse[i]].collision && !edg[EdgeInterse[j]].collision) crossing_black_black++;
				}
			}
			//printf("\n");
		}
		//printf("\n");
	}//else
	//printf("======Tuvimos %d intersecciones Gamma vs Gamma \n", NumIntersec);
	/*if (NumIntersec == 0)
	{
		printf("Checare colisiones con las ligas de la capa profunda \n"); //recuerda que son dos capas que estamos evaluando
	}*/
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ColisionEntreLineasDeDiferenteCapaVsLineasCapaMasProfunda()
{
	Vector3D vect;
	float r1, angulo;

	NumIntersecAro = 0;
	//printf("EdgeIntersecNum=%d vs EdgeIntersecNumAro=%d \n", EdgeIntersecNum, EdgeIntersecNumAro);

	//printf("Las parejas formadas son: \n");  //EdgeIntersecNum = EdgeIntersecNum - 1;
	for (unsigned int i = 0; i < EdgeIntersecNum; i++)
	{
		//printf("i=%d\n", i);
		//for (unsigned int j = i + 1; j <= EdgeIntersecNum; j++)
		edg[EdgeInterse[i]].gammav1 = -1;
		edg[EdgeInterse[i]].gammav2 = -1;
		for (int j = 0; j < EdgeIntersecNumAro; j++)
		{
			//printf("i=%d, j=%d, ", i, j);
			if (Intersect(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInter[j]].vertex1, edg[EdgeInter[j]].vertex2))
			{
				edg[EdgeInterse[i]].collision = 1; edg[EdgeInter[j]].collision = 1;
				NumIntersecAro++;
				report[i].num_collisions++;
				//printf("NumIntersecAro=%d con arista j=%d.  ", NumIntersecAro, j);
				//printf("liga <%d, %d> vs liga Aro <%d, %d>, collision=%d, %d\n", edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInter[j]].vertex1, edg[EdgeInter[j]].vertex2, edg[EdgeInterse[i]].collision, edg[EdgeInter[j]].collision);
				edg[EdgeInterse[i]].gammav1 = edg[EdgeInter[j]].vertex1;
				edg[EdgeInterse[i]].gammav2 = edg[EdgeInter[j]].vertex2;
				//printf("gammav1=%d, gammav2=%d \n", edg[EdgeInterse[i]].gammav1, edg[EdgeInterse[i]].gammav2);

				//printf("originAngle =%f \n", edg[EdgeInterse[i]].originAngle);
			}//if Inersec
		}//for

	}
	//printf("\n");

	printf("====== Tuvimos %d intersecciones (gamma vs beta) \n", NumIntersecAro);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GiraHorario(unsigned int CapaMenosProfunda)
{
	float xtemp, ytemp;
	int count = 0, index = 0;
	float tmpx[2], tmpy[2];
	bool inter = 0;

	//printf("\n Hay %d nodos en el nivel %d\n", NodosPorNivel[CapaMenosProfunda], CapaMenosProfunda);
	for (std::vector<int>::iterator it = PorNivel[CapaMenosProfunda].begin(); it != PorNivel[CapaMenosProfunda].end() - 1; ++it)
	{
		if (count == 0)
		{
			tmpx[1] = v[*it].centro[0]; tmpy[1] = v[*it].centro[1];
			//guarda el indice del primer nodo
			index = *it;
		}

		tmpx[inter] = v[*it + 1].centro[0]; tmpy[inter] = v[*it + 1].centro[1];
		v[*it + 1].centro[0] = tmpx[!inter]; v[*it + 1].centro[1] = tmpy[!inter];
		inter = !inter;
		count++;
	}
	v[index].centro[0] = tmpx[!inter]; v[index].centro[1] = tmpy[!inter];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GiraAntiHorario(unsigned int CapaMenosProfunda)
{
	float xtemp, ytemp;
	int count = 0;
	printf("\n Hay %d nodos en el nivel %d\n", NodosPorNivel[CapaMenosProfunda], CapaMenosProfunda);
	for (std::vector<int>::iterator it = PorNivel[CapaMenosProfunda].begin(); it != PorNivel[CapaMenosProfunda].end(); ++it)
	{
		//printf("count=%d . . . ", count);
		if (count == 0)
		{
			xtemp = v[*it].centro[0]; ytemp = v[*it].centro[1];
			//printf("grabando las coordenadas de v[%d]\n", *it);
		}

		if (count == NodosPorNivel[CapaMenosProfunda] - 1)//es uno menos que 10, es decir, 9
		{
			v[*it].centro[0] = xtemp; v[*it].centro[1] = ytemp;
			//printf("v[%d].x=xtemp, v[%d].y=ytemp\n", *it, *it);
		}
		else
		{
			v[*it].centro[0] = v[*it + 1].centro[0]; v[*it].centro[1] = v[*it + 1].centro[1];
			//printf("v[%d].x=v[%d].x, v[%d].y=v[%d].y\n", *it, *it + 1, *it, *it + 1);
		}
		count++;
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool Verifica(unsigned int CapaMenosProfunda, bool giro, bool coord)
{
	unsigned int ciclo = 0, colisiones[100]; //este tamaño debe cambiar
											 //printf("Verifica: CapaMasProfunda=%d, CapaMenosProfunda=%d\n", CapaMasProfunda, CapaMenosProfunda);
											 //recuerda que aqui habra un ciclo
	while (ciclo < NodosPorNivel[CapaMenosProfunda])
		//while (ciclo < 9)
	{
		printf("ciclo=%d - - - - - -  - -- -  - - - - - - - - - - - - - -\n", ciclo);
		//colisiones[ciclo] = 0;
		ColisionEntreLineasDeDiferenteCapa();
		ColisionEntreLineasDeDiferenteCapaVsLineasCapaMasProfunda();
		//aunque debo guardar los valores de las intersecciones para que en caso 
		//de que no de cero, quede el menor numero de intersecciones
		colisiones[ciclo] = NumIntersec + NumIntersecAro;
		printf("============ colisiones[%d]=%d, ", ciclo, colisiones[ciclo]);
		if (NumIntersec == 0 && NumIntersecAro == 0)
		{
			printf("YA QUEDO \n");
			return 1;
		}
		else
		{
			if (!giro)
			{
				printf("Gira en sentido anti horario \n");
				GiraAntiHorario(CapaMenosProfunda);
			}
			else
			{
				printf("Gira en sentido horario \n");
				GiraHorario(CapaMenosProfunda);
			}
		}
		//getchar();
		ciclo++;
	}//while
	printf("Numero de colisiones por giro ");
	if (!giro) printf("anti horario \n");
	else printf("horario \n");

	//calcula el menor de las colisiones
	ciclo = 0;
	unsigned int menor = 1000;
	for (int i = 0; i < NodosPorNivel[CapaMenosProfunda]; i++)
	{
		printf("colisiones[%d]=%d, ", i, colisiones[i]);
		if (colisiones[i] < menor) { menor = colisiones[i]; ciclo = i; }
	}
	printf("El menor es %d en el ciclo %d \n\n", menor, ciclo);
	if (coord == 0 && giro == 0) { NumColision00 = menor; ciclo00 = ciclo; }
	else if (coord == 0 && giro == 1) { NumColision01 = menor; ciclo01 = ciclo; }
	else if (coord == 1 && giro == 0) { NumColision10 = menor; ciclo10 = ciclo; }
	else { NumColision11 = menor; ciclo11 = ciclo; }

	getchar();
	return 0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void reproduce(unsigned int config, unsigned int ciclo, unsigned int CapaMenosProfunda)
{
	if (config == 0) //llamar nuevamente a AdjustOuter()
	{
		AdjustOuterGraph(0, MaxNivel - 1);
		for (unsigned int i = 0; i<ciclo; i++)     GiraAntiHorario(CapaMenosProfunda);
	}
	else if (config == 1)
	{
		AdjustOuterGraph(0, MaxNivel - 1);
		for (unsigned int i = 0; i<ciclo; i++)     GiraHorario(CapaMenosProfunda);
	}
	else if (config == 2) //ya estan en opposite direction
	{
		//AdjustOuterGraphOpposite(CapaMenosProfunda);
		for (unsigned int i = 0; i<ciclo; i++)     GiraAntiHorario(CapaMenosProfunda);
	}
	else
	{
		//AdjustOuterGraphOpposite(CapaMenosProfunda);
		for (unsigned int i = 0; i<ciclo; i++)     GiraHorario(CapaMenosProfunda);
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Test(unsigned int CapaMasProfunda, unsigned int CapaMenosProfunda)
{
	bool fin = 0;
	printf("Test: CapaMasProfunda=%d, CapaMenosProfunda=%d -  -  -  -  -  -  -  -  -  - -  -  -  -  -  -  -  -  -  - \n", CapaMasProfunda, CapaMenosProfunda);
	if (Verifica(CapaMenosProfunda, 0, 0)) { printf("Proceso terminado ...\n"); fin = 1; }
	else if (Verifica(CapaMenosProfunda, 1, 0)) { printf("Proceso terminado ...\n"); fin = 1; }
	else
	{
		printf("-  -  -  -  -  -  -  -  -  -  AdjustOuterGraphOpposite -  -  -  -  -  -  -  -  -  - ");
		AdjustOuterGraphOpposite(CapaMenosProfunda, 5); //pasalo a booleano y
		if (Verifica(CapaMenosProfunda, 0, 1)) { printf("Proceso terminado ...\n"); fin = 1; }
		else if (Verifica(CapaMenosProfunda, 1, 1)) { printf("Proceso terminado ...\n"); fin = 1; }
	}

	printf("NumColision00=%d, ciclo00=%d \n", NumColision00, ciclo00);
	printf("NumColision01=%d, ciclo01=%d \n", NumColision01, ciclo01);
	printf("NumColision10=%d, ciclo10=%d \n", NumColision10, ciclo10);
	printf("NumColision11=%d, ciclo11=%d \n", NumColision11, ciclo11);
	//config => 0 (00), 1(01), 2(10), 3(11)
	unsigned int  ciclo0, ciclo1, ciclo;
	unsigned int menor0, menor1, menor;
	unsigned int config;
	if (!fin)
	{
		if (NumColision00 < NumColision01) { ciclo0 = ciclo00; menor0 = NumColision00; }
		else { ciclo0 = ciclo01; menor0 = NumColision01; }

		if (NumColision10 < NumColision11) { ciclo1 = ciclo10; menor1 = NumColision10; }
		else { ciclo1 = ciclo11; menor1 = NumColision11; }

		if (menor0 < menor1)
		{
			menor = menor0; ciclo = ciclo0;
			if (NumColision00 < NumColision01) config = 0;
			else config = 1;
		}
		else
		{
			menor = menor1; ciclo = ciclo1;
			if (NumColision10 < NumColision11) config = 2;
			else config = 3;
		}

		printf("menor=%d, ciclo=%d  config=%d \n", menor, ciclo, config);
		getchar();
		printf("Reproduciendo ... \n");
		reproduce(config, ciclo, CapaMenosProfunda);
	}

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GeometryTool()
{
	unsigned int ligak = EdgeInterse[NumEdgeInterse], arista, nivel;
	unsigned int neighb1, neighb2, neighb, target, source;
	printf("La liga %d esta formada por vertices(%d, %d)\n ", ligak, edg[ligak].vertex1, edg[ligak].vertex2);
	printf("Source vertex %d \n ", edg[ligak].vertex1);
	nivel = v[edg[ligak].vertex2].level;
	printf("Hay %d ligas en vertice %d con nivel %d: ", v[edg[ligak].vertex2].liga.size(), edg[ligak].vertex2, nivel);
	source = edg[ligak].vertex1;
	target = edg[ligak].vertex2;
	/*for (int j = 0; j < v[edg[ligak].vertex2].liga.size(); j++)
	{
		arista = v[edg[ligak].vertex2].liga[j];
		printf("<%d(%d),%d(%d)>,",edg[arista].vertex1, v[edg[arista].vertex1].level, edg[arista].vertex2, v[edg[arista].vertex2].level);
		if (edg[arista].vertex1 != edg[ligak].vertex2 && v[edg[arista].vertex1].level==nivel) printf("%d* ", edg[arista].vertex1);
		if (edg[arista].vertex2 != edg[ligak].vertex2 && v[edg[arista].vertex2].level == nivel) printf("%d* ", edg[arista].vertex2);
	}*/
	printf("\n");
	printf("Hay %d nodos en CapaActual+1= %d: \n", NodosPorNivel[CapaActual + 1], CapaActual + 1);
	printf("El primer nodo de este aro es %d, el ultimo nodo es %d\n", PorNivel[CapaActual + 1][0], PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1]-1]);
	if (target == PorNivel[CapaActual + 1][0]) {
		neighb1 = target + 1; 
		neighb2 = PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1];
	}
	else if (target == PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1])
	{
		neighb1 = PorNivel[CapaActual + 1][0];
		neighb2 = target - 1;
	}
	else
	{
		neighb1 = target+1;
		neighb2 = target - 1;
	}
	printf("neighb1=%d, neighb2=%d \n", neighb1, neighb2);

	//calcular la distancia de source a cada uno de los vecinos de target
	Vector3D vect1=Vector3D(v[neighb1].centro[0]- v[source].centro[0], v[neighb1].centro[1] - v[source].centro[1],0);
	Vector3D vect2= Vector3D(v[neighb2].centro[0] - v[source].centro[0], v[neighb2].centro[1] - v[source].centro[1], 0);
	//printf("dist a neighb1=%f, dist a neighb2=%f \n", vect1.length(), vect2.length());
	if (vect1.length() >= vect2.length()) neighb = neighb2;
	else neighb = neighb1;
	//vector from the origin to the chosen neighb 
	//add a unit 
	Vector3D vectNeighb = Vector3D(v[neighb].centro[0], v[neighb].centro[1], v[neighb].centro[2]);
	vectNeighb.unitize();
	despliega_curva = 1;
	LineSegment0.x = v[source].centro[0]; LineSegment0.y = v[source].centro[1];
	LineSegment1.x = v[neighb].centro[0] + vectNeighb.x;
	LineSegment1.y = v[neighb].centro[1] + vectNeighb.y;
	LineSegment2.x = v[target].centro[0]; LineSegment2.y = v[target].centro[1];
	
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Intervalo(unsigned int target, unsigned int Gv1, unsigned int Gv2)
{
	
	unsigned int ligak = EdgeInterse[NumEdgeInterse],start, end, startFirst;
	unsigned int cuenta = 0;
	printf("El primer nodo de este aro es %d, el ultimo nodo es %d\n", PorNivel[CapaActual + 1][0], PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1]);
	bool CreaIntervalo = 0, done=0;
	printf("target=%d, Gv1=%d, Gv2=%d \n", target, Gv1, Gv2);
	for (unsigned int i = PorNivel[CapaActual + 1][0]; i <= PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1]; i++)
	{
		printf("Analizando v[%d]. ", i);
		if (v[i].value)
		{
			printf("Entra a value \n ");
			if (!CreaIntervalo) { 
				start = i; CreaIntervalo = 1; printf("[%d ", start); 
				if (cuenta == 0) {
					startFirst = start; printf("    startFirst=%d ", startFirst); cuenta++;
				}
			}
			else { end = i; CreaIntervalo = 0; 
			printf("checando intervalo [%d, %d] \n", start, end);
			if (target >= start && target <= end)
				if (Gv1 >= start && Gv1 <= end)
					if (Gv2 >= start && Gv2 <= end)
					{
						printf("In [%d, %d]. Llama a rutina <c> \n", start, end); done = 1;  break;
					}
			start = i; CreaIntervalo = 1; printf("[%d ", start);
			}//else
			printf("\n");
		}
		//printf("...cuenta=%d...\n ", cuenta);
		
	}//for i
	if (!done)
	{
		printf("\n Termina ciclo. startFirst=%d, start=%d \n", startFirst, start);
		if (start != startFirst) { printf("Son diferentes, puedo formar un intervalo\n"); 
		
		end = PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1];
		printf("\n checando intervalo [%d, %d]  \n", start, end);
		if (target >= start && target <= end)
			if (Gv1 >= start && Gv1 <= end )
				if (Gv2 >= start && Gv2 <= end)
				{
					printf("In [%d, %d]. Llama a rutina <c> \n", start, end);
				}
				else
				{
					start = PorNivel[CapaActual + 1][0]; end = startFirst;
					printf("\n checando intervalo [%d, %d]  \n", start, end);
					if (target >= start && target <= end)
						if (Gv1 >= start && Gv1 <= end)
							if (Gv2 >= start && Gv2 <= end)
							{
								printf("In [%d, %d]. Llama a rutina <c> \n", start, end);
							}
				}
		}
		else printf("No se puede formar intervalo, ya que solo hay un nodo pivote \n");
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CreateCircumf()
{
	float r1, r2, r12;
	r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda, CApaActual+1
	r2 = (MaxNivel - CapaActual)*5.0 - 4;
	printf("CreateCircumf() < < < < < < < < < < < \n");
	printf("CapaActual=%d con radio %f, CapaActual+1=%d con radio %f \n", 
		CapaActual, r2, CapaActual+1, r1);
	/*for (unsigned int i = 0; i < NumIntersecAro; i++)
	{
		printf(".......Arista....i=%d.....collision=%d\n", i, edg[EdgeInterse[i]].collision);
		printf("%d<%d,%d> choca con arista <%d, %d>\n", EdgeInterse[i], edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInterse[i]].gammav1, edg[EdgeInterse[i]].gammav2);
	}*/
	r12 = r1;
	for (unsigned int i = 0; i < aristasBGcollision; i++)
	{
		edg[CircleEdge[i]].circumf.r = r12 + 5.0 / (aristasBGcollision + 1);
		printf("edge[%d]=<%d,%d> choca con arista <%d, %d>, tendra un r=%f\n", CircleEdge[i], edg[CircleEdge[i]].vertex1, edg[CircleEdge[i]].vertex2, edg[CircleEdge[i]].gammav1, edg[CircleEdge[i]].gammav2, edg[CircleEdge[i]].circumf.r);
		r12 = r12 + 5.0 / aristasBGcollision;

		//vector para interseccion linea-circ
		//Vector3D vect1 = Vector3D(v[edg[CircleEdge[i]].vertex2].centro[0] - v[edg[CircleEdge[i]].vertex1].centro[0], v[edg[CircleEdge[i]].vertex2].centro[1] - v[edg[CircleEdge[i]].vertex1].centro[1], 0);
		Vector3D vect1 = Vector3D(0 - v[edg[CircleEdge[i]].vertex1].centro[0], 0 - v[edg[CircleEdge[i]].vertex1].centro[1], 0);
		vect1.unitize();
		edg[CircleEdge[i]].circumf.Psource.x =  v[edg[CircleEdge[i]].vertex1].centro[0] + (r2-edg[CircleEdge[i]].circumf.r)*vect1.x;
		edg[CircleEdge[i]].circumf.Psource.y = v[edg[CircleEdge[i]].vertex1].centro[1] + (r2 - edg[CircleEdge[i]].circumf.r)*vect1.y;
		//el segmento que va de target a la interseccion de la circumf correspondiente
		vect1 = Vector3D(v[edg[CircleEdge[i]].vertex2].centro[0], v[edg[CircleEdge[i]].vertex2].centro[1], 0);
		vect1.unitize();
		edg[CircleEdge[i]].circumf.Ptarget.x = v[edg[CircleEdge[i]].vertex2].centro[0] + (edg[CircleEdge[i]].circumf.r - r1)*vect1.x;
		edg[CircleEdge[i]].circumf.Ptarget.y = v[edg[CircleEdge[i]].vertex2].centro[1] + (edg[CircleEdge[i]].circumf.r - r1)*vect1.y;
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void EdgeSpheresAdjustment(unsigned int source, unsigned int target, unsigned int Gv1, unsigned int Gv2)
{
	printf("source=%d, target=%d, Gv1=%d, Gv2=%d \n", source, target, Gv1, Gv2);
	Vector3D vect1 = Vector3D(v[target].centro[0]-v[source].centro[0], v[target].centro[1] - v[source].centro[1], v[target].centro[2] - v[source].centro[2]);
	printf("Longitud de arista: %f \n", vect1.length());

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Rota_Arista(unsigned int NumArista, unsigned int source, unsigned int target, unsigned int Gv1, unsigned int Gv2)
{
	Vector3D vect, vectAvanza, Gv;
	printf("source=%d, target=%d, Gv1=%d, Gv2=%d \n", source, target, Gv1, Gv2);
	vect = Vector3D(v[target].centro[0] - v[source].centro[0], v[target].centro[1] - v[source].centro[1], v[target].centro[2] - v[source].centro[2]);
	printf("Longitud de arista: %f \n", vect.length());
	/*vectAvanza.x = v[source].centro[0]; vectAvanza.y = v[source].centro[1];
	edg[NumArista].sphere.push_back(Sphere_type(Vector3D(vectAvanza.x, vectAvanza.y, 0), 0.12));

	vectAvanza.x = (v[Gv1].centro[0] + v[Gv2].centro[0])*0.5; vectAvanza.y = (v[Gv1].centro[1] + v[Gv2].centro[1])*0.5;
	edg[NumArista].sphere.push_back(Sphere_type(Vector3D(vectAvanza.x, vectAvanza.y, 0), 0.12));

	//rotacion arista
	Vector3D plane1, plane2, origen;
	float i1 = 0.0f, r=1;
	plane1 = Vector3D(1, 0, 0); plane2 = Vector3D(0, 1, 0); 

	origen = Vector3D(v[source].centro[0], v[source].centro[1], 0);
	unsigned int mayor, menor, incr, decr;
	printf("El primer nodo de este aro es %d, el ultimo nodo es %d\n", PorNivel[CapaActual + 1][0], PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1]);
	if (Gv1 > Gv2) { mayor = Gv1; menor = Gv2; }
	else { mayor = Gv2; menor = Gv1; }

	if (menor == PorNivel[CapaActual + 1][0]) {
		incr = menor + 1;  decr = mayor 	- 1;
	}
	else if (mayor == PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1]) {
		incr = PorNivel[CapaActual + 1][0] + 1;  decr = menor-1;
	}
	else if (menor == PorNivel[CapaActual + 1][0]) {
		incr = mayor + 1;  decr = PorNivel[CapaActual + 1][NodosPorNivel[CapaActual + 1] - 1] ;
	}
	else { incr = mayor + 1;  decr = menor - 1;}

	printf("incr: %d, decr: %d \n", incr, decr);*/
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void pivot_edge_running(unsigned int NumArista, bool flecha)
{
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	unsigned int collisionOld;
	Vector3D vect; //determina la direccion del vector de movimiento de origin

	vect = edg[NumArista].origin;
	if(flecha) edg[NumArista].originAngle = edg[NumArista].originAngle - 1;
	else edg[NumArista].originAngle = edg[NumArista].originAngle + 1;

	//printf("(pivot running) originAngle=%f, rotations=%d \n", edg[NumArista].originAngle, report[NumEdgeInterse].rotations);
	r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
	i1 = 0.0;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
	i1 = edg[NumArista].originAngle * incr;
	xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
	yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);
	edg[NumArista].origin.x = r1 * xc;
	edg[NumArista].origin.y = r1 * yc;
	edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
	//calcula direccion del vector de movimiento
	vect.x = edg[NumArista].origin.x - vect.x;
	vect.y = edg[NumArista].origin.y - vect.y;
	vect.unitize();
	//deteccion de colision con otras aristas
	collisionOld = edg[NumArista].collisions[edg[NumArista].cont_seguimiento - 1];
	//Detecta_Colision_con_Aristas_Gama(NumArista, vect);
	if (edg[NumArista].collisions[edg[NumArista].cont_seguimiento - 1] < collisionOld && edg[NumArista].cont_seguimiento - 1>1)
	{
		edg[NumArista].num_collisions++;
		printf("Entro al IF. collisionOld=%d, ", collisionOld);
		printf("edg[%d].collisions[%d]=%d, ", NumArista, edg[NumArista].cont_seguimiento - 1, edg[NumArista].collisions[edg[NumArista].cont_seguimiento - 1]);
		printf("edg[%d].num_collisions=%d \n", NumArista, edg[NumArista].num_collisions);
	}
	//ajuste originAngle
	if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
	else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
	if (Calcula_vector_secante(NumArista))
	{
		Crea_mas_vectores_seguimiento(NumArista, r1);
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ext_path_check_pending()
{
	float menor_radio = 100;
	unsigned int arista=0, aristas_ext_path=0;

	printf("- - - void ext_path_check_pending() - - - \n");
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) 
		if (edg[EdgeInterse[i]].new_proc && !edg[EdgeInterse[i]].end) { //aristas pendientes
			                                                                             //pero puede ser una pendiente solitaria
			if (edg[EdgeInterse[i]].radio_path < menor_radio) {
				menor_radio = edg[EdgeInterse[i]].radio_path;
				arista = i; //la arista que comienza tiene el menor radio
			}
			aristas_ext_path++;
	}//for i, if

	if (arista != -1) {
		NumEdgeInterse = arista;
		printf("- - - Continuamos camino pending ext_path pending con EdgeInterse[%d]=%d, hay aristas_ext_path=%d, menor_radio=%f \n", NumEdgeInterse, EdgeInterse[NumEdgeInterse], aristas_ext_path, menor_radio);
		//si radio menor es 0, etonces ...
		ext_path = 1;
		edg[EdgeInterse[NumEdgeInterse]].num_seguidores = aristas_ext_path + 1;
		Calcula_vector_origin_extended_path(EdgeInterse[NumEdgeInterse]);
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ext_path_check(unsigned int NumArista)
{
	float menor_radio = 100;
	unsigned int aristas_ext_path;
	int arista;
	printf("Verifica si hay aristas ext_path \n");
	arista = -1; aristas_ext_path = 0;
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		printf("edg[%d],  radio_path=%f, end=%d. ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path, edg[EdgeInterse[i]].end);
		if (edg[EdgeInterse[i]].collision && EdgeInterse[i] != NumArista) //si es roja y diferente del pivote
			if (edg[EdgeInterse[i]].cont_seguimiento > 0) { //si es follower
				edg[EdgeInterse[i]].ext_path = 1; //marcala como seguidora pendiente
				if (edg[EdgeInterse[i]].radio_path < menor_radio) {
					menor_radio = edg[EdgeInterse[i]].radio_path; 
					arista = i; //la arista que comienza tiene el menor radio
				}
				aristas_ext_path++; //aunque tengas varias, con la del menor radio abarcas a las demas
			}
			else {
				//printf("Requiere nuevo proceso \n"); //puede ser más de una arista, pero solo una considero
			//ext_path = 1;
			}//cont_seguimiento > 0
			printf("\n");
	}//for
	printf("\n");
	if (arista != -1) {
		NumEdgeInterse = arista;
		printf("Continuamos camino con EdgeInterse[%d]=%d, hay aristas_ext_path=%d, menor_radio=%f \n", NumEdgeInterse, EdgeInterse[NumEdgeInterse], aristas_ext_path, menor_radio);
		ext_path = 1;
		edg[EdgeInterse[NumEdgeInterse]].num_seguidores = aristas_ext_path+1;
		Calcula_vector_origin_extended_path(EdgeInterse[NumEdgeInterse]);
	}
	//
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void new_process_check()
{
	float mayor_radio = 0;
	int arista = -1;
	unsigned int Num_aristas_pendientes = 0;

	printf("verificando si hay aristas pendientes \n");
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		printf("edg[%d],  end=%d, ext_path=%d, cont_seguimiento=%d, cont_fin=%d, radio_path=%f \n", EdgeInterse[i],
			edg[EdgeInterse[i]].end, edg[EdgeInterse[i]].ext_path, edg[EdgeInterse[i]].cont_seguimiento,
			edg[EdgeInterse[i]].cont_fin, edg[EdgeInterse[i]].radio_path);
		if (edg[EdgeInterse[i]].collision && !edg[EdgeInterse[i]].end)
		{
			printf("- - pendiente \n");
			edg[EdgeInterse[i]].new_proc = 1;
			Num_aristas_pendientes++;
			arista = i; //Tomare la ultima arista encontrada
		}
		if (edg[EdgeInterse[i]].radio_path > mayor_radio) {
			mayor_radio = edg[EdgeInterse[i]].radio_path;
		}
	} //for
	printf("Hay %d aristas pendientes \n", Num_aristas_pendientes);
	if (Num_aristas_pendientes > 0) {
		new_proc = 1;
		printf("comienzo con arista %d,  edg[%d], con radio mayor %f \n", arista, EdgeInterse[arista], mayor_radio + 1);
		NumEdgeInterse = arista;
		edg[EdgeInterse[NumEdgeInterse]].radio_path = mayor_radio + 0.5;
		arranca = 1;
		Calcula_vector_origin(EdgeInterse[NumEdgeInterse]);
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ext_path_running_pending(unsigned int NumArista, bool flecha) {
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vArista;
	unsigned int collisionOld, aristas_ext_path;
	int arista;
	if ((abs(edg[NumArista].targetAngle - edg[NumArista].originAngle)) > 1.5) {
		if (flecha) edg[NumArista].originAngle = edg[NumArista].originAngle - 1;
		else edg[NumArista].originAngle = edg[NumArista].originAngle + 1;

		printf("(ext_path_running pending) originAngle=%f\n", edg[NumArista].originAngle);
		r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
		i1 = 0.0;
		incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
		i1 = edg[NumArista].originAngle*incr;
		xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
		yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);
		edg[NumArista].origin.x = r1*xc;
		edg[NumArista].origin.y = r1*yc;
		//ajuste originAngle
		if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
		edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
		//. . . . . . . .
		Detecta_Colision_con_Aristas_Gama_Extended_Path_New_Proc(NumArista);
		if (Calcula_vector_secante(NumArista))	Crea_mas_vectores_seguimiento_extended_path_new_proc(NumArista, r1);
		//. . . . . . . .


	}//if target - origin > 1.5
	else {
		printf("FINALIZA ext_path running Pending...............................................................................\n");
		edg[NumArista].end = 1;
		report[NumEdgeInterse].num_collisions = edg[EdgeInterse[NumEdgeInterse]].num_collisions;
		report[NumEdgeInterse].cont_seguimiento = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		report[NumEdgeInterse].radio_path = edg[EdgeInterse[NumEdgeInterse]].radio_path;
		ext_path_new_proc = 0;
	}

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ext_path_running(unsigned int NumArista, bool flecha) {
	float i1, incr, r1;
	double xc, yc;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vArista;
	unsigned int collisionOld, aristas_ext_path;
	int arista;
	
	if ((abs(edg[NumArista].targetAngle - edg[NumArista].originAngle)) > 1.5) {
		if (flecha) edg[NumArista].originAngle = edg[NumArista].originAngle - 1;
		else edg[NumArista].originAngle = edg[NumArista].originAngle + 1;

		printf("(ext_path_running) originAngle=%f\n", edg[NumArista].originAngle);
		r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
		i1 = 0.0;
		incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
		i1 = edg[NumArista].originAngle*incr;
		xc = (double) plane1[0] * 1 * cos(i1*PI) + (double) plane2[0] * 1 * sin(i1*PI);
		yc = (double) plane1[1] * 1 * cos(i1*PI) + (double) plane2[1] * 1 * sin(i1*PI);
		edg[NumArista].origin.x = r1*xc;
		edg[NumArista].origin.y = r1*yc;
		//ajuste originAngle
		if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
		//printf("Target=%f grados, Source=%f grados. Dif=%f  \n", edg[NumArista].targetAngle, edg[NumArista].originAngle, abs(edg[NumArista].targetAngle - edg[NumArista].originAngle));
		edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
		Detecta_Colision_con_Aristas_Gama_Extended_Path(NumArista);
		//ajuste originAngle
		//if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		//else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
		if (Calcula_vector_secante(NumArista))
		{
			Crea_mas_vectores_seguimiento_extended_path(NumArista, r1);
		}
	}//if ext_path running
	else { //ext_path va a finalizar
		edg[NumArista].num_collisions = edg[NumArista].num_collisions + edg[NumArista].collisions[edg[NumArista].cont_seguimiento - 1];
		//informa de las colisiones en el segmento
		printf("edg[%d].collisions[%d]=%d, ", NumArista, edg[NumArista].cont_seguimiento - 1, edg[NumArista].collisions[edg[NumArista].cont_seguimiento - 1]);
		printf("edg[%d].num_collisions=%d \n", NumArista, edg[NumArista].num_collisions);
		//printf("Checando aristas: "); arista = -1; aristas_ext_path = 0;
		//Finalizamos ariste pivote
		printf("FINALIZA ext_path running...............................................................................\n");
		
		edg[NumArista].end = 1;
		//report[NumEdgeInterse].cont_fin = edg[EdgeInterse[NumEdgeInterse]].cont_fin;
		report[NumEdgeInterse].num_collisions = edg[EdgeInterse[NumEdgeInterse]].num_collisions;
		report[NumEdgeInterse].cont_seguimiento = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		report[NumEdgeInterse].radio_path = edg[EdgeInterse[NumEdgeInterse]].radio_path;
		ext_path = 0;
		//
		//checando segmentos y rotaciones
		printf("checando segmentos y rotaciones . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n");
		for (int i = 0; i < EdgeIntersecNum; i++)
			printf("i=%d, segmentos=%d, rot=%d \n",i, report[i].cont_seguimiento, report[i].rotations);
		//
		//new_process_check(NumArista); //arista pendiente ?
	}//else ext_path NOT  running
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void new_proc_running(unsigned int NumArista, bool flecha) {
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	//Vector3D vArista;
	Vector3D vect; //determina la direccion del vector de movimiento de origin
	unsigned int collisionOld;
	//int arista;
	if ((abs(edg[NumArista].targetAngle - edg[NumArista].originAngle)) > 1.5) {
		if (flecha) edg[NumArista].originAngle = edg[NumArista].originAngle - 1;
		else edg[NumArista].originAngle = edg[NumArista].originAngle + 1;

		printf("(new_proc_running) originAngle=%f\n", edg[NumArista].originAngle);
		r1 = (MaxNivel - CapaActual - 1)*5.0 - 4; //radio menor, capa mas profunda
		i1 = 0.0;
		incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
		i1 = edg[NumArista].originAngle*incr;
		xc = plane1[0] * 1 * cos(i1*PI) + plane2[0] * 1 * sin(i1*PI);
		yc = plane1[1] * 1 * cos(i1*PI) + plane2[1] * 1 * sin(i1*PI);
		edg[NumArista].origin.x = r1*xc;
		edg[NumArista].origin.y = r1*yc;
		//edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
		//ajuste originAngle
		if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
		printf("Target=%f grados, Source=%f grados. Dif=%f  \n", edg[NumArista].targetAngle, edg[NumArista].originAngle, abs(edg[NumArista].targetAngle - edg[NumArista].originAngle));
		edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
		Detecta_Colision_con_Aristas_Gama_New_Proc(NumArista);
		//ajuste originAngle
		if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;
		if (Calcula_vector_secante(NumArista)) Crea_mas_vectores_seguimiento_new_proc(NumArista, r1);
	}//if new_proc running
	else {//new_proc va a finalizar
		edg[EdgeInterse[NumEdgeInterse]].num_collisions = edg[EdgeInterse[NumEdgeInterse]].num_collisions + edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1];
		//informa de las colisiones en el segmento
		printf("edg[%d].collisions[%d]=%d, ", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1, edg[EdgeInterse[NumEdgeInterse]].collisions[edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento - 1]);
		printf("edg[%d].num_collisions=%d \n", EdgeInterse[NumEdgeInterse], edg[EdgeInterse[NumEdgeInterse]].num_collisions);
		printf("FINALIZA new_proc .......................................................... \n: ");
		report[NumEdgeInterse].num_collisions = edg[EdgeInterse[NumEdgeInterse]].num_collisions;
		report[NumEdgeInterse].cont_seguimiento = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		report[NumEdgeInterse].radio_path = edg[EdgeInterse[NumEdgeInterse]].radio_path;

		edg[EdgeInterse[NumEdgeInterse]].end = 1; 		//Finalizamos arista pivote
		edg[EdgeInterse[NumEdgeInterse]].cont_fin = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento = 0;
		new_proc = 0;
	}//else ext_path NOT  running
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ReporteGraph(const char *p_filename)
{

	unsigned int count_GammaBlack = 0, count_GammaRed = 0;
	printf("Creando archivo txt %s ...\n", p_filename);
	using namespace std;
	//std::stringstream ss;
	//ss << std::tab;
	ofstream fout(p_filename);
	//printf("\n - - NumNodes=%d, NumEdges=%d, MaxNivel=%d\n", NumNodes, NumEdges, MaxNivel);
	fout << "NumNodes=" << NumNodes << "   NumEdges=" << NumEdges << "   MaxNivel=" << MaxNivel << endl;
	//fout << "SG" << "\t" << "Capas" << "\t" << "Nodos" << "\t" << "Aristas" << "\t" << "Costo" << "\t" << "Tiempo" << "\t\t" << "CostoS" << endl;
	fout << "Layers: " << CapaActual << ","  << CapaActual+1 <<"     Aristas" << "     Collisions" << endl;
	fout << "Edges" << "\t\t" << "Segments" << "\t" << "GammaBlack" << "\t" << "GamaRed" << "\t" << "Type" << "\t" << "Iterations" << "\t" << "radio_path" << endl;
	for (int i = 0; i < EdgeIntersecNum; i++)
	{
		//printf("edg[%d],  radio_path=%f, end=%d. ", EdgeInterse[i], edg[EdgeInterse[i]].radio_path, edg[EdgeInterse[i]].end);		
		fout << "edg[" << EdgeInterse[i] << "]";
		if (edg[EdgeInterse[i]].collision) fout << "  (R) ";
		else fout << "      ";
		fout << "\t" << report[i].cont_seguimiento <<  "\t" << report[i].num_collisions << "\t\t" << report[i].num_collisions_red;
		fout << "\t" << report[i].tipo  << "\t" << report[i].rotations << "\t" << report[i].radio_path;
		//if (i == NumEdgeInterse) fout << count << "  -> ";
		//if (report[i].tipo == 1) fout << "  -> ";
		
		fout << endl;
		count_GammaBlack += report[i].num_collisions;
		count_GammaRed += report[i].num_collisions_red;
	}
	fout << endl;
	fout << "Crossings \t\t\t" << count_GammaBlack / 2 << "\t\t" << count_GammaRed/2 << "\t =" << count_GammaBlack/2 + count_GammaRed/2;
	fout.close();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ReporteGraph_Start(const char *p_filename)
{
	unsigned count = 0;
	float r1;
	//colision con anullos internos
	count = CapaActual + 2; r1 = (MaxNivel - CapaActual - 1)*5.0 - 4;
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
		 //colision con anullos internos END

	printf("Creando archivo txt %s ...\n", p_filename);
	ofstream fout(p_filename);
	count = 0;
	fout << "Layers: " << CapaActual << "," << CapaActual + 1 << "     Aristas" << "     Collisions" << endl;
	fout << "Edges" << "\t\t" << "Segments"<< "\t" <<  "Beta" << "\t" << "Gamma" << "\t" << "inner" << "\t" << "Type" << "\t" << "Iterations" << "\t" << "radio_path" << endl;
	for (int i = 0; i < EdgeIntersecNum; i++)
	{
		fout << "edg[" << EdgeInterse[i] << "]";
		if (edg[EdgeInterse[i]].collision) fout << "  (R) ";
		else fout << "      ";
		fout << "\t" << report[i].cont_seguimiento << "\t" << report[i].num_collisions << "\t" << report[i].num_collisions_red;
		fout << "\t" << report[i].inner_circles;
		fout << "\t" << report[i].tipo << "\t" << report[i].rotations << "\t" << report[i].radio_path;
		
		fout << endl;
		count += report[i].inner_circles;
	}
	fout << endl;
	fout << "Crossings \t\t\t" << NumIntersecAro << "\t" << NumIntersec << "\t\t =" << NumIntersecAro+NumIntersec << endl;
	fout << "\t\t\t\t\t"  << count << "\t=" << NumIntersecAro + NumIntersec + count << endl;
	fout << endl;
	fout.close();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::