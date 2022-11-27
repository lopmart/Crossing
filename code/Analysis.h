#pragma once
//::::::::::::::::::::::::::::::::::::::::::::::::::::::

void checkCollision(float a, float b, float c, float x, float y, float radius)
{
	// Finding the distance of line from center. 
	float dist = (abs(a * x + b * y + c)) / sqrt(a * a + b * b);
	printf("dist=%2.2f ", dist);
	// Checking if the distance is less than,  
	// greater than or equal to radius. 
	if (radius == dist)
		cout << "Touch" << endl;
	else if (radius > dist)
		cout << "Intersect" << endl;
	else
		cout << "Outside" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Circle_Edge_Collision(unsigned int NumArista, float radio)
{
	printf("Edge %d. ", NumArista);
	float Ax, Ay, Bx, By, Cx, Cy, u;
	Ax = v[edg[NumArista].vertex1].centro[0]; Ay = v[edg[NumArista].vertex1].centro[1];
	Bx = v[edg[NumArista].vertex2].centro[0]; By = v[edg[NumArista].vertex2].centro[1];
	Cx = 0; Cy = 0;
	u = ((Cx - Ax) * (Bx - Ax) + (Cy - Ay) * (By - Ay)) / ((Bx - Ax) * (Bx - Ax) + (By - Ay) * (By - Ay));
	printf("u=%2.2f, ", u);
	printf("\n");
	//si u esta entre 0 y 1, entonces atraviesa el circulo
	//if (u > 0 && u < 1) {
	float dx, dy, dr, D, xval, yval, Discriminante;

	float m, result;
	float a, b, c;
	m = (v[edg[NumArista].vertex2].centro[1] - v[edg[NumArista].vertex1].centro[1]) / (v[edg[NumArista].vertex2].centro[0] - v[edg[NumArista].vertex1].centro[0]);
	result = atan(m) * 180 / PI;
	printf("m=%2.2f, arc tg = %f, ", m, result);
	//valores de la ecuacion
	a = m; b = -1;
	c = v[edg[NumArista].vertex2].centro[1] - m * (v[edg[NumArista].vertex2].centro[0]);
	xval = v[edg[NumArista].vertex2].centro[0];
	yval = v[edg[NumArista].vertex2].centro[1];

	float x0 = -a * c / (a * a + b * b);
	float y0 = -b * c / (a * a + b * b);
	float EPS = 0.05;
	float r = radio;

	if (c * c > r * r * (a * a + b * b) + EPS)
		puts("no points");
	else if (abs(c * c - r * r * (a * a + b * b)) < EPS) {
		puts("1 point");
		cout << x0 << ' ' << y0 << '\n';
		pointCircle.push_back(Vector3D(x0, y0, 0));
	}
	else {
		double d = r * r - c * c / (a * a + b * b);
		double mult = sqrt(d / (a * a + b * b));
		double ax, ay, bx, by;
		ax = x0 + b * mult;
		bx = x0 - b * mult;
		ay = y0 - a * mult;
		by = y0 + a * mult;
		puts("2 points");
		cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
		pointCircle.push_back(Vector3D(ax, ay, 0));
		pointCircle.push_back(Vector3D(bx, by, 0));
	}
	//}//u esta entre cero y uno
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Calcula_vector_counterclockwise(unsigned int NumArista, float radio)
{
	Vector3D vect;
	float angulo, originAngle;
	bool EdgRed = 0;

	EdgRed = edg[NumArista].collision;

	printf("Edge %d. ", NumArista);
	//---------------------------------------------vertex1
	vect.z = 0;
	vect.x = v[edg[NumArista].vertex1].centro[0];
	vect.y = v[edg[NumArista].vertex1].centro[1];
	vect.unitize();
	angulo = vect.dotproduct(Vector3D(1, 0, 0));
	originAngle = acos(angulo) * 180 / 3.1415;

	//ajuste para cuadrantes II, IV
	if (vect.x >= 0) {
		if (vect.y >= 0) {}
		else { originAngle = 360 - originAngle; }
	}
	else {
		if (vect.y >= 0) {}
		else { originAngle = 360 - originAngle; }
	}

	//printf("angulo=%f, originAngle= %f grados \n", angulo, originAngle);
	printf("vertex1= %f grados, ", originAngle);
	Angles.push_back(originAngle);
	//edg_angles.push_back(Analiza_Par(originAngle, NumArista, 0));
	edge_angles.push_back(subject(originAngle, NumArista, 0, edg[NumArista].vertex1, EdgRed));
	//---------------------------------------------vertex2
	vect.z = 0;
	vect.x = v[edg[NumArista].vertex2].centro[0];
	vect.y = v[edg[NumArista].vertex2].centro[1];
	vect.unitize();
	angulo = vect.dotproduct(Vector3D(1, 0, 0));
	originAngle = acos(angulo) * 180 / 3.1415;

	//ajuste para cuadrantes II, IV
	if (vect.x >= 0) {
		if (vect.y >= 0) {}
		else { originAngle = 360 - originAngle; }
	}
	else {
		if (vect.y >= 0) {}
		else { originAngle = 360 - originAngle; }
	}

	printf("vertex2= %f grados \n", originAngle);
	Angles.push_back(originAngle);
	//---------------------------------------------ordena angulos
	//edg_angles.push_back(Analiza_Par(originAngle, NumArista, 1));
	edge_angles.push_back(subject(originAngle, NumArista, 1, edg[NumArista].vertex2, EdgRed));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void insert_edg_in_match(unsigned int NumArista)
{
	bool ya_existe = 0;
	for (unsigned int i = 0; i < Num_match; i++)
	{
		if (match[i].NumEdg == NumArista) { ya_existe = 1; break; }
	}
	if (!ya_existe) {
		match[Num_match].NumEdg = NumArista;
		printf("match[%d].NumEdg=%d \n", Num_match, match[Num_match].NumEdg);
		Num_match++;
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
unsigned int index_edg_in_match(unsigned int NumArista)
{
	for (unsigned int i = 0; i < Num_match; i++)
		if (match[i].NumEdg == NumArista) return i;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void List_Clockwise(unsigned int NumArista, unsigned int pos_inicial)
{
	bool stop = 0;
	unsigned int pos = pos_inicial, len = edge_angles.size(), NumCol = 0, NumFollow = 0, temp = 0;
	struct {
		unsigned int NumEdg;// match = 0;
		bool Label;
	} follower[10];
	unsigned int index = 0;
	for (unsigned int i = 0; i < 10; i++) {
		match[i].NumEdg = 0;
		match[i].count = 0;
	}
	printf("------------ Clockwise run\n");
	printf("Comenzamos en posicion  %d, edg[%d] \n", pos, edge_angles[pos].edge);
	while (!stop) {
		if (pos > 0) pos--;
		else pos = len - 1;
		//printf("posicion  %d, edg[%d], red=%d, end=%d ", pos, edge_angles[pos].edge, edge_angles[pos].red, edge_angles[pos].label);
		printf("posicion  %d, edg[%d], ", pos, edge_angles[pos].edge);
		if (edge_angles[pos].red) printf("red, ");
		else printf("black, ");

		if (edge_angles[pos].label) printf("end ");
		else printf("start ");

		if (edge_angles[pos].edge == NumArista && edge_angles[pos].red && edge_angles[pos].label) stop = 1;
		if (!edge_angles[pos].red && edge_angles[pos].label) {
			printf("  X  ");
			NumCol++;
		}
		if (edge_angles[pos].edge != NumArista && edge_angles[pos].red) {
			if (!edge_angles[pos].label) {
				printf("  Follower  ");
				NumFollow++;
				follower[temp].NumEdg = edge_angles[pos].edge;
				follower[temp].Label = 0;
			}
			else {
				follower[temp].NumEdg = edge_angles[pos].edge;
				follower[temp].Label = 1;
			}
			insert_edg_in_match(follower[temp].NumEdg);
			temp++;
		}
		printf("\n");
	}//while
	//tenemos a las aristas
	//for (unsigned int i = 0; i < EdgeIntersecNum; i++)
		//printf("edg[%d], ", EdgeInterse[i]);

	printf("\n");
	//printf("NumCol=%d, NumFollow=%d \n", NumCol, NumFollow);
	for (unsigned int i = 0; i < temp; i++) {
		printf("%d. edg[%d] ", i, follower[i].NumEdg);
		if (follower[i].Label) {
			index = index_edg_in_match(follower[i].NumEdg);
			printf(" end, index match[%d] ", index);
			if (match[index].count == 1) match[index].count++;
		}
		else {
			index = index_edg_in_match(follower[i].NumEdg);
			printf(" start, index match[%d] ", index);
			match[index].count++;
		}
		printf("  match=%d \n", match[index].count);
	}//for
	printf("\n");
	printf("NumCol=%d, NumFollow=%d \n", NumCol, NumFollow);
	for (unsigned int i = 0; i < Num_match; i++) {
		printf(" match[%d]=%d, %d ", i, match[i].NumEdg, match[i].count);
		if (match[i].count == 2) printf("follower completada \n");
		else if (match[i].count == 1) printf("follower pendiente \n");
		else printf("arista pendiente \n");
	}

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void List_Counter(unsigned int NumArista, unsigned int pos_inicial)
{
	bool stop = 0;
	unsigned int pos = pos_inicial, len = edge_angles.size(), NumCol = 0, NumFollow = 0, temp = 0;
	struct {
		unsigned int NumEdg;// match = 0;
		bool Label;
	} follower[10];
	unsigned int index = 0;
	for (unsigned int i = 0; i < 10; i++) {
		match[i].NumEdg = 0;
		match[i].count = 0;
	}
	printf("------------ Counter Clockwise run\n");
	printf("Comenzamos en posicion  %d, edg[%d] \n", pos, edge_angles[pos].edge);
	while (!stop) {
		if (pos < len - 1) pos++;
		else pos = 0;
		//printf("posicion  %d, edg[%d], red=%d, end=%d ", pos, edge_angles[pos].edge, edge_angles[pos].red, edge_angles[pos].label);
		printf("posicion  %d, edg[%d], ", pos, edge_angles[pos].edge);
		if (edge_angles[pos].red) printf("red, ");
		else printf("black, ");

		if (edge_angles[pos].label) printf("end ");
		else printf("start ");

		if (edge_angles[pos].edge == NumArista && edge_angles[pos].red && edge_angles[pos].label) stop = 1;
		if (!edge_angles[pos].red && edge_angles[pos].label) {
			printf("  X  ");
			NumCol++;
		}
		if (edge_angles[pos].edge != NumArista && edge_angles[pos].red) {
			if (!edge_angles[pos].label) {
				printf("  Follower  ");
				NumFollow++;
				follower[temp].NumEdg = edge_angles[pos].edge;
				follower[temp].Label = 0;
			}
			else {
				follower[temp].NumEdg = edge_angles[pos].edge;
				follower[temp].Label = 1;
			}
			insert_edg_in_match(follower[temp].NumEdg);
			temp++;
		}
		printf("\n");
	}//while
	 //tenemos a las aristas
	 //for (unsigned int i = 0; i < EdgeIntersecNum; i++)
	 //printf("edg[%d], ", EdgeInterse[i]);

	printf("\n");
	//printf("NumCol=%d, NumFollow=%d \n", NumCol, NumFollow);
	for (unsigned int i = 0; i < temp; i++) {
		printf("%d. edg[%d] ", i, follower[i].NumEdg);
		if (follower[i].Label) {
			index = index_edg_in_match(follower[i].NumEdg);
			printf(" end, index match[%d] ", index);
			if (match[index].count == 1) match[index].count++;
		}
		else {
			index = index_edg_in_match(follower[i].NumEdg);
			printf(" start, index match[%d] ", index);
			match[index].count++;
		}
		printf("  match=%d \n", match[index].count);
	}//for
	printf("\n");
	printf("NumCol=%d, NumFollow=%d \n", NumCol, NumFollow);
	for (unsigned int i = 0; i < Num_match; i++) {
		printf(" match[%d]=%d, %d ", i, match[i].NumEdg, match[i].count);
		if (match[i].count == 2) printf("follower completada \n");
		else if (match[i].count == 1) printf("follower pendiente \n");
		else printf("arista pendiente \n");
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Algoritmo_Aristas()
{
	unsigned int red_inicio[10], i = 0, j = 0;

	printf("------------ Algoritmo_Aristas con %d elementos\n", edge_angles.size());
	printf("red angle  edge  start/end vertex \n");
	for (auto x : edge_angles) {
		printf("%d  %3.1f   %d    %d      %d ", x.red, x.angle, x.edge, x.label, x.end);
		//if (x.red && !x.label) printf("->");
		if (x.red) {
			if (!x.label) { printf("->");  red_inicio[i] = j; i++; }
			else printf("<-");
		}

		printf("\n"); j++;
	}//for
	printf("\n %d aristas rojas: \n", i);
	for (unsigned int j = 0; j < i; j++) {
		printf("edg[%d]. posicion  %d \n", edge_angles[red_inicio[j]].edge, red_inicio[j]);
	}

	bool si_esta = 0;
	unsigned int arista_analizada;
	//recorremos la lista desde alguna de las aristas
	if (i > 0) { //hay al menos una arista roja
		arista_analizada = 0;
		List_Counter(edge_angles[red_inicio[arista_analizada]].edge, red_inicio[arista_analizada]);
		if (Num_match < i - 1) {
			printf("arista pendiente \n"); si_esta = 0;
			for (unsigned int j = 0; j < i; j++) {
				if (j != arista_analizada) {
					for (unsigned int k = 0; k < Num_match; k++)
						if (edge_angles[red_inicio[j]].edge == match[k].NumEdg) si_esta = 1;
					if (!si_esta) printf("(%d) \n", edge_angles[red_inicio[j]].edge);
				}//j!=1
				si_esta = 0;
			}//for j
		}//if
		Num_match = 0;
		List_Clockwise(edge_angles[red_inicio[arista_analizada]].edge, red_inicio[arista_analizada]);
		if (Num_match < i - 1) {
			printf("arista pendiente \n"); si_esta = 0;
			for (unsigned int j = 0; j < i; j++) {
				if (j != arista_analizada) {
					for (unsigned int k = 0; k < Num_match; k++)
						if (edge_angles[red_inicio[j]].edge == match[k].NumEdg) si_esta = 1;
					if (!si_esta) printf("(%d) \n", edge_angles[red_inicio[j]].edge);
				}//j!=1
				si_esta = 0;
			}//for j
		}//if
		Num_match = 0;
	}//i>0
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Global_Process_Pending(unsigned int NumArista, bool flecha)
{
	float mayor_radio = 0;
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	//unsigned int Num_aristas_pendientes = 0;

	if ((abs(edg[NumArista].targetAngle - edg[NumArista].originAngle)) > 1.5) {
		if (flecha) edg[NumArista].originAngle = edg[NumArista].originAngle - 1;
		else edg[NumArista].originAngle = edg[NumArista].originAngle + 1;

		//printf("(new_proc_running) originAngle=%f\n", edg[NumArista].originAngle);
		r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
		i1 = 0.0;
		incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI
		i1 = edg[NumArista].originAngle * incr;
		xc = plane1[0] * 1 * cos(i1 * PI) + plane2[0] * 1 * sin(i1 * PI);
		yc = plane1[1] * 1 * cos(i1 * PI) + plane2[1] * 1 * sin(i1 * PI);
		edg[NumArista].origin.x = r1 * xc;
		edg[NumArista].origin.y = r1 * yc;
		printf("\n .  . . . . . Global_Process_Pending( ),  NumArista=%d, edg[NumArista].cont_seguimiento - 1=%d \n", NumArista, edg[NumArista].cont_seguimiento - 1);
		edg[NumArista].seguimiento[edg[NumArista].cont_seguimiento - 1] = edg[NumArista].origin;
		printf("edg[%d].seguimiento[%d] = (%2.2f, %2.2f) \n", NumArista, edg[NumArista].cont_seguimiento - 1, edg[NumArista].origin.x, edg[NumArista].origin.y);
		Detecta_Colision_con_Aristas_Gama_New_Proc(NumArista);

		if (edg[NumArista].originAngle > 360) edg[NumArista].originAngle = 0;
		else if (edg[NumArista].originAngle < 0) edg[NumArista].originAngle = 360;

		if (Calcula_vector_secante(NumArista)) Crea_mas_vectores_seguimiento_new_proc(NumArista, r1);
	}//if new_proc running
	else {//new_proc va a finalizar
		printf("FINALIZA new_proc .......................................................... \n: ");
		//
		report[NumEdgeInterse].num_collisions = edg[EdgeInterse[NumEdgeInterse]].num_collisions;
		report[NumEdgeInterse].cont_seguimiento = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		report[NumEdgeInterse].radio_path = edg[EdgeInterse[NumEdgeInterse]].radio_path;

		//
		edg[EdgeInterse[NumEdgeInterse]].end = 1; 		//Finalizamos arista pivote
		edg[EdgeInterse[NumEdgeInterse]].cont_fin = edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento;
		edg[EdgeInterse[NumEdgeInterse]].cont_seguimiento = 0;
		Num_aristas_pendientes--; //decrementa pending edges
		printf("\n < < < < < < Global_Process_Pending(%d),  Num_aristas_pendientes=%d \n", NumArista, Num_aristas_pendientes);

		new_proc = 0;
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checa_colisiones_segmentos_pivote_vs_gamma(int NumAr)
{
	Vector3D vect, vect1, vect2;
	float r1, nodo1X, nodo1Y, nodo2X, nodo2Y;
	double dot;

	r1 = (MaxNivel - CapaActual - 1) * 5.0 - 1; //radio intermedio
	edg[NumAr].num_collisions = 0;
	for (unsigned int i = 1; i < edg[NumAr].cont_fin; i++) //comienza de i=1 para formar segm con i=0
	{
		edg[NumAr].collisions[i - 1] = 0; //inicializa contador de colisiones por segmento
		//segmento i-1, i
		vect.x = edg[NumAr].stop[i - 1].x;
		vect.y = edg[NumAr].stop[i - 1].y;
		vect.z = edg[NumAr].stop[i - 1].z;

		vect1.x = edg[NumAr].stop[i].x;
		vect1.y = edg[NumAr].stop[i].y;
		vect1.z = edg[NumAr].stop[i].z;

		//vect.unitize(); vect1.unitize();
		//nodo1X = r1 * vect.x; nodo1Y = r1 * vect.y;
		//nodo2X = r1 * vect1.x; nodo2Y = r1 * vect1.y;

		cout << "\n segmento i=" << i - 1 << "-" << i << ", edg[" << NumAr << "].cont_fin=" << edg[NumAr].cont_fin << "  x=" << vect.x << " y=" << vect.y << " z=" << vect.z << "\n";
		//cout << "     j=";
		for (unsigned int j = 0; j < EdgeIntersecNum; j++)
		{
			//cout << j  << " edg[" << EdgeInterse[j] << "] ";
			if (EdgeInterse[j] != NumAr)
				if (Intersect_one_vertex(edg[EdgeInterse[j]].vertex1, edg[EdgeInterse[j]].vertex2,
					vect.x, vect.y, vect1.x, vect1.y)) {
					//cout << " (collision) "; 
					edg[NumAr].collisions[i - 1]++;
					edg[NumAr].num_collisions++;
				}
			//cout << ",";
		}//for j
		cout << "edg[" << NumAr << "].collisiones[" << i - 1 << "]=" << edg[NumAr].collisions[i - 1] << "\n";
	}//for
	//el ultimo segmento
	vect.x = edg[NumAr].stop[edg[NumAr].cont_fin - 1].x;
	vect.y = edg[NumAr].stop[edg[NumAr].cont_fin - 1].y;
	vect.z = edg[NumAr].stop[edg[NumAr].cont_fin - 1].z;

	vect1.x = edg[NumAr].seguimiento[edg[NumAr].cont_fin - 1].x;
	vect1.y = edg[NumAr].seguimiento[edg[NumAr].cont_fin - 1].y;
	vect1.z = edg[NumAr].seguimiento[edg[NumAr].cont_fin - 1].z;
	cout << " Ultimo segmento \n ";
	edg[NumAr].collisions[edg[NumAr].cont_fin - 1] = 0; //inicializa contador de colisiones por segmento
	for (unsigned int j = 0; j < EdgeIntersecNum; j++)
	{
		//cout << j << " edg[" << EdgeInterse[j] << "] ";
		if (EdgeInterse[j] != NumAr)
			if (Intersect_one_vertex(edg[EdgeInterse[j]].vertex1, edg[EdgeInterse[j]].vertex2,
				vect.x, vect.y, vect1.x, vect1.y)) {
				//cout << " (collision) "; 
				edg[NumAr].collisions[edg[NumAr].cont_fin - 1]++;
				edg[NumAr].num_collisions++;
			}
		cout << ",";
	}//for j
	cout << "edg[" << NumAr << "].collisiones[" << edg[NumAr].cont_fin - 1 << "]=" << edg[NumAr].collisions[edg[NumAr].cont_fin - 1] << "\n";
	cout << "edg[" << NumAr << "].num_collisions =" << edg[NumAr].num_collisions << "\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Global_Process_todo_el_circulo() {
	int NumAr = EdgeInterse[NumEdgeInterse];
	for (int i = 0; i < edg[0].cont_seguimiento; i++) {
		edg[NumAr].stop[i].x = edg[0].seguimiento[i].x;
		edg[NumAr].stop[i].y = edg[0].seguimiento[i].y;
	}
	edg[NumAr].cont_fin = edg[0].cont_seguimiento;
	printf("edg[NumAr].cont_fin=%d\n", edg[NumAr].cont_fin);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Global_Process()
{
	char tecla; int NumAr = EdgeInterse[NumEdgeInterse];
	float ang = 360.0 / edg[0].cont_seguimiento;
	int PuntoSegmentoInicio, PuntoSegmentoFin;
	bool flecha = sentido_rot;

	//printf("ang=%f. flecha=%d \n", ang, flecha);
	//printf("edg[%d].originAngle=%f. edg[%d].targetAngle=%f. \n", NumAr, edg[NumAr].originAngle, NumAr, edg[NumAr].targetAngle);

	if (flecha) PuntoSegmentoInicio = floor(edg[NumAr].originAngle / ang);
	else PuntoSegmentoInicio = ceil(edg[NumAr].originAngle / ang); //anti horario
	//ceil((edg[NumAr].originAngle / ang)+1.0); 
	if (flecha) PuntoSegmentoFin = ceil(edg[NumAr].targetAngle / ang);
	else PuntoSegmentoFin = floor(edg[NumAr].targetAngle / ang); //anti horario
	//si te pasas del numero de segmentos
	if (PuntoSegmentoInicio >= edg[0].cont_seguimiento) PuntoSegmentoInicio = 0;
	if (PuntoSegmentoFin >= edg[0].cont_seguimiento) PuntoSegmentoFin = 0;

	//printf("PuntoSegmentoInicio= %d, PuntoSegmentoFin=%d\n", PuntoSegmentoInicio, PuntoSegmentoFin);


	edg[NumAr].stop[0].x = v[edg[NumAr].vertex1].centro[0]; //v[edg[NumAr].gammav1].centro[0];
	edg[NumAr].stop[0].y = v[edg[NumAr].vertex1].centro[1]; // v[edg[NumAr].gammav1].centro[1];
	//printf("stop[%d] = seguimiento[vertex1]\n", 0);

	edg[NumAr].stop[1].x = edg[0].seguimiento[PuntoSegmentoInicio].x;
	edg[NumAr].stop[1].y = edg[0].seguimiento[PuntoSegmentoInicio].y;
	//printf("stop[%d] = seguimiento[%d]\n", 1, PuntoSegmentoInicio);

	int j = 2;
	//llena el recorrido
	if (flecha) { //horario
		//printf("Sentido Horario \n");
		if (PuntoSegmentoInicio > PuntoSegmentoFin) {
			//printf("Sentido Horario: PuntoSegmentoInicio > PuntoSegmentoFin \n");
			for (int i = PuntoSegmentoInicio - 1; i >= PuntoSegmentoFin; i--) {
				edg[NumAr].stop[j].x = 1.0 * edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}
		}
		else {
			//printf("Sentido Horario: PuntoSegmentoFin > PuntoSegmentoInicio \n");
			for (int i = PuntoSegmentoInicio - 1; i >= 0; i--) {
				edg[NumAr].stop[j].x = edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}//llegas al primer punto de los segmentos
			//y sigues decrementando hasta llegar al nodo destino
			for (int i = edg[0].cont_seguimiento - 1; i >= PuntoSegmentoFin; i--) {
				edg[NumAr].stop[j].x = edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}
		}
	}//sentido horario
	else { //anti horario
		if (PuntoSegmentoFin > PuntoSegmentoInicio) {
			for (int i = PuntoSegmentoInicio + 1; i < PuntoSegmentoFin + 1; i++) {
				edg[NumAr].stop[j].x = edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}
		}
		else {
			//printf("j= ");
			for (int i = PuntoSegmentoInicio + 1; i < edg[0].cont_seguimiento; i++) {
				edg[NumAr].stop[j].x = edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}
			//llega al final, ahora hasta alcanza el Punto final del destino
			for (int i = 0; i < PuntoSegmentoFin + 1; i++) {
				edg[NumAr].stop[j].x = edg[0].seguimiento[i].x;
				edg[NumAr].stop[j].y = edg[0].seguimiento[i].y;
				//printf("stop[%d] = seguimiento[%d]\n", j, i);
				j++;
			}

		}//else

	}//else anti horario
	//finalmente unimos con el nodo destino
	edg[NumAr].stop[j].x = v[edg[NumAr].vertex2].centro[0];
	edg[NumAr].stop[j].y = v[edg[NumAr].vertex2].centro[1];
	//printf("stop[%d] = seguimiento[vertex2]\n", j);
	j++;
	edg[NumAr].cont_fin = j;
	//if (imprime) { printf("edg[NumAr].cont_fin=%d\n", edg[NumAr].cont_fin); imprime = 0; }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Start_Running_Pending_Edge(unsigned int NumArista, bool flecha)
{
	float mayor_radio = 0;
	float xc, yc, i1, incr, r1;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	unsigned int inter, count;

	//encontrar el mayor radio y marcar las aristas pendientes
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		if (edg[EdgeInterse[i]].collision && !edg[EdgeInterse[i]].end)
		{
			edg[EdgeInterse[i]].new_proc = 1;
			Num_aristas_pendientes++;

		}
		if (edg[EdgeInterse[i]].radio_path > mayor_radio)
			mayor_radio = edg[EdgeInterse[i]].radio_path;
	}//for
	 //Elegí una arista pendiente con <x>
	printf("\n > > > > > > Start_Running_Pending_Edge(%d),  con radio mayor %f, Num_aristas_pendientes=%d \n", NumArista, mayor_radio + 1, Num_aristas_pendientes);
	edg[EdgeInterse[NumEdgeInterse]].radio_path = mayor_radio + 0.5;
	arranca = 1;
	printf("NumEdgeInterse=%d, NumArista=%d\n", NumEdgeInterse, NumArista);
	Calcula_vector_origin(NumArista);
	//comienza rotacion arista pendiente------------------------------------------------------------------
	inter = sentido_rot; //sentido de la rotacion 1 clock, 0 counter clock
	while (new_proc) {
		Global_Process_Pending(NumArista, inter);
		report[NumEdgeInterse].rotations++;
		//_getch();
	}//while
	report[NumEdgeInterse].tipo = 4;
	printf("report[%d].rotations: %d \n", NumEdgeInterse, report[NumEdgeInterse].rotations);
	unsigned int rotPivot = report[NumEdgeInterse].rotations;
	ext_path_check_pending(); //verifica si hay aristas ext_path pending
	//ext path BEGIN---------------------------------------------------------------------------------------------------
	count = 0;
	if (Num_aristas_pendientes > 0) {
		new_proc = 1; ext_path_new_proc = 1;
		//Hay que elegir la arista con menor radio_path
		printf("Comienza el otro while - - - - - - - - - - - - - - - -\n");
		if (edg[EdgeInterse[NumEdgeInterse]].targetAngle < edg[EdgeInterse[NumEdgeInterse]].originAngle) inter = 1;
		else inter = 0;
		while (ext_path_new_proc) {
			ext_path_running_pending(EdgeInterse[NumEdgeInterse], inter);
			report[NumEdgeInterse].rotations++;
		}
		report[NumEdgeInterse].tipo = 5;
		printf("ext_path report[%d].rotations: %d \n", NumEdgeInterse, report[NumEdgeInterse].rotations);
		//agregas las rotaciones de su pivote
		report[NumEdgeInterse].rotations += rotPivot;
	} //aun seguire en pending edges
	//ext path END---------------------------------------------------------------------------------------------------
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Start_Running() {
	unsigned int index = NumEdgeInterse, inter, count = 0;
	pivot_running = 0; ext_path = 0;
	//printf("\n Comenzamos proceso con arista edg[%d], EdgeInterse[%d], index=%d - - - - - - - - - - - - - -\n ", EdgeInterse[NumEdgeInterse], NumEdgeInterse, index);
	pivot_running = 1; //comienza proceso ------------------------------------------------------------------
	Global_Process();
	count++;
	report[NumEdgeInterse].rotations = count;
	unsigned rotPivot = count;
	report[NumEdgeInterse].tipo = 1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Dist_Edge_Origin(unsigned int NumArista, unsigned int index)
{
	float Ax, Ay, Bx, By, Cx, Cy, u, Px, Py;
	Vector3D P;

	Ax = v[edg[NumArista].vertex1].centro[0]; Ay = v[edg[NumArista].vertex1].centro[1];
	Bx = v[edg[NumArista].vertex2].centro[0]; By = v[edg[NumArista].vertex2].centro[1];
	Cx = 0; Cy = 0;
	u = ((Cx - Ax) * (Bx - Ax) + (Cy - Ay) * (By - Ay)) / ((Bx - Ax) * (Bx - Ax) + (By - Ay) * (By - Ay));
	//printf("u=%2.2f, ", u);
	//si u esta entre 0 y 1, entonces atraviesa el circulo
	if (u > 0 && u < 1) {
		P.x = Ax + u * (Bx - Ax); 	P.y = Ay + u * (By - Ay); 	P.z = 0;
		//calculando distancia de (0,0) a (Px, Py) 
		dist[index] = P.length();
		//printf("distancia a P desde C es %f \n ", dist[index]);
	}
	else if (u > 1)
	{
		P.vectbuild(Vector3D(0, 0, 0), Vector3D(Bx, By, 0));
		dist[index] = P.length();
		//printf("distancia a C desde B es %f \n ", dist[index]);
	}
	else {
		P.vectbuild(Vector3D(0, 0, 0), Vector3D(Ax, Ay, 0));
		dist[index] = P.length();
		//printf("distancia a C desde A es %f \n ", dist[index]);
	}
	//printf("\n");

}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
float Circle_Edge_Cross_v1(unsigned int NumArista, float radio) {
	float t, dist = 0.0;
	Vector3D A, B, M, P;
	Vector3D temp;

	A = Vector3D(v[edg[NumArista].vertex1].centro[0], v[edg[NumArista].vertex1].centro[1], 0.0);
	B = Vector3D(v[edg[NumArista].vertex2].centro[0], v[edg[NumArista].vertex2].centro[1], 0.0);
	M = B - A;

	//recuerda que P es el centro del circulo en el origen
	temp = P - A;
	t = temp.dotproduct(M) / M.dotproduct(M);
	//printf("M.x=%2.2f, M.y=%2.2f, M.z=%2.2f, t=%2.2f ", M.x, M.y, M.z, t);
	//solo ajusto para que aceptar a los menores o iguales a 0.96
	if (t < 0.99) {
		//A.x = A.x + t * M.x; A.y = A.y + t * M.y; A.z = A.z + t * M.z;
		P.x = P.x - (A.x + t * M.x); P.y = P.y - (A.y + t * M.y); P.z = 0.0;
		dist = P.length();
		//printf(" dist=%2.2f  ", dist);
		//printf("  <----");
		edg[NumArista].collision = 1;
		NumIntersecAro++;

	}
	//printf("\n");
	return dist;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void detecta_colision_con_circulos_interiores() {
	float distance = 0.0, r1;
	unsigned int temp = 0, tempSum = 0;
	NumIntersecAro = 0;
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
		distance = Circle_Edge_Cross_v1(EdgeInterse[i], r1);

		if (edg[EdgeInterse[i]].collision) {
			//printf(". . .edg[%d]  NumIntersecAro = %d, ", EdgeInterse[i], NumIntersecAro);
			BetaCircAdj = NumIntersecAro;
			temp = 0; // NumIntersecAro;
			for (unsigned int j = CapaActual + 2; j < MaxNivel; j++) {
				r1 = (MaxNivel - j) * 5.0 - 4; //radio menor, capa mas profunda
				if (distance < r1) {
					temp = temp + 2;

					//printf(". . . temp=%d, ", temp);
				}

			}//for j
			tempSum = tempSum + temp;
			//printf(". . . tempSum=%d . . . ", tempSum);
		}//arista marcada
		//printf("\n");
	}//for i

	//printf("\n");
	NumIntersecAro = NumIntersecAro + tempSum;
	BetaCirc = tempSum;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Circle_Edge_Cross(unsigned int NumArista, float radio)
{

	float Ax, Ay, Bx, By, Cx, Cy, u, Px, Py;
	Vector3D P, V;

	Ax = v[edg[NumArista].vertex1].centro[0]; Ay = v[edg[NumArista].vertex1].centro[1];
	Bx = v[edg[NumArista].vertex2].centro[0]; By = v[edg[NumArista].vertex2].centro[1];
	Cx = 0; Cy = 0;
	u = ((Cx - Ax) * (Bx - Ax) + (Cy - Ay) * (By - Ay)) / ((Bx - Ax) * (Bx - Ax) + (By - Ay) * (By - Ay));
	printf("u=%2.2f, ", u);
	//si u esta entre 0 y 1, entonces atraviesa el circulo
	if (u > 0 && u < 1) {
		P.x = Ax + u * (Bx - Ax); 	P.y = Ay + u * (By - Ay); 	P.z = 0;
		//calculando distancia de (0,0) a (Px, Py) 
		printf("distancia a P desde C es %f \n ", P.length());
	}
	else if (u > 1)
	{
		V.vectbuild(Vector3D(0, 0, 0), Vector3D(Bx, By, 0));
		printf("distancia a C desde B es %f \n ", V.length());
	}
	else {
		V.vectbuild(Vector3D(0, 0, 0), Vector3D(Ax, Ay, 0));
		printf("distancia a C desde A es %f \n ", V.length());
	}
	printf("\n");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void collision_ultimo_segmento(unsigned int NumArista1, unsigned int NumSegm1, unsigned int report1, unsigned int NumArista2, unsigned int NumSegm2, unsigned int report2)
{
	Vector3D vect0, vect1, vect2, vect3;
	printf("edg[%d] (%d)  ", NumArista1, NumSegm1);
	printf("vs edg[%d] (%d) \n", NumArista2, NumSegm2);

	if (NumSegm1 > 0) { //el ultimo segmento
		vect0 = edg[NumArista1].stop[report[report1].cont_seguimiento - 1];
		vect1.x = v[edg[NumArista1].vertex2].centro[0];
		vect1.y = v[edg[NumArista1].vertex2].centro[1];
		//------------------------------------------
		if (NumSegm2 > 0) {
			for (int h = 0; h < NumSegm2 - 1; h++) { //segmentos de arista j
				printf("   segm h=%d,  ", h);
				vect2 = edg[NumArista2].stop[h];
				vect3 = edg[NumArista2].stop[h + 1];
				//:   :   :   :   :   :   :   :   :
				if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
				{
					report[report1].num_collisions_red++; report[report2].num_collisions_red++;
					printf(" (cross) ");
				}//if
				 //:   :   :   :   :   :   :   :   :
			}//for h, segmentos de arista j
			printf("\n");
			//the last segment
			printf("   Last segm  edg[%d] \n", NumArista2);
			vect2 = edg[NumArista2].stop[report[report2].cont_seguimiento - 1];
			vect3.x = v[edg[NumArista2].vertex2].centro[0];
			vect3.y = v[edg[NumArista2].vertex2].centro[1];
			//:   :   :   :   :   :   :   :   :
			if (edg[NumArista1].vertex2 != edg[NumArista2].vertex2)
				if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
				{
					report[report1].num_collisions_red++; report[report2].num_collisions_red++;
					printf(" (cross) ");
				}//if
				 //the last segment
		}//if (NumSegm2 > 0)
		else {
			printf("   Just the edg[%d] \n", NumArista2);
			vect2.x = v[edg[NumArista2].vertex1].centro[0]; vect2.y = v[edg[NumArista2].vertex1].centro[1];
			vect3.x = v[edg[NumArista2].vertex2].centro[0]; vect3.y = v[edg[NumArista2].vertex2].centro[1];
			//:   :   :   :   :   :   :   :   :
			if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
			{
				report[report1].num_collisions_red++; report[report2].num_collisions_red++;
				printf(" (cross) ");
			}//if
			 //:   :   :   :   :   :   :   :   :
		}//else, if (NumSegm2 > 0)
		 //------------------------------------------
	}//if (NumSegm1 > 0)
	else {

	}	//else, if (NumSegm1 > 0)
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void collision_segmentos(unsigned int NumArista1, unsigned int NumSegm1, unsigned int report1, unsigned int NumArista2, unsigned int NumSegm2, unsigned int report2)
{
	Vector3D vect0, vect1, vect2, vect3;

	printf("edg[%d] (%d)  ", NumArista1, NumSegm1);
	printf("vs edg[%d] (%d) \n", NumArista2, NumSegm2);

	if (NumSegm1 > 0) {
		for (int k = 0; k < NumSegm1 - 1; k++) { //segmentos de arista i
			printf("   segm k=%d \n ", k);
			vect0 = edg[NumArista1].stop[k];
			vect1 = edg[NumArista1].stop[k + 1];
			//------------------------------------------
			if (NumSegm2 > 0) {
				for (int h = 0; h < NumSegm2 - 1; h++) { //segmentos de arista j
					printf("   segm h=%d,  ", h);
					vect2 = edg[NumArista2].stop[h];
					vect3 = edg[NumArista2].stop[h + 1];
					//:   :   :   :   :   :   :   :   :
					if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
					{
						report[report1].num_collisions_red++; report[report2].num_collisions_red++;
						printf(" (cross) ");
					}//if
					 //:   :   :   :   :   :   :   :   :
				}//for h, segmentos de arista j
				printf("\n");
				//the last segment
				printf("   Last segm  edg[%d] \n", NumArista2);
				vect2 = edg[NumArista2].stop[report[report2].cont_seguimiento - 1];
				vect3.x = v[edg[NumArista2].vertex2].centro[0];
				vect3.y = v[edg[NumArista2].vertex2].centro[1];
				//:   :   :   :   :   :   :   :   :
				if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
				{
					report[report1].num_collisions_red++; report[report2].num_collisions_red++;
					printf(" (cross) ");
				}//if
				//the last segment
			}//if (NumSegm2 > 0)
			else {
				printf("   Just the edg[%d] \n", NumArista2);
				vect2.x = v[edg[NumArista2].vertex1].centro[0]; vect2.y = v[edg[NumArista2].vertex1].centro[1];
				vect3.x = v[edg[NumArista2].vertex2].centro[0]; vect3.y = v[edg[NumArista2].vertex2].centro[1];
				//:   :   :   :   :   :   :   :   :
				if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
				{
					report[report1].num_collisions_red++; report[report2].num_collisions_red++;
					printf(" (cross) ");
				}//if
				 //:   :   :   :   :   :   :   :   :
			}//else, if (NumSegm2 > 0)
			//------------------------------------------
		}//for k, segmentos de arista i
	}//if (NumSegm1 > 0)
	else {
		printf("   Just the edg[%d] \n", NumArista1);
		vect0.x = v[edg[NumArista1].vertex1].centro[0]; vect0.y = v[edg[NumArista1].vertex1].centro[1];
		vect1.x = v[edg[NumArista1].vertex2].centro[0]; vect1.y = v[edg[NumArista1].vertex2].centro[1];
		//------------------------------------------
		if (NumSegm2 > 0) {
			for (int h = 0; h < NumSegm2 - 1; h++) { //segmentos de arista j
				printf("   segm h=%d,  ", h);
				vect2 = edg[NumArista2].stop[h];
				vect3 = edg[NumArista2].stop[h + 1];
				//:   :   :   :   :   :   :   :   :
				if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
				{
					report[report1].num_collisions_red++; report[report2].num_collisions_red++;
					printf(" (cross) ");
				}//if
				 //:   :   :   :   :   :   :   :   :
			}//for h, segmentos de arista j
			printf("\n");
		}//if (NumSegm2 > 0)
		else {
			printf("   Just the edg[%d] \n", NumArista2);
			vect2.x = v[edg[NumArista2].vertex1].centro[0]; vect2.y = v[edg[NumArista2].vertex1].centro[1];
			vect3.x = v[edg[NumArista2].vertex2].centro[0]; vect3.y = v[edg[NumArista2].vertex2].centro[1];
			//:   :   :   :   :   :   :   :   :
			if (Intersect_two_vertex(vect0.x, vect0.y, vect1.x, vect1.y, vect2.x, vect2.y, vect3.x, vect3.y))
			{
				report[report1].num_collisions_red++; report[report2].num_collisions_red++;
				printf(" (cross) ");
			}//if
			 //:   :   :   :   :   :   :   :   :
		}//else, if (NumSegm2 > 0)
		//------------------------------------------
	}//else, if (NumSegm1 > 0)
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void colision_entre_segmentos(int NumArista1, int NumArista2) {
	/*int iniciaArista2 = 0, j; //report[0].num_collisions_red
	bool flag = 0;*/

	for (int i = 0; i < edg[NumArista1].cont_fin - 1; i++) {
		for (int j = 0; j < edg[NumArista2].cont_fin - 1; j++) {
			//printf("      segm %d-%d de arista %d vs segm %d-%d de arista %d  ", i, i + 1, NumArista1, j, j + 1, NumArista2);
			if (Intersect_two_vertex(edg[NumArista1].stop[i].x, edg[NumArista1].stop[i].y,
				edg[NumArista1].stop[i + 1].x, edg[NumArista1].stop[i + 1].y,
				edg[NumArista2].stop[j].x, edg[NumArista2].stop[j].y,
				edg[NumArista2].stop[j + 1].x, edg[NumArista2].stop[j + 1].y)) {
				if (edg[NumArista1].stop[i + 1].x == edg[NumArista2].stop[j + 1].x && edg[NumArista1].stop[i + 1].y == edg[NumArista2].stop[j + 1].y)
				{
					//printf("   nodo destino compartido. . . ");
				}
				else {
					NumColision00++;
					//printf(" (collision %d) \n", NumColision00); 
					break;
				}
			}
			//printf("\n");
		}//for j
		//printf("\n");
	}//for i
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void colision_entre_segmentos_negros(int NumArista) {
	for (unsigned int j = 0; j < EdgeIntersecNum; j++) {
		if (!edg[EdgeInterse[j]].collision) {  //si es arista gamma negra
			for (int i = 0; i < edg[NumArista].cont_fin - 1; i++) {
				if (Intersect_two_vertex(edg[NumArista].stop[i].x, edg[NumArista].stop[i].y,
					edg[NumArista].stop[i + 1].x, edg[NumArista].stop[i + 1].y,
					v[edg[EdgeInterse[j]].vertex1].centro[0], v[edg[EdgeInterse[j]].vertex1].centro[1],
					v[edg[EdgeInterse[j]].vertex2].centro[0], v[edg[EdgeInterse[j]].vertex2].centro[1]))
				{
					if (edg[NumArista].stop[i + 1].x == v[edg[EdgeInterse[j]].vertex2].centro[0] && edg[NumArista].stop[i + 1].y == v[edg[EdgeInterse[j]].vertex2].centro[1])
					{
						//printf(". . . nodo destino compartido. . . ");
					}
					else {
						NumColision01++;
						//printf(" (collision %d) \n", NumColision01); 
						break;
					}
				}
			}
		}
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void colision_entre_segmentos_negros_only() {
	for (unsigned int i = 0; i < EdgeIntersecNum - 1; i++)
	{
		if (!edg[EdgeInterse[i]].collision)  //si es arista gamma negra
			for (unsigned int j = i + 1; j <= EdgeIntersecNum - 1; j++)
			{
				if (!edg[EdgeInterse[j]].collision)  //si es arista gamma negra
					if (Intersect(edg[EdgeInterse[i]].vertex1, edg[EdgeInterse[i]].vertex2, edg[EdgeInterse[j]].vertex1, edg[EdgeInterse[j]].vertex2))
					{
						NumColision10++;
						//printf("collision entre aristas edg[EdgeInterse[%d] vs edg[EdgeInterse[%d] NumIntersec=%d\n", i, j, NumColision10);
					}
			}
		//printf("\n");
	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void colision_entre_aristas() {
	//cout << "- - - - - - - -   colision_entre_aristas() \n";
	NumColision00 = 0;
	//cout << "- - - - - - - -   colision_entre_segmentos(   ) \n";
	if(Num_Aristas_Procesar - 1 > 0)
	for (int i = 0; i < Num_Aristas_Procesar - 1; i++)
		for (int j = i + 1; j < Num_Aristas_Procesar; j++) {
			//printf("Checar segmentos de arista %d(%d) vs segmentos de arista %d(%d) - - -  \n", i, EdgeInterse[aristas_a_procesar[i]], j, EdgeInterse[aristas_a_procesar[j]]);
			colision_entre_segmentos(EdgeInterse[aristas_a_procesar[i]], EdgeInterse[aristas_a_procesar[j]]);
		}
	//printf("gammHat,gammaHat=%d \n", NumColision00);
	//para aristas gamma negras
	NumColision01 = 0;
	//cout << "- - - - - - - -   colision_entre_segmentos_negros(   ) \n";
	if (Num_Aristas_Procesar  > 0)
	for (int i = 0; i < Num_Aristas_Procesar; i++) {
		//printf("Checar segmentos de arista %d(%d) vs aristas gamma negras - - -  \n", i, EdgeInterse[aristas_a_procesar[i]]);
		colision_entre_segmentos_negros(EdgeInterse[aristas_a_procesar[i]]);
	}
	//printf("gammHat,gamma=%d \n", NumColision01);

	//colision entre aristas negras only
	//printf("colision entre aristas negras only - - -  \n");
	NumColision10 = 0;
	//cout << "- - - - - - - -   colision_entre_segmentos_negros_only(   ) \n";
	colision_entre_segmentos_negros_only();
	//printf("gamm,gamma=%d \n", NumColision10);

	//printf("recorrido=%d, NumColision00=%d, NumColision01=%d, NumColision10=%d, Num Segmentos: <", recorido_a_realizar, NumColision00, NumColision01, NumColision10);
	/*if (Num_Aristas_Procesar > 0)
	for (int i = 0; i < Num_Aristas_Procesar; i++) {
		printf("edg[ EdgeInterse[aristas_a_procesar[%d]] ]=%d, ", i, edg[EdgeInterse[aristas_a_procesar[i]]].cont_fin-1);
		printf("%d, ", edg[EdgeInterse[aristas_a_procesar[i]]].cont_fin - 1);
	}
	printf(" > \n");*/
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Genera_Reporte()
{
	Vector3D vect, vect1, vect0, vect2, vect3;
	printf("Detecta colisiones con ...aristas  negras... \n");
	for (int i = 0; i < EdgeIntersecNum; i++) {
		printf("\n i=%d. ---------------------------------------------------------\n", i);
		//report[i].num_collisions = 0;
		if (edg[EdgeInterse[i]].collision) {
			if (report[i].cont_seguimiento > 0) {
				printf("Red i=%d. edg[%d]...\n", i, EdgeInterse[i]);
				for (int j = 0; j < report[i].cont_seguimiento - 1; j++) { //segmentos de arista roja
					//printf("segm %d  ", j);
					vect = edg[EdgeInterse[i]].stop[j];
					vect1 = edg[EdgeInterse[i]].stop[j + 1];
					//.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
					for (int k = 0; k < EdgeIntersecNum; k++)
						if (!edg[EdgeInterse[k]].collision) {                          //vs aristas negras
							printf("vs k=%d edg[%d]...", k, EdgeInterse[k]);
							if (Intersect_one_vertex(edg[EdgeInterse[k]].vertex1, edg[EdgeInterse[k]].vertex2,
								vect.x, vect.y, vect1.x, vect1.y))
							{
								report[i].num_collisions++; report[k].num_collisions++;
								printf("report[%d].num_collisions=%d, report[%d].num_collisions=%d\n", i, report[i].num_collisions, k, report[k].num_collisions);
							}
							else printf("  NO ");
						}//for k
					//.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
					//printf("\n");
				}//for j

				 //verifica colisiones con el segmento final
				printf("segm final de stop[%d] a vertex2=%d\n        ", report[i].cont_seguimiento - 1, edg[EdgeInterse[i]].vertex2);
				//.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
				for (int k = 0; k < EdgeIntersecNum; k++)
					if (!edg[EdgeInterse[k]].collision) {                      //vs aristas negras
						printf("vs k=%d edg[%d]...", k, EdgeInterse[k]);
						vect = edg[EdgeInterse[i]].stop[report[i].cont_seguimiento - 1];
						vect1.x = v[edg[EdgeInterse[i]].vertex2].centro[0];
						vect1.y = v[edg[EdgeInterse[i]].vertex2].centro[1];
						if (edg[EdgeInterse[i]].vertex2 != edg[EdgeInterse[k]].vertex2) {
							if (Intersect_one_vertex(edg[EdgeInterse[k]].vertex1, edg[EdgeInterse[k]].vertex2,
								vect.x, vect.y, vect1.x, vect1.y))
							{
								report[i].num_collisions++; report[k].num_collisions++;
								printf("report[%d].num_collisions=%d, report[%d].num_collisions=%d\n", i, report[i].num_collisions, k, report[k].num_collisions);
							}
							else printf("  NO ");
						}
					}//for k
				 //.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
				//printf("\n");
			}//if (report[i].cont_seguimiento > 0)
			else {
				//arista roja sin segmentos
				vect.x = v[edg[EdgeInterse[i]].vertex1].centro[0]; vect.y = v[edg[EdgeInterse[i]].vertex1].centro[1];
				vect1.x = v[edg[EdgeInterse[i]].vertex2].centro[0]; vect1.y = v[edg[EdgeInterse[i]].vertex2].centro[1];
				//.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
				for (int k = 0; k < EdgeIntersecNum; k++)
					if (!edg[EdgeInterse[k]].collision) {              //vs aristas negras
						printf("vs k=%d edg[%d]...", k, EdgeInterse[k]);
						if (Intersect_one_vertex(edg[EdgeInterse[k]].vertex1, edg[EdgeInterse[k]].vertex2,
							vect.x, vect.y, vect1.x, vect1.y))
						{
							report[i].num_collisions++; report[k].num_collisions++;
							printf("report[%d].num_collisions=%d, report[%d].num_collisions=%d\n", i, report[i].num_collisions, k, report[k].num_collisions);
						}
						else printf("  NO ");
					}//for k
					 //.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
					 //No hay segmento final, ni de rojas, ni de negras
					 //.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
			}//else if (report[i].cont_seguimiento > 0)
		}//if (edg[EdgeInterse[i]].collision)
		else {

		}//else if (edg[EdgeInterse[i]].collision)
	}//for i

	//verificando valores
	//for (int i = 0; i < EdgeIntersecNum; i++)
		//printf("report[%d].num_collisions=%d,\n", i, report[i].num_collisions);
	//Colisiones de aristas gamma black - black

	//
	printf("\n = = = = = = = = = = = = = = = = = = = = \n");
	printf("\n Detecta colisiones con ...aristas  rojas... \n");
	for (int i = 0; i < EdgeIntersecNum - 1; i++) {
		if (edg[EdgeInterse[i]].collision)
			for (int j = i + 1; j < EdgeIntersecNum; j++) {
				if (edg[EdgeInterse[j]].collision)
					collision_segmentos(EdgeInterse[i], report[i].cont_seguimiento, i, EdgeInterse[j], report[j].cont_seguimiento, j);
			}//for j
	}//for i
	printf("\n Segmento final... \n");
	for (int i = 0; i < EdgeIntersecNum - 1; i++) {
		if (edg[EdgeInterse[i]].collision)
			for (int j = i + 1; j < EdgeIntersecNum; j++) {
				if (edg[EdgeInterse[j]].collision)
					collision_ultimo_segmento(EdgeInterse[i], report[i].cont_seguimiento, i, EdgeInterse[j], report[j].cont_seguimiento, j);
			}
	}
	//
	printf("NumEdgeInterse=%d \n", NumEdgeInterse);
	//colision con anullos internos

	ReporteGraph("reportes/GDraw.txt");
	/*for (int i = 0; i < EdgeIntersecNum; i++) {
		printf("edg[%d], ", EdgeInterse[i]);
		printf("num_collisions=%d, ", edg[EdgeInterse[i]].num_collisions);
		printf("cont_seguimiento=%d, ", edg[EdgeInterse[i]].cont_seguimiento);
		printf("cont_fin=%d, ", edg[EdgeInterse[i]].cont_fin);
		printf("radio_path=%f, ", edg[EdgeInterse[i]].radio_path);
		printf("end=%d, ", edg[EdgeInterse[i]].end);
		printf("\n");
	}*/
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void inicia_antes_de_recorridos()
{
	printf("\n GiroMin=%d \n", GiroMin);
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) edg[EdgeInterse[i]].collision = 0;
	for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) edg[EdgeInter[i]].collision = 0;
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) { //reporte
		report[i].tipo = 0;
		report[i].num_collisions = 0;
		report[i].num_collisions_red = 0;
		report[i].cont_seguimiento = 0;
		report[i].rotations = 0;
		report[i].radio_path = 0.0f;
		report[i].inner_circles = 0;
	}
	crossing_black_black = 0;

	//colision con anillos internos
	//printf("Estoy en Capas: %d, %d \n", CapaActual, CapaActual + 1);
	for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
		//printf("%d. edg[%d]. ", i, EdgeInterse[i]);
		Dist_Edge_Origin(EdgeInterse[i], i);
		Calcula_vector_origin(EdgeInterse[i]);
		edg[EdgeInterse[i]].cont_seguimiento = 0;
		edg[EdgeInterse[i]].cont_fin = 0;
		edg[EdgeInterse[i]].num_seguidores = 0;
		edg[EdgeInterse[i]].num_collisions = 0;
		edg[EdgeInterse[i]].ext_path = 0;
		Calcula_vector_Arista(EdgeInterse[i]);
	}

	//para el primer caso si es arista gamma roja
	printf("EdgeInterse[%d]=%d.   ", NumEdgeInterse, EdgeInterse[NumEdgeInterse]);
	if (edg[EdgeInterse[0]].collision) Calcula_vector_origin(EdgeInterse[0]);
	sentido_rot = !sentido_rot;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checa_angulos_de_particiones() {
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	Vector3D vect, vect1;
	float xc, yc, i1, incr, xcR, ycR;
	float r1, rpath;
	float angle, angle1;

	r1 = (MaxNivel - CapaActual - 1) * 5.0 - 4; //radio menor, capa mas profunda
	//el vector de la arista correspondiente
	//vect = edg[NumAr].vect;

	double angulo, result;

	i1 = 0.0;
	incr = 2.0 / 360; //360 nodos para la vuelta entera 2 PI, si tuviera mas nodos habria que ajustar
	//i1 = edg[NumAr].AngleLeft * incr;
	rpath = 1.0;
	r1 = r1 + rpath;

	for (int i = 0; i < 10; i++) {
		xc = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
		yc = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);

		vect.z = 0; vect.x = xc; vect.y = yc;
		vect.unitize();
		angle = vect.dotproduct(Vector3D(1, 0, 0));
		angle1 = (double)acos(angle) * 180 / 3.1415;
		//ajuste para cuadrantes II, IV
		if (vect.x >= 0) {
			if (vect.y >= 0) {}
			else { angle1 = 360 - angle1; }
		}
		else {
			if (vect.y >= 0) {}
			else { angle1 = 360 - angle1; }
		}
		printf("angle1= %f grados, ", angle1);

		i1 += (10 + CapaActual) * incr;
		xcR = plane1[0] * 1 * cos(i1 * PI) + (double)plane2[0] * 1 * sin(i1 * PI);
		ycR = plane1[1] * 1 * cos(i1 * PI) + (double)plane2[1] * 1 * sin(i1 * PI);

	}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void graba_particiones_en_seguimiento(float rpath) {
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
	r1 = r1 + rpath;
	//NumParticiones = 360 / (10 + CapaActual);
	temp = 360.0 / (10.0 + CapaActual);
	NumParticiones = temp;
	//para 20OuterPlanar
	if (CapaActual > 12) { NumParticiones = 13; temp = 27.692308; }
	else if (CapaActual > 9) { NumParticiones = 15; temp = 24.0; }
	else if (CapaActual > 4) { NumParticiones = 20; temp = 18.0; }
	else { NumParticiones = 25; temp = 14.4; }

	//para 25OuterPlanar
	/*if (CapaActual > 17) { NumParticiones = 13; temp = 27.692308; }
	else if (CapaActual > 14) { NumParticiones = 15; temp = 24.0; }
	else if (CapaActual > 10) { NumParticiones = 20; temp = 18.0; }
	else if (CapaActual > 5) { NumParticiones = 25; temp = 14.4; }
	else { NumParticiones = 30; temp = 12.0; }*/

	//para 35OuterPlanar
	/*if (CapaActual >30) { NumParticiones = 13; temp = 27.692308; }
	else if (CapaActual > 26) { NumParticiones = 15; temp = 24.0; }
	else if (CapaActual > 20) { NumParticiones = 20; temp = 18.0; }
	else if (CapaActual > 16) { NumParticiones = 25; temp = 14.4; }
	else if (CapaActual > 10) { NumParticiones = 30; temp = 12.0; }
	else if (CapaActual > 4) { NumParticiones = 35; temp = 10.28; }
	else { NumParticiones = 40; temp = 9.0; }*/

	//para 50OuterPlanar
	/*if (CapaActual > 42) { NumParticiones = 13; temp = 27.692308; }
	else if (CapaActual > 39) { NumParticiones = 15; temp = 24.0; }
	else if (CapaActual > 35) { NumParticiones = 20; temp = 18.0; }
	else if (CapaActual > 31) { NumParticiones = 25; temp = 14.4; }
	else if (CapaActual > 25) { NumParticiones = 30; temp = 12.0; }
	else if (CapaActual > 19) { NumParticiones = 35; temp = 10.28; }
	else if (CapaActual > 14) { NumParticiones = 40; temp = 9.0; }
	else if (CapaActual > 7) { NumParticiones = 45; temp = 8.0; }
	else { NumParticiones = 50; temp = 7.2; }*/


	//printf("%d  particiones de %f grados cada una \n", NumParticiones, temp);
	edg[NumAr].cont_seguimiento = 0;

	double param, resultSin, resultCos;
	param = 0; //paticiones
	puntoSegmento = 0;
	for (int i = 0; i < NumParticiones; i++)
	{
		resultSin = sin(param * PI / 180);
		resultCos = cos(param * PI / 180);
		edg[NumAr].seguimiento[puntoSegmento].x = r1 * resultCos;
		edg[NumAr].seguimiento[puntoSegmento].y = r1 * resultSin;
		edg[NumAr].seguimiento[puntoSegmento].z = 0;
		//printf("%d, ", puntoSegmento);
		puntoSegmento++;
		//printf("The sine, cos of %f degrees is %f, %f.\n", param, resultSin, resultCos);
		param = param + temp; //grados
	}
	edg[NumAr].cont_seguimiento = puntoSegmento;
	//printf("\n edg[%d].cont_seguimiento= %d \n", NumAr, edg[NumAr].cont_seguimiento);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void proceso_calcula_crossings_Todo() {
	time_t start_t, end_t;
	double total_time;
	unsigned int NumColisionNew, count = 0, index;
	unsigned int NumColision_Min = 100, recorrido_Min = 0;
	unsigned int NumColision00_Min = 0, NumColision01_Min = 0, NumColision10_Min = 0;
	unsigned int NumSegmentos = 0;
	float r1;

	printf("Creando archivo txt %s ...\n", "reporteOriginal.txt");
	ofstream fout("reporteOriginal.txt");
	count = 0;
	fout << "Layers " << "\tg,g" << "\tg,b" << "\tg,bb" << "\ttime" << "\tg,g" << "\tg,b" << "\tg,bb" << "\ttime" << "\tg,g" << "\tg,b" << "\tgH,gH" << "\tgH,g" << "\ttime" << "\tsegm" << endl;
	start_t = clock();
	for (unsigned int h = 0; h < 20; h++) { //hasta 20 para oterplanar20
		NumColision_Min = 100, recorrido_Min = 0; count = 0;
		printf("************************<1>   \n");
		//
		crossing_black_black = 0; Giro = 0;
		//for (Giro = 0; Giro < NodosPorNivel[CapaActual]; Giro++) //-
		//{ //-
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) edg[EdgeInterse[i]].collision = 0;
		for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) edg[EdgeInter[i]].collision = 0;
		//GiraHorario(CapaActual); //-
		ColisionEntreLineasDeDiferenteCapa();
		detecta_colision_con_circulos_interiores();
		NumColisionNew = NumIntersec + NumIntersecAro;

		if (NumColisionNew < MinNumColisions) { GiroMin = Giro; MinNumColisions = NumColisionNew; }//-

	//} //-
		end_t = clock();
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;
		printf("\nElapsed time : %0.3f \n", total_time);
		report[0].num_collisions = NumIntersec; 	report[0].cont_seguimiento = BetaCircAdj;
		report[0].rotations = BetaCirc; 	report[0].t_original = total_time;
		fout << CapaActual << "," << CapaActual + 1 << "\t" << report[0].num_collisions << "\t" << report[0].cont_seguimiento << "\t" << report[0].rotations << "\t" << total_time;
		printf("************************< . >   \n");

		printf("************************< SPACE >   \n");

		end_t = clock();
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;
		printf("\nElapsed time : %0.3f \n", total_time);
		report[0].num_collisions = NumIntersec; 	report[0].cont_seguimiento = BetaCircAdj;//-
		report[0].rotations = BetaCirc; 	report[0].t_original = total_time;//-
		fout << "\t" << report[0].num_collisions << "\t" << report[0].cont_seguimiento << "\t" << report[0].rotations << "(" << GiroMin << ")" << "\t" << total_time;


		recorido_a_realizar = 0;

		//printf("   \n");
		//recorrido minimo
		fout << "\t" << NumColision10_Min << "\t" << 0 << "\t" << NumColision00_Min << "\t" << NumColision01_Min << "\t" << total_time << "\t" << NumSegmentos << endl;
		printf("Decrementa Capa Actual \n ");
		if (CapaActual > 0) CapaActual--;
		printf("CapaActual=%d \n", CapaActual);
		GiroMin = 0; MinNumColisions = 100;
		GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
		GeneraLineasEnCapa(CapaActual + 1);
	}
	fout << endl;
	fout.close();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void save_in_file() {
	printf("Creando archivo txt %s ...\n", "random_report.txt");
	ofstream fout("random_report.txt");
		for (unsigned int h = 0; h < 49; h++) {
		for (unsigned int i = 0; i < randomRow[h].size(); i++) {
			fout << "randomRow[" << h << "].push_back(" << randomRow[h][i] << ");   ";
		}
		fout << "\n";
	}

	fout << endl;
	fout.close();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
void proceso_Todo() {
	time_t start_t, end_t;
	time_t start_c, end_c;
	double total_time = 0.0, total_time1 = 0.0, total_time2 = 0.0, total_timec = 0.0, t_acum = 0.0, t2_acum = 0.0;
	unsigned int NumColisionNew, count = 0, index;
	unsigned int NumColision_Min = 100, recorrido_Min = 0;
	unsigned int NumColision00_Min=0, NumColision01_Min=0, NumColision10_Min=0;
	unsigned int NumSegmentos = 0;
	int indexRandom, MaxRows;
	float r1; int NumSegmentosTotal = 0;
	int NumPathsTotal = 0;
	

	printf("Creando archivo txt %s ...\n", "reporte.txt");
	ofstream fout("reporte.txt");
	count = 0;
	fout << "Layers " << "\tg,g" << "\tg,b" << "\tg,b'" << "\tcrossings" << "\ttime1" <<  "\tg,g" << "\tg,b" << "\tg,b'" << "\tphase1" << "\tg,g" << "\tg',g'" << "\tg',g" << "\tphase2" << "\ttime2" << "\tedges" << "\tpaths" << "\tpathsTot" << "\tsegm" << "\tsegmTot" << endl;

	for (unsigned int h = 0; h < 19; h++) {
		printf("........................................................................................................h=%d ...\n", h);
		NumColision_Min = 100, recorrido_Min = 0; count = 0;
		NumSegmentosTotal = 0; NumPathsTotal = 0;
		start_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		crossing_black_black = 0; Giro = 0; GiroMin = Giro;
		//

		//
		report[0].num_collisions = NumIntersec; 	report[0].cont_seguimiento = BetaCircAdj;
		report[0].rotations = BetaCirc;
		fout << CapaActual << "-" << CapaActual + 1 << "\t" << report[0].num_collisions << "\t" << report[0].cont_seguimiento << "\t" << report[0].rotations << "\t" << report[0].num_collisions + report[0].cont_seguimiento + report[0].rotations;

		//-for (Giro = 1; Giro < NodosPorNivel[CapaActual]; Giro++) //-
		//-{ //-
			for (unsigned int i = 0; i < EdgeIntersecNum; i++) edg[EdgeInterse[i]].collision = 0;
			for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) edg[EdgeInter[i]].collision = 0;
			//-GiraHorario(CapaActual); //-
			ColisionEntreLineasDeDiferenteCapa();
			detecta_colision_con_circulos_interiores();
			NumColisionNew = NumIntersec + NumIntersecAro;

			//-if (NumColisionNew < MinNumColisions) { GiroMin = Giro; MinNumColisions = NumColisionNew; }//-

		//-}//-

		end_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;
		total_time1 = total_time;

		report[0].t_original = total_time;

		printf("************************< . >   \n");
		start_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*for (Giro = 0; Giro < GiroMin + 1; Giro++)
		{ 
			GiraHorario(CapaActual); 
		}*/
		end_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;
		t_acum += total_time1 + total_time;
		fout << "\t" << total_time1 + total_time; 


		start_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inicia_antes_de_recorridos();
		ColisionEntreLineasDeDiferenteCapa();
		detecta_colision_con_circulos_interiores();
		end_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;
		report[0].num_collisions = NumIntersec; 	report[0].cont_seguimiento = BetaCircAdj;//-
		report[0].rotations = BetaCirc; 	report[0].t_original = total_time;//-
		fout << "\t" << report[0].num_collisions << "\t" << report[0].cont_seguimiento << "\t" << report[0].rotations << "(" << GiroMin << ")" << "\t" << report[0].num_collisions + report[0].cont_seguimiento + report[0].rotations;// << "\t" << total_time;

		start_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Num_Aristas_Procesar = 0; r1 = 0.0;
		for (unsigned int i = 0; i < EdgeIntersecNum; i++) {
			if (edg[EdgeInterse[i]].collision) //si es arista gamma roja 
			{
				edg[EdgeInterse[i]].radio_path = r1;
				aristas_a_procesar[Num_Aristas_Procesar] = i;
				Num_Aristas_Procesar++;
				r1 = r1 + 0.5;
			}//if
		}//for

		count = pow(2, Num_Aristas_Procesar);
		for (int j = 0; j < Num_Aristas_Procesar; j++) {
			sentido_de_recorridos[0][j] = 0;
		}
		//llenas tabla
		for (int i = 1; i < count; i++) {
			for (int j = 0; j < Num_Aristas_Procesar; j++) {
				index = pow(2, j);
				if (i % index == 0) sentido_de_recorridos[i][j] = !sentido_de_recorridos[i - 1][j];
				else  sentido_de_recorridos[i][j] = sentido_de_recorridos[i - 1][j];
			}
		}
		recorido_a_realizar = 0; //para contabilizar recorridos, primer renglon para todas las aristas
		end_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_time = (double)(end_t - start_t) / (double)CLK_TCK;


		start_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_timec = 0.0;
		NumSegmentos = 0;
		MaxRows = pow(2, Num_Aristas_Procesar);//renglones
		printf("\n Num_Aristas_Procesar=%d, MaxRows =%d count=%d\n", Num_Aristas_Procesar, MaxRows, count);
		//printf("NumSegmentos: ", NumSegmentos);
		while (recorido_a_realizar < pow(2, Num_Aristas_Procesar)) {
		//-while (recorido_a_realizar < Num_Aristas_Procesar) {//escogere n renglones
			r1 = 1.0;
			if (recorido_a_realizar < randomRow[CapaActual].size()) indexRandom = randomRow[CapaActual][recorido_a_realizar];
			//-else indexRandom = rand() % MaxRows; //aleatoriamente un renglon
			//**indexRandom = rand() % MaxRows; //**
			//**randomRow[CapaActual].push_back(indexRandom);//**
			//randomRow[CapaActual][recorido_a_realizar]= indexRandom;

			//printf("randomRow[%d][%d]=%d, ", CapaActual, recorido_a_realizar, indexRandom);//**
			
			for (int j = 0; j < Num_Aristas_Procesar; j++) { //por columna, cada arista
				graba_particiones_en_seguimiento(r1);
				sentido_rot = sentido_de_recorridos[recorido_a_realizar][j];
				//-sentido_rot = sentido_de_recorridos[indexRandom][j];
				//
				NumEdgeInterse = aristas_a_procesar[j];
				Calcula_vector_origin(EdgeInterse[NumEdgeInterse]);
				/*if (CapaActual == 4 && j == 0) {
					printf(". . . Checking... recorido_a_realizar=%d\n", recorido_a_realizar); Checa_Capa1(4);
					imprime = 1;
				}*/
				Start_Running();
				
				NumSegmentos += edg[EdgeInterse[NumEdgeInterse]].cont_fin - 1;
				//despues de hacer el recorrido de la arista 0, hay que extender los puntos de cont_seguimiento[ ]
				r1 = r1 + 0.5;
				NumPathsTotal++;
			}
			
			//printf("%d, ", NumSegmentos);
			start_c = clock();////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			colision_entre_aristas(); //NumColision00, NumColision01, NumColision10
			end_c = clock();////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			total_timec += (double)(end_c - start_c) / (double)CLK_TCK;  //t5 en excel

			if (NumColision00 + NumColision01 + NumColision10 < NumColision_Min) {
				NumColision_Min = NumColision00 + NumColision01 + NumColision10;
				recorrido_Min = recorido_a_realizar;
				//-recorrido_Min = indexRandom;
				NumColision00_Min = NumColision00; NumColision01_Min = NumColision01;
				NumColision10_Min = NumColision10;
			}
			checaCircles = 1;
			recorido_a_realizar++;
		}
		//while
		//if (CapaActual < 7) Checa_Capa1(4);
		recorido_a_realizar = 0;
		end_t = clock(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		total_time2 = (double)(end_t - start_t) / (double)CLK_TCK; //t6 en excel
		
		NumSegmentosTotal += NumSegmentos;
		//printf("\nNumSegmentosTotal=%d\n", NumSegmentosTotal);
		//llamo nuevamente al recorrido Min elegido para que se grafique
		recorido_a_realizar = recorrido_Min;
		r1 = 1.0; NumSegmentos = 0;
		for (int j = 0; j < Num_Aristas_Procesar; j++) {
			graba_particiones_en_seguimiento(r1);
			sentido_rot = sentido_de_recorridos[recorido_a_realizar][j];
			NumEdgeInterse = aristas_a_procesar[j];
			Calcula_vector_origin(EdgeInterse[NumEdgeInterse]);
			Start_Running();
			NumSegmentos += edg[EdgeInterse[NumEdgeInterse]].cont_fin - 1;
			//despues de hacer el recorrido de la arista 0, hay que extender los puntos de cont_seguimiento[ ]
			r1 = r1 + 0.5;
		}

		//recorrido minimo

		t2_acum += total_timec + total_time2;
		fout << "\t" << NumColision10_Min << "\t" << NumColision00_Min << "\t" << NumColision01_Min << "\t" << NumColision10_Min + NumColision00_Min + NumColision01_Min << "\t" << total_timec + total_time2; // << "\t" << NumSegmentos;
		fout << "\t" << EdgeIntersecNum << "\t" << Num_Aristas_Procesar << "\t" << NumPathsTotal << "\t" << NumSegmentos << "\t" << NumSegmentosTotal <<endl;

		if (CapaActual > 0) 	CapaActual--;

			printf("CapaActual=%d \n", CapaActual);
			GiroMin = 0; MinNumColisions = 100;
			//if(CapaActual + 1  == 4) Checa_Capa1(CapaActual + 1);
			
			GeneraLineasEntreCapas(CapaActual + 1, CapaActual);
			GeneraLineasEnCapa(CapaActual + 1);
			//
			for (unsigned int i = 0; i < EdgeIntersecNum; i++) edg[EdgeInterse[i]].collision = 0;
			for (unsigned int i = 0; i < EdgeIntersecNumAro; i++) edg[EdgeInter[i]].collision = 0;
			ColisionEntreLineasDeDiferenteCapa();
			detecta_colision_con_circulos_interiores();

			MinNumColisions = NumIntersec + NumIntersecAro;
			NumEdgeInterse = 0; aristasBGcollision = 0;

	}
	fout << endl;
	fout.close();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::
// 