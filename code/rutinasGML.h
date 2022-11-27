//#pragma once
void GMLread(const char *p_filename)
{
	using namespace std;
	fstream f;
	char c;
	char str[20], labelValue[20], label[20];
	int count = 0;
	int value, idValue, idValue2, idValue3, nivel, colorValue = 0;
	unsigned int bien = 0, i, j;

	NumNodes = 0; NumEdges = 0;
	f.open(p_filename, ios::in);
	printf("Openning %s\n", p_filename);
	while (!f.eof())
	{
		f >> str;
		//printf(": %s ", str); 
		//getch();
		if (bien == 0) //comenzando a leer el archivo
		{
			if (strcmp(str, "graph") == 0)
			{
				//printf("Al primer intento tengo a %s =\n", str); 
				bien++;
				//getch();
			}
			else
			{
				printf("%s + ", str);
				while (strcmp(str, "graph") != 0)
				{
					f >> str;
					printf("%s, ", str);
					//getch();
				}
				printf("\n despues del while se tiene a %s \n", str); bien++;
			}
		}
		else
		{
			if (bien == 1)
				if (strcmp(str, "node") == 0)
				{
					//printf("\n efectivamente: %s \n", str);
					//-printf("%s", str);
					f >> str;
					if (strcmp(str, "[") == 0)
					{
						//-printf("%s", str);
						while (strcmp(str, "]") != 0)
						{
							f >> str;
							if (strcmp(str, "id") == 0)
							{
								f >> idValue;
								//-printf("%d", idValue);
							}//if
							if (strcmp(str, "label") == 0)
							{
								f >> labelValue;
							}//if
							if (strcmp(str, "level") == 0)
							{
								f >> nivel;
								if (nivel > MaxNivel) MaxNivel = nivel;
								NodosPorNivel[nivel - 1]++;
							}//if
							if (strcmp(str, "color") == 0)
							{
								f >> colorValue;
							}//if
						}//while
						 //-printf("%s = ", str); //]
						j = 0; //strcpy(label, "\0");
						strcpy_s(label, sizeof label, "\0");
						//printf("%d, ", strlen(label));
						//--printf("%s, %d", labelValue, nivel);

						v[NumNodes] = Vertex(labelValue, NumNodes, nivel);
						if (colorValue == 1) { v[NumNodes].color[0] = 1; v[NumNodes].color[1] = 0; v[NumNodes].color[2] = 0; }
						else { v[NumNodes].color[0] = 0; v[NumNodes].color[1] = 0; v[NumNodes].color[2] = 0; }
						//printf("v[%d]=(%s, %d, %d)\n", NumNodes, labelValue, NumNodes, nivel);
						PorNivel[nivel - 1].push_back(NumNodes); //*
						NumNodes++;
						//-printf("\n");
					} //if de if(strcmp(str,"[")==0)
				}//if de if(strcmp(str,"node")==0)
				else if (strcmp(str, "edge") == 0)
				{
					//printf("%s", str); //edge
					f >> str;
					if (strcmp(str, "[") == 0)
					{
						//printf("%s", str); //[
						while (strcmp(str, "]") != 0)
						{
							f >> str;
							if (strcmp(str, "source") == 0)
							{
								f >> idValue;
								//printf("%d", idValue);
							}//if
							if (strcmp(str, "target") == 0)
							{
								f >> idValue2;
								//printf("%d", idValue2);
							}//if
							if (strcmp(str, "value") == 0)
							{
								f >> idValue3;
								//printf("%d", idValue3);
							}//if
						}//while
						edg[NumEdges] = Edge(idValue, idValue2);
						//--printf("edg[%d]=(%d, %d)...", NumEdges, idValue, idValue2);
						//relacion con los dos nodos de esta arista
						v[idValue].liga.push_back(NumEdges);
						v[idValue2].liga.push_back(NumEdges);
						NumEdges++;
					}//if(strcmp(str,"[")==0)
				}//else if(strcmp(str,"edge")==0)
		}//else de if(bien==0)

	}//while
	f.close();
	printf("\n - - NumNodes=%d, NumEdges=%d, MaxNivel=%d\n", NumNodes, NumEdges, MaxNivel);
	//getchar();
}


void GMLwrite(char *p_filename)
{
	printf("Creando archivo GML %s\n", p_filename);
	using namespace std;
	//int a;
	ofstream fout(p_filename);
	fout << "graph " << endl;
	fout << "[" << endl;
	for (int i = 0; i < NumNodes; i++)
	{
		fout << "node " << endl;
		fout << "[" << endl;
		fout << "id " << i << endl;
		fout << "]" << endl;
	}
	for (unsigned int i = 0; i < NumEdges; i++)
	{
		fout << "edge " << endl;
		fout << "[" << endl;
		fout << "source " << edg[i].vertex1 << endl;
		fout << "target " << edg[i].vertex2 << endl;
		fout << "]" << endl;
	}
	fout << "]" << endl;
	fout.close();
}

void AdjustOuterGraph(unsigned int comienzaNivel, unsigned int terminaNivel)
{
	unsigned short Num = 0, circunf = 1; //para el nivel 1, circunf=9
	float radio, i1;
	double incr;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };

	int k = terminaNivel;
	while (k >= 0)
	{
		i1 = 0.0;
		incr = 2.0 / NodosPorNivel[k];

		radio = 1.0*(circunf);
		for (std::vector<int>::iterator it = PorNivel[k].begin(); it != PorNivel[k].end(); ++it) //para el nivel k
		{
			
			v[*it].centro[0] = plane1[0] * radio*cos(i1*PI) + plane2[0] * radio*sin(i1*PI);
			v[*it].centro[1] = plane1[1] * radio*cos(i1*PI) + plane2[1] * radio*sin(i1*PI);
			//printf("v[%d].centro[%d]=%f \n", *it, 0, v[*it].centro[0]);
			i1 += incr;
		}
		circunf=circunf+5;
		//printf("k=%d, radio=%2.2f, NodosPorNivel[%d]=%d \n", k, radio, k, NodosPorNivel[k]);
		k = k - 1;
	}//for k
}

void AdjustOuterGraphSize(unsigned int Nivel, float Size)
{
	for (std::vector<int>::iterator it = PorNivel[Nivel].begin(); it != PorNivel[Nivel].end(); ++it)
	{
		v[*it].centro[0] = Size * v[*it].centro[0];
		v[*it].centro[1] = Size * v[*it].centro[1];
	}

}

void AdjustOuterGraphOneLayer(unsigned int Nivel, unsigned int Size)
{
	unsigned short Num = 0, circunf = Size;
	float radio, i1;
	double incr;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	int k = Nivel;
	//no requiero el while, ya que es una sola capa
	i1 = 0.0;
	incr = 2.0 / NodosPorNivel[k];
	radio = 0.1*circunf;
	for (std::vector<int>::iterator it = PorNivel[k].begin(); it != PorNivel[k].end(); ++it)
	{
		v[*it].centro[0] = plane1[0] * radio*cos(i1*PI) + plane2[0] * radio*sin(i1*PI);
		v[*it].centro[1] = plane1[1] * radio*cos(i1*PI) + plane2[1] * radio*sin(i1*PI);
		i1 += incr;
	}
	//circunf++;
}

void ScaleRingLayer(unsigned int Nivel, float Size)
{
	int k = Nivel;
	for (std::vector<int>::iterator it = PorNivel[k].begin(); it != PorNivel[k].end(); ++it)
	{
		v[*it].centro[0] = Size*v[*it].centro[0];
		v[*it].centro[1] = Size*v[*it].centro[1];
	}
}
void AdjustOuterGraphOpposite(unsigned int terminaNivel, unsigned int size)
{
	unsigned short Num = 0, circunf = size; // MaxNivel - terminaNivel; //para el nivel 1, circunf=9
	float radio, i1;
	double incr;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };
	int k = terminaNivel;
	//no requiero el while, ya que es una sola capa
	i1 = 0.0;
	incr = 2.0 / NodosPorNivel[k];
	radio = 0.1*(circunf);
	for (std::vector<int>::iterator it = PorNivel[k].begin(); it != PorNivel[k].end(); ++it) //para el nivel k
	{
		v[*it].centro[0] = -1 * (plane1[0] * radio*cos(i1*PI) + plane2[0] * radio*sin(i1*PI));
		v[*it].centro[1] = plane1[1] * radio*cos(i1*PI) + plane2[1] * radio*sin(i1*PI);
		i1 += incr;
	}
	//circunf++;
	//k = k - 1;
}
void AdjustOuterGraphSG(unsigned int NumSubGrafo, unsigned int comienzaNivel, unsigned int terminaNivel)
{
	unsigned short Nume = 0, circunf = 9, k;
	float radio, i1;
	double incr;
	GLfloat plane1[] = { 1.0, 0.0 };
	GLfloat plane2[] = { 0.0, 1.0 };

	printf("\n ===NumSubGrafo= %d, Begin=%d, End=%d \n", NumSubGrafo, NivelBegin[NumSubGrafo], NivelEnd[NumSubGrafo]);
	printf("======NumSubGrafo= %d, comienzaNivel=%d, terminaNivel=%d \n", NumSubGrafo, comienzaNivel, terminaNivel);
	k = comienzaNivel;
	
	Nume = Num[NumSubGrafo];
	//if (comienzaNivel > 0)
	//for (unsigned int i = NivelBegin[NumSubGrafo]-1; i < comienzaNivel; i++) 
	//Num += NodosPorNivel[i];

	//printf("Comenzare en "); printf("vSG[%d][%d].level=%d vs k+1=%d\n", NumSubGrafo, Num, vSG[NumSubGrafo][Num].level, k + 1);

	for (unsigned int i = comienzaNivel; i <= terminaNivel; i++)
	{
		i1 = 0.0;
		incr = 2.0 / NodosPorNivel[i];
		radio = 0.1*circunf;
		printf("radio=%f, NodosPorNivel[%d]=%d.......................\n", radio, i, NodosPorNivel[i]);
		//printf("Num=%d : : : : : : : . : : \n", Num);
		for (unsigned int j = 0; j < NodosPorNivel[i]; j++)
		{
			printf("vSG[%d][%d].level=%d \n", NumSubGrafo, j+Nume, vSG[NumSubGrafo][j+Nume].level);
			vSG[NumSubGrafo][j + Nume].centro[0] = radio*(plane1[0] * radio*cos(i1*3.1415 / 1) + plane2[0] * radio*sin(i1*3.1415 / 1));
			vSG[NumSubGrafo][j + Nume].centro[1] = radio*(plane1[1] * radio*cos(i1*3.1415 / 1) + plane2[1] * radio*sin(i1*3.1415 / 1));
			i1 += incr;
		}
		Nume += NodosPorNivel[i];
		circunf--;
	}

}
