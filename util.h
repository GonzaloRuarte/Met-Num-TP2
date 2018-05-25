#ifndef PPMLOADER_UTIL_H
#define PPMLOADER_UTIL_H

using namespace std;

void cargarDataSetEnMatriz(string pathAlDataSet, vector<vector<double>>* dataSet, vector<uint>* labelsX);
void verificarMatrizAImagen(string pathImagen, int cantidadDeImagenes, int alto, int ancho, vector<vector<double>>* dataSet);
bool obtenerParametros(int argc, char * argv[], string *metodo, string *trainSet, string *testSet, string *classif);
void convertirMatrizAImagen(string pathImagen, int cantidadDeImagenes, vector<vector<double>>* dataSet);
ofstream getFlujo(const string& nombreArchivo);
/*
 * nombreArchivo:   el nombre del test, pero sin el '.in' ni el '.out' el metodo se va a encargar de cargar los datos de ambos.
 * dataSet:         la matriz que contendra las imagenes de entrada.
 * labels:          los labels de las imagenes cargadas en dataSet.
 * autovalores:     los 15 autovalores de mayor magnitud de la matriz de covarianza Mx (v1, v2, ..., v15) ordenados decrecientemente.
 * */
void cargarTest(string nombreArchivo, vector<vector<double>> *dataSet, vector<uint> *labels, vector<double> *autovalores);
string int2stringConCantidadDigitos(int cantidadDigitos, int i);

void cargarTest_aleat(vector<vector<double>> *dataSet, bool repeticion);

void cargarSet(string nombreArchivo, vector<vector<double>> *dataSet, vector<uint> *labels);
void cargarSet(string nombreArchivo, vector<vector<double>> *dataSet);
void guardarClasificacion(string nombreArchivo, vector<uint> &clasificacion);

#endif //PPMLOADER_UTIL_H
