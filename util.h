#ifndef PPMLOADER_UTIL_H
#define PPMLOADER_UTIL_H

using namespace std;

void cargarDataSetEnMatriz(string pathAlDataSet, vector<vector<double>>* dataSet, vector<uint>* labelsX);
void verificarMatrizAImagen(string pathImagen, int cantidadDeImagenes, int alto, int ancho, vector<vector<double>>* dataSet);
bool obtenerParametros(int argc, char * argv[], string *metodo, string *trainSet, string *testSet, string *classif);

#endif //PPMLOADER_UTIL_H
