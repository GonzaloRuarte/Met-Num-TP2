//
// Created by christian on 25/05/18.
//

#ifndef PAGERANK_CLASIFICADOR_H
#define PAGERANK_CLASIFICADOR_H


#include <stdlib.h>
#include <vector>
#include <iostream>
#include <map>
#include <random>
#include <queue>
#include <tuple>
#include <cmath>
#include <functional>

#include "util.h"

using namespace std;

vector<vector<double> > PCATecho (vector<vector<double> > trainX, uint alpha);
vector<vector<double> > multMat(const vector<vector<double> >& mat1, const vector<vector<double> >& mat2);
void multMatEsc(vector<vector<double> > &mat, double escalar); //para multiplicar una matriz por un escalar. Afecta a la matriz parámetro.
vector<uint> vectorDeKnns (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<vector<double> >& testY, uint k);
clase_t Knn (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<double>& newImg, uint k);
clase_t Knn_sim (const vector<pair<double,clase_t> >& vecNormas, uint k);
vector<pair<double,clase_t> > vector_de_distancias(const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<double>& newImg);


#endif //PAGERANK_CLASIFICADOR_H
