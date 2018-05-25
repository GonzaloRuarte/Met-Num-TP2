#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <tuple>
#include <stdlib.h>

#include "util.h"
#include "clasificador.h"

#include <cmath>
#include <chrono>

using namespace std;

int main(int argc, char * argv[]) {
    string metodo, nombreTrainSet, nombreTestSet, nombreClassif;
    int kdeKnn = 1;
    int alpha = 31;

    if (!obtenerParametros(argc, argv, &metodo, &nombreTrainSet, &nombreTestSet, &nombreClassif)) {
        cout << "Modo de uso: tp2 -m <method> -i <train_set> -q <test_set> -o <classif>\n";
    } else {
        vector<vector<double>> *trainSet = new vector<vector<double>>;
        vector<uint> *labelsTrainSet = new vector<uint>;
        vector<vector<double>> *testSet = new vector<vector<double>>;
        vector<uint> clasificacion;

        cargarSet(nombreTrainSet, trainSet, labelsTrainSet);
        cargarSet(nombreTestSet,  testSet);

        int met = stoi(metodo);

        switch(met) {
            case 0: // kNN
                clasificacion = vectorDeKnns(*trainSet, *labelsTrainSet, *testSet, kdeKnn);
                break;
            case 1: // PCA + kNN
                vector<vector<double>> V = PCATecho(*trainSet, alpha);
                vector<vector<double> > trainSetTemp = multMat(*trainSet, V);
                vector<vector<double> > testSetTemp = multMat(*testSet, V);
                clasificacion = vectorDeKnns(trainSetTemp, *labelsTrainSet, testSetTemp, kdeKnn);
                break;
                //default:
        }
        guardarClasificacion(nombreClassif, clasificacion);
        delete trainSet;
        delete labelsTrainSet;
        delete testSet;

    }
    return 0;
}

