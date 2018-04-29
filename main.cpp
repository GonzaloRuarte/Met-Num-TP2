#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <algorithm>
#include "rdtsc.h"

using namespace std;

vector<double> calcularMedias(vector<vector<double> > imgs) {
	vector<double> res;
	double acum = 0;
	for (uint i = 0; i < 112*92; i++) {//itero sobre la cantidad de variables (cantidad de pixeles de cada imagen)
		for (uint j = 0; i < imgs.size(); j++){//itero sobre la cantidad de imagenes (cantidad de muestras de cada variable)
			acum += imgs[j][i];//acumulo todos los valores
		}
		res[i] = acum/imgs.size();//divido por la cantidad de imagenes para obtener la media de cada variable y la guardo en el lugar correspondiente de res
		acum = 0;
		
	}
	return res;
}




int main(int argc, char * argv[]) {

    if (argc != 3) {
        cout << "Modo de uso: tp2 archivo p\n";
    } else {
        string nombreArchivo = argv[1];

       
    }

    return 0;
}
