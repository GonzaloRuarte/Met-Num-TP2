#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <algorithm>
#include "rdtsc.h"

using namespace std;

vector<double> calcularMedias(vector<vector<double> > imgs) {
	vector<double> res (112*92) ;
	double acum;
	for (uint i = 0; i < 112*92; i++) {//itero sobre la cantidad de variables (cantidad de pixeles de cada imagen)
		acum = 0;
		for (uint j = 0; j < imgs.size(); j++){//itero sobre la cantidad de imagenes (cantidad de muestras de cada variable)
			acum += imgs[j][i];//acumulo todos los valores
		}
		res[i] = acum/imgs.size();//divido por la cantidad de imagenes para obtener la media de cada variable y la guardo en el lugar correspondiente de res
		
		
	}
	return res;
}

vector<vector<double> > obtenerX(vector<vector<double> > imgs, vector<double> medias){
	vector<vector<double> > res (112*92) ;
	for (uint i = 0; i < 112*92; i++) {
		for (uint j = 0; j < imgs.size(); j++){
			imgs[j][i] -= medias[i];//le resto la media correspondiente a cada variable
		}
	}
	return res;
}


vector<vector<double> > calcularMx (vector<vector<double> > imgs) {
	vector<double> medias = calcularMedias(imgs);
	vector<vector<double> > X = obtenerX(imgs,medias);
	vector<vector<double> > res (112*112*92*92);
	double acum;
	for (uint i = 0; i < 112*92; i++){
		for (uint j = 0; j < 112*92; j++){
			acum = 0;
			for (uint k = 0; k < imgs.size(); k++){
				acum += X[k][i]*X[k][j]; //calculo de la sumatoria de productos para calcular varianza/covarianza
			}
			res[i][j] = acum/(imgs.size()-1);
		}
	}
	return res;
}

vector<vector<double> > trasponer(vector<vector<double> > mat){
	vector<vector<double> > res (mat.size());
	for (uint i = 0; i<mat.size();i++) {
		for (uint j = 0; j<mat.size();j++) {
			res[i][j] = mat[j][i];
		}
	}
	return res;
	
}

pair<double,vector<double> > metodoPotencia(vector<vector<double> > mat) {
	pair<double,vector<double> > res;
	//recordar normalizar autovector
	return res;
}



vector<vector<double> > multMatEsc(vector<vector<double> > mat, double escalar) {//para multiplicar una matriz por un escalar
	vector<vector<double> > res (mat.size());
	for (uint i = 0; i<mat.size();i++) {
		for (uint j = 0; j<mat.size();j++) {
			res[i][j] = mat[i][j]*escalar;
		}
	}
	return res;
}

vector<vector<double> > multVec(vector<double> vec1) {//para generar una matriz a partir de un vector y su traspuesto
	vector<vector<double> > res (vec1.size());
	for (uint i = 0; i<vec1.size();i++) {
		for (uint j = 0; j<vec1.size();j++) {
			res[i][j] = vec1[i]*vec1[j];
		}
	}
	return res;
}



vector<vector<double> > sumMat(vector<vector<double> > mat1, vector<vector<double> > mat2) {//suma de matrices
	vector<vector<double> > res (mat1.size());
	for (uint i = 0; i<mat1.size();i++) {
		for (uint j = 0; j<mat1.size();j++) {
			res[i][j] = mat1[i][j]+mat2[i][j];
		}
	}
	return res;
}

vector< pair<double,vector<double> > > deflacion(vector<vector<double> > mat) {
	vector< pair<double,vector<double> > > res (mat.size());
	pair<double,vector<double> > autovTemp;
	vector<vector<double> > matrizTemp = mat;
	for (uint i = 0; i < mat.size(); i++){
		autovTemp = metodoPotencia(matrizTemp);
		res[i] = autovTemp;
		matrizTemp = sumMat(matrizTemp, multMatEsc(multVec(autovTemp.second),autovTemp.first*(-1)));

	}
	return res;
}

vector<vector<double> > generarP(vector<vector<double> > mat){
	vector<vector<double> > res (mat.size());
	vector< pair<double,vector<double> > > autovectores = deflacion(mat);
	for (uint i = 0; i< mat.size(); i++){
		res[i] = autovectores[i].second;

	}
	return res;

}

vector<double> restaVec(vector<double> vec1, vector<double> vec2) {
	vector<double> res (vec1.size());
	for (uint i = 0; i < vec1.size(); i++){
		res[i] = vec1[i]-vec2[i];
	}
	return res;
}

double norma2(vector<double> vec){//no tomo raiz porque no hace falta en nuestro caso
	double acum = 0;
	for (uint i = 0; i < vec.size(); i++) {
		acum+= vec[i]*vec[i];
	}
	return acum;
}

uint Knn (vector<vector<double> > trainX, vector<uint> labelsX, vector<double> newImg, uint k) {//las labels podrian no ser un int pero lo dejo asi en un principio
	vector<pair<double,uint> > vecNormas (trainX.size());
	vector<pair<double,uint> > sorted (trainX.size());
	double temp;
	for (uint i = 0; i < trainX.size(); i++) {
		temp = norma2(restaVec(trainX[i],newImg));
		vecNormas[i].first = temp;
		vecNormas[i].second = labelsX[i];
	}
	
	for(uint i = 0; i < trainX.size(); i++) {//sort
		double min = 255*255*92*112;
		uint temp;
		for(uint j = 0; j < vecNormas.size(); j++) {
			if(vecNormas[j].first < min) {
				min=vecNormas[j].first;
				temp = j;
			}
		}
		sorted.push_back(vecNormas[temp]);
		vecNormas.erase(vecNormas.begin()+temp);
	}
	pair<uint,int> masRepetido;
	masRepetido.second = 0;
	for (uint i = 0; i < k; i++) {
		int repetidoTemp = 0;
		for (uint j = 0; j < k; j++) {
			if(sorted[j].second == i){
				repetidoTemp++;
			}
		}
		if (repetidoTemp >= masRepetido.second) {
			masRepetido.second = repetidoTemp;
			masRepetido.first = i;
		}
	}
	

	return masRepetido.first;
}


int main(int argc, char * argv[]) {

    if (argc != 3) {
        cout << "Modo de uso: tp2 archivo p\n";
    } else {
        string nombreArchivo = argv[1];

    }

    return 0;
}
