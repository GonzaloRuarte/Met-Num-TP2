#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <algorithm>
#include "rdtsc.h"
#include "util.h"
#include <cmath>

using namespace std;

vector<double> calcularMedias(const vector<vector<double> > imgs) {
	uint m = imgs[0].size();
	vector<double> res (m) ;
	double acum;
	for (uint i = 0; i < m; i++) {//itero sobre la cantidad de variables (cantidad de pixeles de cada imagen)
		acum = 0;
		for (uint j = 0; j < imgs.size(); j++){//itero sobre la cantidad de imagenes (cantidad de muestras de cada variable)
			acum += imgs[j][i];//acumulo todos los valores
		}
		res[i] = acum/imgs.size();//divido por la cantidad de imagenes para obtener la media de cada variable y la guardo en el lugar correspondiente de res
		
		
	}
	return res;
}

vector<vector<double> > obtenerX(vector<vector<double> > imgs, vector<double> medias){
	vector<vector<double> > res (0);
	uint n = imgs.size();
	uint m = imgs[0].size();
	for (uint j = 0; j < n; j++){
		vector<double> temp (m);
		for (uint i = 0; i < m; i++) {
			temp[i] = imgs [j][i] - medias[i];//le resto la media correspondiente a cada variable
		}
		res.push_back(temp);
	}
	return res;
}

vector<vector<double> > calcularMx (const vector<vector<double> >* imgs) {
	vector<double> medias = calcularMedias(*imgs);
	vector<vector<double> > X = obtenerX(*imgs,medias);
	vector<vector<double> > res (0);
	const size_t& n = imgs->size();
	uint m = (*imgs)[0].size();
	/*double covar_ij;
	for (uint i = 0; i < 112*92; i++){
        covar_ij = 0;
        for (uint k = 0; k < imgs->size(); k++)
            covar_ij += (*imgs)[k][i]*(*imgs)[k][i]; //calculo de la sumatoria de productos para calcular varianza
        res[i][i] = covar_ij/(n-1);
		for (uint j = i+1; j < 112*92; j++){ //como la matriz es simetrica basta calcular la mitad superior.
            covar_ij = 0;
			for (uint k = 0; k < imgs->size(); k++)
                covar_ij += (*imgs)[k][i]*(*imgs)[k][j]; //calculo de la sumatoria de productos para calcular covarianza
			res[i][j] = covar_ij/(n-1);
			res[j][i] = covar_ij/(n-1);
		}
	}*/
	double acum;
	for (uint i = 0; i < m; i++){
		vector<double> temp (m);
		for (uint j = i; j < m; j++){
			acum = 0;
			for (uint k = 0; k < n; k++){
				acum += X[k][i]*X[k][j]; //calculo de la sumatoria de productos para calcular varianza/covarianza
			}
			temp[j] = acum/(n-1);
		}
		res.push_back(temp);
	}
	return res;
}

vector<vector<double> > trasponer(const vector<vector<double> > mat){
	vector<vector<double> > res = mat;
	for (uint i = 0; i<mat.size();i++)
		for (uint j = 0; j<mat.size();j++)
			res[i][j] = mat[j][i];
	return res;
	
}

double norma1(const vector<double> &v){
    double res = 0;
    for(size_t i = 0; i < v.size(); ++i)
        res += abs(v[i]);
    return res;
}

double norma2(const vector<double> &vec){//no tomo raiz porque no hace falta en nuestro caso
    double acum = 0;
    for (uint i = 0; i < vec.size(); i++) {
        acum+= vec[i]*vec[i];
    }
    return acum;
}

vector<double> mult_matr_por_vect(const vector<vector<double> > &M, const vector<double> &v){
    const size_t& n = v.size();
    vector<double> res(n);
    for(size_t i = 0; i < n; ++i) {
        res[i] = 0;
        for (size_t k = 0; k < n; ++k)
            res[i] += M[i][k]*v[k];
    }
    return res;
}

void normalizar1(vector<double>& v){     //Según norma 1
    double norma = norma1(v);
    for(size_t i = 0; i < v.size(); ++i)
        v[i] = v[i]/norma;
}

void normalizar2(vector<double>& v){     //Según norma 2
    double norma = sqrt(norma2(v));
    for(size_t i = 0; i < v.size(); ++i)
        v[i] = v[i]/norma;
}

vector<double> restaVec(const vector<double> &vec1, const vector<double> &vec2) {
    vector<double> res (vec1.size());
    for (uint i = 0; i < vec1.size(); i++)
        res[i] = vec1[i]-vec2[i];
    return res;
}

pair<double,vector<double> > metodoPotencia(const vector<vector<double> > &M) {
    const size_t& n = M[0].size();
	pair<double,vector<double> > res;
	double& autovalor = res.first;
    double autovalor_temp;
    vector<double>& autovector = res.second;
    autovector = vector<double>(n,0);
    autovector[0] = 1;   //Empieza siendo e1.
    double norma = 1;
    vector<double> autovector_temp;
    double norma_temp;
    //Cálculo del autovalor:
    double diferencia = 1;
    while(diferencia >= 0.000001){
        autovector_temp = mult_matr_por_vect(M, autovector);
        norma_temp = norma1(autovector_temp);
        autovalor_temp = norma_temp/norma;
        autovector = mult_matr_por_vect(M, autovector_temp);
        norma = norma1(autovector);
        autovalor = norma/norma_temp;
        diferencia = abs(autovalor - autovalor_temp);
    }
    size_t i = 0;
    while(autovector[i] == 0 || autovector_temp[i] == 0)
        ++i;
    if(abs(autovector[i] + autovector_temp[i]) < abs(autovector[i] - autovector_temp[i]))   //Si tienen signos distintos...
        autovalor *= (-1);                                                                  //el autovalor es negativo.
    //Cálculo del autovector:
    diferencia = 1;
    while(diferencia >= 0.001){
        autovector_temp = mult_matr_por_vect(M, autovector);
        autovector = mult_matr_por_vect(M, autovector_temp);
        if(autovalor < 0)       //Si el autovalor es negativo los vectores están en sentidos opuestos.
            autovector = mult_matr_por_vect(M, autovector);    //Ahora están en el mismo sentido.
        normalizar1(autovector);
        normalizar1(autovector_temp);
        diferencia = norma1(restaVec(autovector,autovector_temp));
    }
    normalizar2(autovector);
	return res;
}


vector<vector<double> > multMatEsc(const vector<vector<double> > &mat, double escalar) {//para multiplicar una matriz por un escalar
	vector<vector<double> > res = mat;
	for (uint i = 0; i<mat.size();i++) {
		for (uint j = 0; j<mat.size();j++) {
			res[i][j] = mat[i][j]*escalar;
		}
	}
	return res;
}

vector<vector<double> > multVec(const vector<double> &vec1) {//para generar una matriz a partir de un vector y su traspuesto
	vector<vector<double> > res (0); 
	uint n = vec1.size();
	for (uint i = 0; i<n;i++) {
		vector<double> temp (n);
		for (uint j = 0; j<n;j++) {
			temp[j] = vec1[i]*vec1[j];
		}
		res.push_back(temp);
	}
	return res;
}


vector<vector<double> > sumMat(const vector<vector<double> > &mat1, const vector<vector<double> > &mat2) {//suma de matrices
	vector<vector<double> > res = mat1;
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
	for (uint i = 0; i < mat.size(); i++){
		autovTemp = metodoPotencia(mat);
		res[i] = autovTemp;
		mat = sumMat(mat, multMatEsc(multVec(autovTemp.second),autovTemp.first*(-1)));
	}
	return res;
}

vector<vector<double> > generarP(const vector<vector<double> > &mat){
	vector<vector<double> > res (mat.size());
	vector< pair<double,vector<double> > > autovectores = deflacion(mat);
	for (uint i = 0; i< mat.size(); i++){
		res[i] = autovectores[i].second;
	}
	return res;
}


uint Knn (vector<vector<double> > trainX, vector<uint> labelsX, vector<double> newImg, uint k) {//las labels podrian no ser un uint pero lo dejo asi en un principio
	vector<pair<double,uint> > vecNormas (trainX.size());
	vector<pair<double,uint> > sorted (0);
	double temp;
	for (uint i = 0; i < labelsX.size(); i++) {
		for(uint j = 0; j<10; j++){
			temp = norma2(restaVec(trainX[i*10+j],newImg));
			vecNormas[i*10+j].first = temp;
			vecNormas[i*10+j].second = labelsX[i];
		}
	}
	for(uint i = 0; i < k; i++) {//sort de menor a mayor segun las normas
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
	for (uint i = 0; i < labelsX.size(); i++) {//calculo del mas repetido de los k vecinos mas cercanos
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

vector<vector<double> > multMat( vector<vector<double> > mat1, vector<vector<double> > mat2) {
	vector<vector<double> > res (0);
	uint n = mat1.size();
	uint m = mat2[0].size();
	for (uint i = 0; i < mat1.size(); i++){
		vector<double> temp (n);
		for (uint j = 0; j < m; j++){
			double acum = 0;
			for (uint k = 0; k < n; k++){
				acum+= mat1[i][k]+mat2[k][j];
			}
			temp[j] = acum;
		}
		res.push_back(temp);	

	}
	return res;
}

vector<vector<double> > PCA (vector<vector<double> >* trainX, uint alpha) {
	cout << 1 << endl;
	vector<vector<double> > Mx = calcularMx(trainX);
	cout << 1 << endl;
	vector<vector<double> > V = trasponer(generarP(Mx));
	cout << 1 << endl;
	for (uint i = 0; i < 92*112; i++){
		V[i].erase(V[i].begin()+alpha, V[i].end());
	}
	return multMat(*trainX,V);
}

int main(int argc, char * argv[]) {
    string metodo, trainSet, testSet, classif;

/*    if (!obtenerParametros(argc, argv, &metodo, &trainSet, &testSet, &classif)) {
        cout << "Modo de uso: tp2 -m <method> -i <train_set> -q <test_set> -o <classif>\n";
    } else {*/
		vector<vector<double>>* dataSet = new vector<vector<double> >;
		vector<uint>* labelsX = new vector<uint > (41);
		

		cargarDataSetEnMatriz("./ImagenesCaras",dataSet, labelsX);
		uint x = Knn(*dataSet,*labelsX,(*dataSet)[71],1);
		vector<vector<double>> asd = PCA(dataSet,6);
		delete labelsX;
		delete dataSet;

    /*}*/

    return 0;
}
