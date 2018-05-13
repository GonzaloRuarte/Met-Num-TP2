#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <tuple>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "rdtsc.h"
#include "util.h"
#include <cmath>
#include <chrono>

using namespace std;

typedef uint clase_t;

typedef struct {
	clase_t clase;
	double precision;
	double recall;
	double f1;
} resultados;


vector<double> calcularMedias(const vector<vector<double> >& imgs) {
    const unsigned long& m = imgs[0].size();
	vector<double> res (m) ;
	double acum;
	for (uint j = 0; j < m; j++) {//itero sobre la cantidad de variables (cantidad de pixeles de cada imagen)
		acum = 0;
		for (uint i = 0; i < imgs.size(); i++){//itero sobre la cantidad de imagenes (cantidad de muestras de cada variable)
			acum += imgs[i][j];//acumulo todos los valores
		}
		res[j] = acum/imgs.size();//divido por la cantidad de imagenes para obtener la media de cada variable y la guardo en el lugar correspondiente de res
	}
	return res;
}

vector<vector<double> > obtenerX(const vector<vector<double> >& imgs, const vector<double>& medias){
    const unsigned long& n = imgs.size();
    const unsigned long& m = imgs[0].size();
	vector<vector<double> > res (n, vector<double>(m));
	for (uint i = 0; i < n; i++){
		for (uint j = 0; j < m; j++){
			res[i][j] = imgs [i][j] - medias[j];//le resto la media correspondiente a cada variable
		}
	}
	return res;
}

vector<vector<double> > obtenerXt(const vector<vector<double> >& imgs, const vector<double>& medias){
    const unsigned long& n = imgs.size();
    const unsigned long& m = imgs[0].size();
    vector<vector<double> > res (m, vector<double>(n));
    for (uint i = 0; i < m; i++){
        for (uint j = 0; j < n; j++){
            res[i][j] = imgs [j][i] - medias[i];//le resto la media correspondiente a cada variable
        }
    }
    return res;
}



vector<vector<double> > calcularMx (const vector<vector<double> >& imgs) {
	vector<double> medias = calcularMedias(imgs);
	vector<vector<double> > X = obtenerX(imgs,medias);
	const size_t& n = imgs.size();
	const size_t& m = imgs[0].size();
	vector<vector<double> > res(m, vector<double>(m));
	double covar_ij;
	double& var_i = covar_ij;
	for (uint i = 0; i < m; i++){
		var_i = 0;
		for (uint k = 0; k < n; k++){
			var_i += X[k][i]*X[k][i]; //calculo de la sumatoria de productos para calcular varianza
		}
		res[i][i] = var_i/(n-1);
		for (uint j = i+1; j < m; j++){ //como la matriz es simetrica basta calcular la mitad superior.
			covar_ij = 0;
			for (uint k = 0; k < n; k++){
				covar_ij += X[k][i]*X[k][j];//calculo de la sumatoria de productos para calcular covarianza
			}	
			res[i][j] = covar_ij/(n-1);
			res[j][i] = covar_ij/(n-1);
		}
	}
/*	double acum;
	for (uint i = 0; i < m; i++){
		for (uint j = i; j < m; j++){
			acum = 0;
			for (uint k = 0; k < n; k++){
				acum += X[k][i]*X[k][j]; //calculo de la sumatoria de productos para calcular varianza/covarianza
			}
			res[i][j] = acum/(n-1);
			res[j][i] = acum/(n-1);
		}
	}*/
	return res;
}

vector<vector<double> > calcularMxTecho (const vector<vector<double> >& X) { //como traspongo X aca el m de aca es el n de la funcion de arriba y viceversa
	const size_t& n = X.size();
	const size_t& m = X[0].size();
	vector<vector<double> > res(m, vector<double>(m));
	double covar_ij;
	double& var_i = covar_ij;
	for (uint i = 0; i < m; i++){
		var_i = 0;
		for (uint k = 0; k < n; k++){
			var_i += X[k][i]*X[k][i]; //calculo de la sumatoria de productos para calcular varianza
		}
		res[i][i] = var_i/(m-1);
		for (uint j = i+1; j < m; j++){ //como la matriz es simetrica basta calcular la mitad superior.
			covar_ij = 0;
			for (uint k = 0; k < n; k++){
				covar_ij += X[k][i]*X[k][j];//calculo de la sumatoria de productos para calcular covarianza
			}	
			res[i][j] = covar_ij/(m-1);
			res[j][i] = covar_ij/(m-1);
		}
	}
	return res;
}



vector<vector<double> > trasponer(const vector<vector<double> >& mat){
    const unsigned long& n = mat.size();
    const unsigned long& m = mat[0].size();
	vector<vector<double> > res (m, vector<double>(n));
	for (uint i = 0; i<n;i++)
		for (uint j = 0; j<m;j++)
			res[j][i] = mat[i][j];
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

double producto_interno(const vector<double> &v, const vector<double> &v2){
    double res = 0;
    if(v.size() != v2.size())
      std::cout << "producto_interno: los vectores no son del mismo tamaño" << std::endl;
    else{
      for(size_t i = 0; i < v.size(); ++i)
        res += v[i]*v2[i];
    }
  return res;
}

double normalizar1(vector<double>& v){     //Según norma 1. Devuelve la norma 1 que tenía el vector.
    double norma = norma1(v);
    for(size_t i = 0; i < v.size(); ++i)
        v[i] = v[i]/norma;
  return norma;
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


vector<double> multVecEsc(const vector<double> &vec1, const double &esc) {
    vector<double> res (vec1.size());
    for (uint i = 0; i < vec1.size(); i++)
        res[i] = vec1[i]*esc;
    return res;
}

pair<double,vector<double> > powMethod(const vector<vector<double> > &M) {
	const size_t& n = M[0].size();
	pair<double,vector<double> > res2;
	uniform_real_distribution<double> unif(0.0,1.0);
	mt19937 re(std::random_device{}());
	auto generator = std::bind(unif, re);
	vector<double> autovector = vector<double>(n);
	vector<double> autovector_temp;
	generate(autovector.begin(), autovector.end(), generator);
	double autovalor = sqrt(norma2(autovector));
		
	while(abs(sqrt(norma2(mult_matr_por_vect(M,autovector)))-sqrt(norma2(multVecEsc(autovector,autovalor)))) >0.00001){
		autovector_temp = mult_matr_por_vect(M,autovector);
		autovalor = norma2(autovector_temp)/norma2(autovector);
		autovector = autovector_temp;
		normalizar2(autovector);
		
	}
	return make_pair(autovalor,autovector);
}

//ofstream salida("Comparación_de_algoritmos.txt", ios_base::out);
pair<double,vector<double> > metodoPotencia(const vector<vector<double> > &M) {
    const size_t& n = M[0].size();
    //algoritmo 1:
/*    pair<double,vector<double> > res;
	double& autovalor = res.first;
    double autovalor_temp;
    vector<double>& autovector = res.second;
    autovector = vector<double>(n,0);
    autovector[0] = 1;   //Empieza siendo e1.
    vector<double> autovector_temp;
    //Cálculo del autovalor:
    double diferencia = 1;
    float cantidad_iteraciones = 0;
    auto t1 = chrono::system_clock::now();
    while(diferencia >= 0.000001){
        autovector_temp = mult_matr_por_vect(M, autovector);
        autovalor_temp = normalizar1(autovector_temp);    //En este paso debo ignorar la normalización
        autovector = mult_matr_por_vect(M, autovector_temp);
        autovalor = normalizar1(autovector);              //Debiese dividir por la norma de autovector_temp, pero es 1.
        diferencia = abs(autovalor - autovalor_temp);
        ++cantidad_iteraciones;
    }
    auto t2 = chrono::system_clock::now();
    size_t i = 0;
    while(autovector[i] == 0 || autovector_temp[i] == 0)
        ++i;
    if(autovector[i] * autovector_temp[i] < 0) {   //Si tienen signos distintos...
        autovalor *= (-1);                         //el autovalor es negativo.
        autovector = mult_matr_por_vect(M, autovector);   //Autovector está en el mismo sentido que autovector_temp.
    }
    //Cálculo del autovector:
    diferencia = norma1(restaVec(autovector,autovector_temp));
    while(diferencia >= 0.001){
        autovector_temp = mult_matr_por_vect(M, autovector);
        normalizar1(autovector_temp);
        autovector = mult_matr_por_vect(M, autovector_temp);
        if(autovalor < 0)       //Si el autovalor es negativo los vectores están en sentidos opuestos.
            autovector = mult_matr_por_vect(M, autovector);    //Ahora están en el mismo sentido.
        normalizar1(autovector);
        diferencia = norma1(restaVec(autovector,autovector_temp));
    }
    normalizar2(autovector);*/
    //algoritmo 2:
    pair<double,vector<double> > res2;
    double& autovalor2 = res2.first;
    double autovalor2_temp;
    vector<double>& autovector2 = res2.second;
    autovector2 = vector<double>(n,0);
    autovector2[0] = 1;   //Empieza siendo e1.
    vector<double> autovector2_temp;
    //Cálculo del autovalor:
    double diferencia2 = 1;
    float cantidad_iteraciones2 = 0;
    //auto t3 = chrono::system_clock::now();
    while(diferencia2 >= 0.00000001){
        autovector2_temp = mult_matr_por_vect(M, autovector2);
        autovalor2_temp = producto_interno(autovector2, autovector2_temp); //autovector está normalizado.
        normalizar2(autovector2_temp);
        autovector2 = mult_matr_por_vect(M, autovector2_temp);
        autovalor2 = producto_interno(autovector2_temp, autovector2); //autovector_temp está normalizado.
        normalizar2(autovector2);
        diferencia2 = abs(autovalor2-autovalor2_temp);
        ++cantidad_iteraciones2;
    }
    //auto t4 = chrono::system_clock::now();
    if(autovalor2 < 0){
        autovector2 = mult_matr_por_vect(M, autovector2);
        normalizar2(autovector2);
        diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        while(diferencia2 > 0.001){
            autovector2_temp = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        }
    }else{
        diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        while(diferencia2 > 0.001){
            autovector2_temp = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2_temp);
            normalizar2(autovector2_temp);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        }
    }
/*    salida.open("Comparación_de_algoritmos.txt", ios::app);
    salida << "Iteraciones: " << cantidad_iteraciones/cantidad_iteraciones2 << ";\t Ciclos de clock: " << double((t2-t1).count())/(t4-t3).count() << endl;
    salida.close();*/
	
    return res2;
}


void multMatEsc(vector<vector<double> > &mat, double escalar) {//para multiplicar una matriz por un escalar. Afecta a la matriz parámetro.
	for (uint i = 0; i<mat.size();i++)
		for (uint j = 0; j<mat[i].size();j++)
			mat[i][j] *= escalar;
}

vector<vector<double> > multVec(const vector<double> &vec1) {//para generar una matriz a partir de un vector y su traspuesto
	const size_t& n = vec1.size();
    vector<vector<double> > res(n, vector<double>(n));
	for (uint i = 0; i<n;i++)
		for (uint j = 0; j<n;j++)
			res[i][j] = vec1[i]*vec1[j];
	return res;
}


void sumMat(vector<vector<double> > &mat1, const vector<vector<double> > &mat2) {//suma de matrices
	for (uint i = 0; i<mat1.size();i++)
		for (uint j = 0; j<mat1.size();j++)
			mat1[i][j] += mat2[i][j];
}

vector< pair<double,vector<double> > > deflacion(vector<vector<double> > mat, uint alpha) {
	vector< pair<double,vector<double> > > res (alpha);
    for (uint i = 0; i < alpha; i++){
        res[i] = metodoPotencia(mat);
        vector<vector<double> > v_x_vt = multVec(res[i].second);    //v*vt
        multMatEsc(v_x_vt,res[i].first*(-1));                       //-lambda_i*(v*vt)
        sumMat(mat, v_x_vt);
	}
	return res;
}

vector<vector<double> > generarV(const vector<vector<double> > &mat, uint alpha){
	vector<vector<double> > res(mat[0].size(), vector<double>(alpha));
	vector< pair<double,vector<double> > > autovectores = deflacion(mat,alpha);
	for (uint j = 0; j < alpha; ++j)
	    for(uint i = 0; i < mat[0].size(); ++i)
		    res[i][j] = autovectores[j].second[i];
	return res;
}


clase_t Knn (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<double>& newImg, uint k) {//las labels podrian no ser un uint pero lo dejo asi en un principio
	vector<pair<double,clase_t> > vecNormas (trainX.size());
    //vector<pair<double,clase_t> > sorted(k);
    /*
    for (uint i = 0; i < labelsX.size(); i++) { //De las 41 labels
        for(uint j = 0; j < 10; j++) { //Leo las 10 imagenes y tomo la distancia.
            vecNormas[i].first  = norma2(restaVec(trainX[i],newImg));
            vecNormas[i].second = labelsX[i];
        }
    }
    make_heap(vecNormas.begin(),vecNormas.end(), greater<pair<double, clase_t> >()); //hago min-heap (para sacar los primeros k minimos.
    for (int i = 0; i < k; i++) {
        sorted[i].first = vecNormas[0].first;
        sorted[i].second = vecNormas[0].second;
        pop_heap(vecNormas.begin(), vecNormas.end(), greater<pair<double, clase_t> >());
        vecNormas.pop_back();
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
    */
	for (uint i = 0; i < labelsX.size(); i++) {
        vecNormas[i].first  = norma2(restaVec(trainX[i],newImg));
		vecNormas[i].second = labelsX[i];
	}
    vector<pair<double,clase_t> >::const_iterator primero = vecNormas.cbegin();
    vector<pair<double,clase_t> >::const_iterator k_esimo = primero+k;
	priority_queue<pair<double,clase_t> > heap(primero, k_esimo);  //Creo un max_heap con los primeros k elementos de vecNormas
	for(size_t i = k; i < vecNormas.size(); ++i){
        if(vecNormas[i] < heap.top()){  //Si el i-ésimo elemento es más chico que el más grande del heap...
            heap.pop();                 //entonces saco al más grande...
            heap.push(vecNormas[i]);    //y meto al i-ésimo elemento.
        }                               //De esta forma me quedo con los k elementos más chicos.
	}
    priority_queue<clase_t> posibles;
	for(int i = k-1; i >= 0; --i){
	    posibles.push(heap.top().second);
	    heap.pop();
	}
/*	for(uint i = 0; i < k; i++) {//sort de menor a mayor segun las normas
		double min = 255*255*92*112;
		clase_t temp;
		for(uint j = 0; j < vecNormas.size(); j++) {
			if(vecNormas[j].first < min) {
				min=vecNormas[j].first;
				temp = j;
			}
		}
		sorted.push_back(vecNormas[temp]);
		vecNormas.erase(vecNormas.begin()+temp);
	}*/
	pair<clase_t,uint> masRepetido;
	masRepetido.second = 0;
    pair<clase_t,uint> comparador;
	while(!posibles.empty()){
        comparador = make_pair(posibles.top(), 1);
        posibles.pop();
	    while(!posibles.empty() && comparador.first == posibles.top()) {
            ++comparador.second;
            posibles.pop();
        }
        if(comparador.second > masRepetido.second)
            masRepetido = comparador;
	}
/*	vector<pair<double, clase_t> > sorted;
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
	}*/
	return masRepetido.first;
}

vector<uint> vectorDeKnns (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<vector<double> >& testY, uint k){
	const unsigned long& n = testY.size();
	vector<uint> res (n);
	for (uint i = 0; i < n; i++) {
		res[i] = Knn(trainX,labelsX,testY[i],k);
	}
	return res;
}


double accuracy (const vector<clase_t>& labelsY, vector<uint>& vectorDeKnns) {
	const unsigned long& n = labelsY.size(); //me di cuenta que si vamos a usar como label el numero de fila dentro de la matriz entonces el labelsY de mucho no nos sirve pero bueno
	double acum = 0;
	for (uint i = 0; i < n; i++) {
		if (vectorDeKnns[i] == labelsY[i]){
			acum++;
		}
	}
	
	
	return acum/n;
}	

double precision(const vector<clase_t>& labelsY, clase_t clase, vector<uint>& vectorDeKnns) {
    const unsigned long& n = labelsY.size();
	double truepositives = 0, positives = 0;
	for (uint i = 0; i < n; i++) {
		if (vectorDeKnns[i] == clase){ //si el Knn dio igual a la clase que estoy procesando entonces sumo 1 a los elementos recuperados
			positives++;
			if (vectorDeKnns[i] == labelsY[i]){ //si el Knn ademas dio bien el resultado sumo 1 a los true positives
				truepositives++;
			}
		}
	}
	double res = 0;
	if (positives>0){res=truepositives/positives;}
	return res;
}

double recall(const vector<clase_t>& labelsY, clase_t clase, vector<uint>& vectorDeKnns) {
    const unsigned long& n = labelsY.size();
	double truepositives = 0, relevants = 0;
	for (uint i = 0; i < n; i++) {
		if (labelsY[i] == clase){//si el elemento que estoy analizando pertenece a la clase que estoy procesando, sumo 1 a los relevantes
		 	relevants++;
			if (vectorDeKnns[i] == clase){ //si el Knn ademas dio bien el resultado sumo 1 a los true positives
				truepositives++;
			}
		}
	}
	return truepositives/relevants;
}

vector<vector<double> > multMat(const vector<vector<double> >& mat1, const vector<vector<double> >& mat2) {
    const unsigned long& n = mat1.size();
    const unsigned long& m = mat2[0].size();
    const unsigned long& l = mat2.size();
	vector<vector<double> > res (n, vector<double>(m, 0));
	for (uint i = 0; i < n; i++)
		for (uint j = 0; j < m; j++)
			for (uint k = 0; k < l; k++)
				res[i][j] += mat1[i][k]*mat2[k][j];
	return res;
}

vector<vector<double> > PCA (vector<vector<double> > trainX, uint alpha) {
    const unsigned long& m = trainX[0].size();
	vector<vector<double> > Mx = calcularMx(trainX);
	vector<vector<double> > V = generarV(Mx,alpha);
	/*vector<vector<double> > H = trasponer(V); //codigo para verificar las componentes principales
	for (uint i = 0; i < H.size(); i++){
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] < 0){
				H[i][j] = -H[i][j];
			}
		}
		double min = H[i][0];
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] < min){
				min = H[i][j];
			}
		}
		for (uint j = 0; j < H[0].size(); j++){
			H[i][j] += min;
		}
		double max = H[i][0];
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] > max){
				max = H[i][j];
			}
		}
		for (uint j = 0; j < H[0].size(); j++){
			H[i][j] *= 255/max;
		}
	}
	convertirMatrizAImagen("./salidaVtraspuesta", alpha, &H);*/
	/*for (uint i = 0; i < m; i++){ //esto no hace falta por ahora porque la V se calcula ya con alpha columnas
		V[i].erase(V[i].begin()+alpha, V[i].end());
	}*/
	return V; //devuelvo la V, recordar multiplicar fuera de la funcion
}

vector<vector<double> > PCATecho (vector<vector<double> > trainX, uint alpha) {
    const unsigned long& m = trainX[0].size();
	vector<double> medias = calcularMedias(trainX);
	vector<vector<double> > Xt = obtenerXt(trainX,medias);
	vector<vector<double> > Mx = calcularMxTecho(Xt);
	vector<vector<double> > P = generarV(Mx,alpha);
	vector<vector<double> > V = multMat(Xt,P);
	V = trasponer(V);
	for (uint i = 0; i < V.size(); i++){
		normalizar2(V[i]);
	}
	V = trasponer(V);
 	/*vector<vector<double> > H = trasponer(V); //codigo para verificar las componentes principales
	for (uint i = 0; i < H.size(); i++){
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] < 0){
				H[i][j] = -H[i][j];
			}
		}
		double min = H[i][0];
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] < min){
				min = H[i][j];
			}
		}
		for (uint j = 0; j < H[0].size(); j++){
			H[i][j] += min;
		}
		double max = H[i][0];
		for (uint j = 0; j < H[0].size(); j++){
			if (H[i][j] > max){
				max = H[i][j];
			}
		}
		for (uint j = 0; j < H[0].size(); j++){
			H[i][j] *= 255/max;
		}
	}
    	convertirMatrizAImagen("./salidaVtraspuesta", alpha, &H);*/
	/*for (uint i = 0; i < m; i++){ //esto no hace falta por ahora porque la V se calcula ya con alpha columnas
		V[i].erase(V[i].begin()+alpha, V[i].end());
	}*/
	return V; //devuelvo la V, recordar multiplicar fuera de la funcion
}

//------------------------- Escritura de las estadisticas -------------------------//

void escribirEstadisiticas(string nombreArchivo, vector<pair<vector<resultados >,double> > &estadisticas, uint kDekfold,uint alpha, bool varioAlpha) { //si varioAlpha es false es porque estoy variando el kdeKnn
    vector<resultados>* estadistica;
    int i;
	if(varioAlpha){
		i=alpha;
	}else{
		i= 1;
	}
    for (vector<pair<vector<resultados >,double> >::iterator it = estadisticas.begin() ; it != estadisticas.end(); ++it) {
        vector<resultados >& estadistica = it->first;
        string accuracy = to_string(it->second);
        ofstream salida = getFlujo(nombreArchivo + "_" + to_string(i));
        string precision = "";
        string recall = "";
        string f1="";
        for (vector<resultados>::iterator it = estadistica.begin() ; it != estadistica.end(); ++it) {
            precision += to_string(it->precision) + "\t";
            recall += to_string(it->recall) + "\t";
            f1 += to_string(it->f1) + "\t";
        }
        salida << accuracy << endl;
        salida << precision << endl;
        salida << recall << endl;
        salida << f1 << endl;
        //cout << precision << endl;
        salida.close();
	if(varioAlpha){
		i-=20; //le resto lo que fui variando el alpha
	}else{
        	i+=20; //le resto lo que fui variando el kDeKnn
	}
    }


}
//------------------------- Escritura de las estadisticas -------------------------//

vector<pair<vector<resultados >,double> > kFold (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k, uint kdeKnn, uint alpha, bool conPCA) {
//codigo para calcular la cantidad de imagenes por persona suponiendo que las muestras son balanceadas y la cantidad de clases
	uint imagenesPorPersona = 0;
	int imagenesPPparagenerador = 0;
	uint cantidadDeClases = 0;
	uint personasTemp = 0;
	for( uint i = 0; i < labelsX.size(); i++){
		if (labelsX[i] != personasTemp){
			personasTemp = labelsX[i];
			cantidadDeClases++;
		}
		if(personasTemp == 1){
			imagenesPorPersona++;
			imagenesPPparagenerador++;
		}
	}
//***********************************************************************//
	vector<int> folds(imagenesPorPersona); //to store the random numbers
    for(uint i = 0; i < imagenesPorPersona; ++i)
        folds[i] = i;
    random_shuffle(folds.begin(), folds.end());
/*	random_device rd; //seed generator
	mt19937_64 generator{rd()}; //generator initialized with seed from rd
	uniform_int_distribution<> dist{0, imagenesPPparagenerador-1}; //the range is inclusive, so this produces numbers in range [0, 10)
	for(uint i=0; i<imagenesPorPersona; ++i) {  //REPITE NÚMEROS
		folds.push_back( dist(generator) );
	} // la idea es que voy a tener muestras balanceadas, entonces para cada persona voy a tener la misma cantidad de imagenes en test y en train
		// como cada persona tiene 10 imagenes, el k puede ser 1, 2, 5 o 10, k = 1 no tiene mucho sentido*/
	uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
	vector<pair<vector<resultados >,double> > res;
	for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
		vector<vector<double> > trainXTemp;
		vector<vector<double> > testYTemp;
		vector<uint> labelsXTemp;
		vector<uint> labelsYTemp;
		for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
			for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
				uint temp = (j*imagenesPorPersona)+folds[u];
				if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
					testYTemp.push_back(trainX[temp]);
					labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
				} else{ 
					trainXTemp.push_back(trainX[temp]);
					labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
				}

			}
		}
		//tengo armado el train y el test para este fold
		if (conPCA){
			vector<vector<double>> V = PCATecho(trainXTemp,alpha); //92*112 = 10304, puse 6 para probar
			uint size_V = V[0].size();
			vector<pair<vector<resultados >,double> > resVariandoAlphaParaUnFold;
			for(uint h = 0; h <= 5; ++h) {//esto sirve para iterar el alpha (voy borrando columnas de la matriz V dependiendo del h)
				for (uint i = 0; i < V.size(); i++){ //borro las columnas de V que necesito borrar para variar el alpha
					V[i].erase(V[i].begin()+size_V-h, V[i].end());//originalmente va h*20 
				}
				vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
				vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
				vector<resultados > resultadosTemp (0);
				//for(uint y = 0; y < 328; ++y) { //este for seria para variar el kDeKnn
				vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,kdeKnn);
				for (uint j = 1; j <= cantidadDeClases; j++){ //itero sobre las clases
					resultados resXClase;
					resXClase.precision = precision(labelsYTemp,j,vectordeKnns);
					resXClase.recall = recall(labelsYTemp,j,vectordeKnns);
					if (resXClase.precision+resXClase.recall > 0){
						resXClase.f1 = 2.0*resXClase.precision*resXClase.recall/(resXClase.precision+resXClase.recall);
					}
					else {
						resXClase.f1 = 0;
					}
					resXClase.clase = j;
					resultadosTemp.push_back(resXClase);
				}//entonces en el vector la posicion 0 corresponde a la clase 1 y asi sucesivamente
				double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
				
				resVariandoAlphaParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp));
				//} //aca terminaria el for que varia el kDeKnn
			
			}
			escribirEstadisiticas("./Resultados/ResultadosVariandoAlpha", resVariandoAlphaParaUnFold,i,alpha,true); //true es que estoy variando el alpha
		//res.push_back(make_pair(resultadosTemp,accuracyTemp)); //comente esto porque no importa mucho lo que devuelve, solo queremos escribir los resultados en archivos, luego de experimentar hay que cambiar esto
		}
		else{//sin PCA
			vector<pair<vector<resultados >,double> > resVariandoKSinPCAParaUnFold;
			for(uint y = 1; y < kdeKnn; y+=20) { //este for seria para variar el kDeKnn
				vector<resultados > resultadosTemp (0);
				vector<uint> vectordeKnns = vectorDeKnns(trainXTemp,labelsXTemp,testYTemp,y);
				for (uint j = 1; j <= cantidadDeClases; j++){ //itero sobre las clases
					resultados resXClase;
					resXClase.precision = precision(labelsYTemp,j,vectordeKnns);
					resXClase.recall = recall(labelsYTemp,j,vectordeKnns);
					if (resXClase.precision+resXClase.recall > 0){
						resXClase.f1 = 2.0*resXClase.precision*resXClase.recall/(resXClase.precision+resXClase.recall);
					}
					else {
						resXClase.f1 = 0;
					}
					resXClase.clase = j;
					resultadosTemp.push_back(resXClase);
				}

//**********************Codigo para la matriz de confusion **************//
				uint vector_size = vectordeKnns.size();
				vector< vector<uint> > matrizConfusion (cantidadDeClases, vector<uint> (cantidadDeClases,0));
				for (uint elem = 0; elem < vector_size; ++elem){
						++matrizConfusion[labelsYTemp[elem]-1][vectordeKnns[elem]-1]; //el -1 es porque las clases van de 1 a 10 y aca la matriz se indexa de 0 a 9
//la idea es, indexo a la fila segun la clase del elemento, y luego indexo a la columna segun lo que devolvio el knn
//una columna sumada contiene la cantidad de veces que el knn eligio esa clase como resultado
//idealmente la diagonal es la que va a tener numeros mas grandes
//al mirar una fila i vemos que devolvio el knn para la clase i+1
//al mirar una columna j vemos a que clase pertenecia lo que el knn determino que era de la clase j+1
				}
//**********************Codigo para la matriz de confusion END**************//
				double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
				resVariandoKSinPCAParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp)); //los resultados entran dependiendo del k, los de k=1 van primeros y los de k =321 van ultimos
			} //aca terminaria el for que varia el kDeKnn

			escribirEstadisiticas("./Resultados/ResultadosVariandoKSinPCA", resVariandoKSinPCAParaUnFold,i,0,false); //false es que estoy variando elkdeKnn, pongo 0 porque no importa en este caso ya que no estoy variando el alpha
		//res.push_back(make_pair(resultadosTemp,accuracyTemp)); //comente esto porque no importa mucho lo que devuelve, solo queremos escribir los resultados en archivos, luego de experimentar hay que cambiar esto
		}
	}
	return res;	
}






int main(int argc, char * argv[]) {
    string metodo, trainSet, testSet, classif;
/*    salida.open("Comparación_de_algoritmos.txt", ios_base::app);
    salida << "Rendimiento relativo, algor.1/algor.2" << endl;
    salida.close();*/

/*    if (!obtenerParametros(argc, argv, &metodo, &trainSet, &testSet, &classif)) {
        cout << "Modo de uso: tp2 -m <method> -i <train_set> -q <test_set> -o <classif>\n";
    } else {*/
		/*vector<vector<double>>* dataSet = new vector<vector<double> >;
		vector<uint>* labelsX = new vector<uint >;*/


		//------- cargamos los datos de uno de los tests en la funcion cargarTest esta la explicacion de que hace-------------------//
        vector<vector<double>>* dataSetTest = new vector<vector<double> >;
        vector<uint>* labelsTest = new vector<uint>;
        vector<double>* autovaloresTest = new vector<double>;
        cargarTest("./tests/testFullBig", dataSetTest, labelsTest, autovaloresTest);
        //------- cargamos los datos de uno de los tests en la funcion cargarTest esta la explicacion de que hace-------------------//

		//cargarDataSetEnMatriz("./ImagenesCarasRed",dataSet, labelsX);
		vector<pair<vector<resultados >,double> > dasdsa = kFold(*dataSetTest,*labelsTest,5,328,328,false); //328 es la cantidad de imagenes de trainX por cada fold en el testBigFull que a su vez coincide con el rango de la matriz de covarianza y por ende la cantidad maxima de autovalores y autovectores que podemos conseguir con el metodoPotencia antes de que empiece a colgarse porque los autovalores empiezan a ser 0
       		//escribirEstadisiticas("./pruebaEstadisticas", dasdsa);
		/*delete labelsX;
		delete dataSet;*/
		delete dataSetTest;
		delete labelsTest;
		delete autovaloresTest;

    /*}*/

    return 0;
}

