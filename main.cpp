#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <tuple>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "rdtsc.h"
#include "util.h"
#include <cmath>
#include <chrono>

using namespace std;

typedef struct {
	uint clase;
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
	    for (uint k = 0; k < n; k++)
		    var_i += imgs[k][i]*imgs[k][i]; //calculo de la sumatoria de productos para calcular varianza
        res[i][i] = var_i/(n-1);
		for (uint j = i+1; j < m; j++){ //como la matriz es simetrica basta calcular la mitad superior.
            covar_ij = 0;
			for (uint k = 0; k < imgs.size(); k++)
			    covar_ij += imgs[k][i]*imgs[k][j]; //calculo de la sumatoria de productos para calcular covarianza
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
    auto t3 = chrono::system_clock::now();
    while(diferencia2 >= 0.000001){
        autovector2_temp = mult_matr_por_vect(M, autovector2);
        autovalor2_temp = producto_interno(autovector2, autovector2_temp); //autovector está normalizado.
        normalizar2(autovector2_temp);
        autovector2 = mult_matr_por_vect(M, autovector2_temp);
        autovalor2 = producto_interno(autovector2_temp, autovector2); //autovector_temp está normalizado.
        normalizar2(autovector2);
        diferencia2 = abs(autovalor2-autovalor2_temp);
        ++cantidad_iteraciones2;
    }
    auto t4 = chrono::system_clock::now();
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

vector<vector<double> > generarP(const vector<vector<double> > &mat, uint alpha){
	vector<vector<double> > res (alpha);
	vector< pair<double,vector<double> > > autovectores = deflacion(mat,alpha);
	for (uint i = 0; i< alpha; i++){
		res[i] = autovectores[i].second;
	}
	return res;
}


uint Knn (const vector<vector<double> >& trainX, const vector<uint>& labelsX, const vector<double>& newImg, uint k) {//las labels podrian no ser un uint pero lo dejo asi en un principio
	vector<pair<double,uint> > vecNormas (trainX.size());
	vector<pair<double,uint> > sorted (k);
	double temp;
	for (uint i = 0; i < labelsX.size(); i++) {
		temp = norma2(restaVec(trainX[i],newImg));
		vecNormas[i].first = temp;
		vecNormas[i].second = labelsX[i];
	}
    vector<pair<double,uint> >::const_iterator primero = vecNormas.cbegin();
    vector<pair<double,uint> >::const_iterator k_esimo = primero+k;
	priority_queue<pair<double,uint> > heap(primero, k_esimo);  //Creo un max_heap con los primeros k elementos de vecNormas
	for(size_t i = k; i < vecNormas.size(); ++i){
        if(vecNormas[i] < heap.top()){  //Si el i-ésimo elemento es más chico que el más grande del heap...
            heap.pop();                 //entonces saco al más grande...
            heap.push(vecNormas[i]);    //y meto al i-ésimo elemento.
        }                               //De esta forma me quedo con los k elementos más chicos.
	}
	for(int i = k-1; i >= 0; --i){
	    sorted[i] = heap.top();
	    heap.pop();
	}
/*	for(uint i = 0; i < k; i++) {//sort de menor a mayor segun las normas
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
	}*/
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

double accuracy (const vector<vector<double> >& trainX, const vector<uint>& labelsX, const vector<vector<double> >& testY, const vector<uint>& labelsY, uint k) {
	const unsigned long& n = labelsY.size(); //me di cuenta que si vamos a usar como label el numero de fila dentro de la matriz entonces el labelsY de mucho no nos sirve pero bueno
	double acum = 0;
	for (uint i = 0; i < n; i++) {
		if (Knn(trainX,labelsX,testY[i],k) == labelsY[i]){
			acum++;
		}
	}
	return acum/n;
}	

double precision(const vector<vector<double> >& trainX, const vector<uint>& labelsX, const vector<vector<double> >& testY, const vector<uint>& labelsY, uint clase, uint k) {
    const unsigned long& n = labelsY.size();
	double truepositives = 0, positives = 0;
	for (uint i = 0; i < n; i++) {
		uint Knnres = Knn(trainX,labelsX,testY[i],k);
		if (Knnres == clase){ //si el Knn dio igual a la clase que estoy procesando entonces sumo 1 a los elementos recuperados
			positives++;
			if (Knnres == labelsY[i]){ //si el Knn ademas dio bien el resultado sumo 1 a los true positives
				truepositives++;
			}
		}
	}
	return truepositives/positives;
}

double recall(const vector<vector<double> >& trainX, const vector<uint>& labelsX, const vector<vector<double> >& testY, const vector<uint>& labelsY, uint clase, uint k) {
    const unsigned long& n = labelsY.size();
	double truepositives = 0, relevants = 0;
	double res;
	for (uint i = 0; i < n; i++) {
		if (labelsY[i] == clase){//si el elemento que estoy analizando pertenece a la clase que estoy procesando, sumo 1 a los relevantes
		 	relevants++;
			if (Knn(trainX,labelsX,testY[i],k) == clase){ //si el Knn ademas dio bien el resultado sumo 1 a los true positives
				truepositives++;
			}
		}
	}
	res = truepositives/relevants;
	return res;
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
	uint m = trainX[0].size();
	vector<vector<double> > Mx = calcularMx(trainX);
	vector<vector<double> > V = trasponer(generarP(Mx,alpha));
    convertirMatrizAImagen("./salidaVtraspuesta", 10, &V);
	/*for (uint i = 0; i < m; i++){ //esto no hace falta por ahora porque la V se calcula ya con alpha columnas
		V[i].erase(V[i].begin()+alpha, V[i].end());
	}*/
	return V; //devuelvo la V, recordar multiplicar fuera de la funcion
}
vector<pair<vector<resultados >,double> > kFold (vector<vector<double> > trainX, vector<uint> labelsX, uint k, uint kdeKnn, uint alpha) {
	uint imagenesPorPersona = 10; //esto podria variar si cambiamos el trainX
	int imagenesPPparagenerador = 10; //es el mismo numero de arriba pero lo uso para generar numeros aleatorios
	uint cantidadDeClases = 41; // idem arriba
	vector<int> folds; //to store the random numbers
	random_device rd; //seed generator
	mt19937_64 generator{rd()}; //generator initialized with seed from rd
	uniform_int_distribution<> dist{0, imagenesPPparagenerador-1}; //the range is inclusive, so this produces numbers in range [0, 10)
	for(uint i=0; i<imagenesPorPersona; ++i) {
		folds.push_back( dist(generator) );
	} // la idea es que voy a tener muestras balanceadas, entonces para cada persona voy a tener la misma cantidad de imagenes en test y en train
		// como cada persona tiene 10 imagenes, el k puede ser 1, 2, 5 o 10, k = 1 no tiene mucho sentido
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
		vector<vector<double>> V = PCA(trainXTemp,6);
		trainXTemp = multMat(trainXTemp,V);
		testYTemp = multMat(testYTemp,V);
		vector<resultados > resultadosTemp;
		for (uint j = 1; j < cantidadDeClases; j++){ //itero sobre las clases
			resultados resXClase;
			resXClase.precision = precision(trainXTemp,labelsXTemp,testYTemp,labelsYTemp,j,kdeKnn);
			resXClase.recall = recall(trainXTemp,labelsXTemp,testYTemp,labelsYTemp,j,kdeKnn);
			resXClase.f1 = 2.0*resXClase.precision*resXClase.recall/(resXClase.precision+resXClase.recall);
			resXClase.clase = j;
			resultadosTemp.push_back(resXClase);
		}//entonces en el vector la posicion 0 corresponde a la clase 1 y asi sucesivamente
		double accuracyTemp = accuracy(trainXTemp,labelsXTemp,testYTemp,labelsYTemp,kdeKnn);
		res.push_back(make_pair(resultadosTemp,accuracyTemp));
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
		vector<vector<double>>* dataSet = new vector<vector<double> >;
		vector<uint>* labelsX = new vector<uint > (41);


		cargarDataSetEnMatriz("./ImagenesCarasRed",dataSet, labelsX);

		vector<pair<vector<resultados >,double> > dasdsa = kFold(*dataSet,*labelsX,5,2,6);


		delete labelsX;
		delete dataSet;

    /*}*/

    return 0;
}
