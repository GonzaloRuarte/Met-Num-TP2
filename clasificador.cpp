//
// Created by christian on 25/05/18.
//

#include "clasificador.h"


using namespace std;

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
        cout << "producto_interno: los vectores no son del mismo tamaño" << endl;
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

pair<double,vector<double> > metodoPotencia(const vector<vector<double> > &M) {
    const size_t& n = M[0].size();
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
    while(diferencia2 >= 0.0001 or cantidad_iteraciones2 < 30){
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
            normalizar2(autovector2);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
            ++cantidad_iteraciones2;
        }
    }else{
        diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
        while(diferencia2 > 0.001){
            autovector2_temp = mult_matr_por_vect(M, autovector2);
            normalizar2(autovector2_temp);
            autovector2 = mult_matr_por_vect(M, autovector2_temp);
            normalizar2(autovector2);
            diferencia2 = norma1(restaVec(autovector2, autovector2_temp));
            ++cantidad_iteraciones2;
        }
    }
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

vector<pair<double,clase_t> > vector_de_distancias(const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<double>& newImg){
    vector<pair<double,clase_t> > vecNormas (trainX.size());
    for (uint i = 0; i < labelsX.size(); i++) {
        vecNormas[i].first  = norma2(restaVec(trainX[i],newImg));
        vecNormas[i].second = labelsX[i];
    }
    return vecNormas;
}

clase_t Knn_sim (const vector<pair<double,clase_t> >& vecNormas, uint k){
    vector<pair<double,clase_t> >::const_iterator primero = vecNormas.cbegin();
    vector<pair<double,clase_t> >::const_iterator k_esimo = primero+k;
    priority_queue<pair<double,clase_t> > heap(primero, k_esimo);  //Creo un max_heap con los primeros k elementos de vecNormas
    for(size_t i = k; i < vecNormas.size(); ++i){
        if(vecNormas[i] < heap.top()){  //Si el i-ésimo elemento es más chico que el más grande del heap...
            heap.pop();                 //entonces saco al más grande...
            heap.push(vecNormas[i]);    //y meto al i-ésimo elemento.
        }                               //De esta forma me quedo con los k elementos más chicos.
    }

    map<clase_t, pair<uint, double> > candidatos;   //La 1° coordenada del significado cuenta la cantidad de apariciones.
    while(!heap.empty()){                           //La 2° coordenada es la distancia más grande de todas las distancias que hay para una misma clase.
        pair<clase_t, pair<uint, double> > clave_signif = make_pair(heap.top().second, make_pair(1, heap.top().first));
        pair<map<clase_t, pair<uint, double> >::iterator, bool> it_bool = candidatos.insert(clave_signif);
        if(not it_bool.second)  //Si la clase ya estaba en el map, debo únicamente sumar 1 a la cant. de apariciones.
            ++it_bool.first->second.first;
        heap.pop();
    }   //Cada clase queda asociada a la 1° distancia con la que entró, que es la más grande.
    auto it = candidatos.begin();
    clase_t clase = it->first;
    uint cant_apari = it->second.first;
    double max_dist = it->second.second;
    while (it != candidatos.end()){
        if(cant_apari < it->second.first || (cant_apari == it->second.first && max_dist > it->second.second)){
            clase = it->first;
            cant_apari = it->second.first;
            max_dist = it->second.second;
        }
        ++it;
    }
    return clase;
}

clase_t Knn (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<double>& newImg, uint k) {//las labels podrian no ser un uint pero lo dejo asi en un principio
    vector<pair<double,clase_t> > vecNormas (trainX.size());
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

    map<clase_t, pair<uint, double> > candidatos;   //La 1° coordenada del significado cuenta la cantidad de apariciones.
    while(!heap.empty()){                           //La 2° coordenada es la distancia más grande de todas las distancias que hay para una misma clase.
        pair<clase_t, pair<uint, double> > clave_signif = make_pair(heap.top().second, make_pair(1, heap.top().first));
        pair<map<clase_t, pair<uint, double> >::iterator, bool> it_bool = candidatos.insert(clave_signif);
        if(not it_bool.second)  //Si la clase ya estaba en el map, debo únicamente sumar 1 a la cant. de apariciones.
            ++it_bool.first->second.first;
        heap.pop();
    }   //Cada clase queda asociada a la 1° distancia con la que entró, que es la más grande.
    auto it = candidatos.begin();
    clase_t clase = it->first;
    uint cant_apari = it->second.first;
    double max_dist = it->second.second;
    while (it != candidatos.end()){
        if(cant_apari < it->second.first || (cant_apari == it->second.first && max_dist > it->second.second)){
            clase = it->first;
            cant_apari = it->second.first;
            max_dist = it->second.second;
        }
        ++it;
    }
    return clase;
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
    return V; //devuelvo la V, recordar multiplicar fuera de la funcion
}

vector<uint> vectorDeKnns (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, const vector<vector<double> >& testY, uint k){
    const unsigned long& n = testY.size();
    vector<uint> res (n);
    for (uint i = 0; i < n; i++) {
        res[i] = Knn(trainX,labelsX,testY[i],k);
    }
    return res;
}

