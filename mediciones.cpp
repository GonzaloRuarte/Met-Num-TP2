//
// Created by christian on 25/05/18.
//

#include "mediciones.h"

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

//------------------------- Escritura de las estadisticas -------------------------//
void escribirTiempos(string nombreArchivo, vector<vector<unsigned long> > &tiempos, bool conPCA, bool varioAlpha, int variacion, uint kdeKnninit, uint alpha) { //si varioAlpha es false es porque estoy variando el kdeKnn
    //vector<unsigned long>* tiempo;
    int i;
    if(varioAlpha){
        i=alpha;
    }else{
        i= kdeKnninit;
    }
    for (vector<vector<unsigned long> >::iterator it = tiempos.begin() ; it != tiempos.end(); ++it) {
        vector<unsigned long>& tiempo = *it;
        ofstream salida;
        string valor_parametro = int2stringConCantidadDigitos(4, i);
        if (conPCA){
            if(varioAlpha){
                salida = getFlujo(nombreArchivo + "K_" + to_string(kdeKnninit) + "Alpha_" + valor_parametro);
            } else{
                salida = getFlujo(nombreArchivo + "Alpha_" + to_string(alpha) + "K_" + valor_parametro);
            }
        } else{
            salida = getFlujo(nombreArchivo + "K_" + valor_parametro);
        }
        for (vector<unsigned long>::iterator it2 = tiempo.begin() ; it2 != tiempo.end(); ++it2) {
            salida << *it2 << endl;
        }
        salida.close();
        i+=variacion; //le sumo lo que fui variando
    }


}

void escribirEstadisiticas(string nombreArchivo, vector<pair<vector<resultados >,double> > &estadisticas, uint kDekfold,uint alpha, bool varioAlpha, int variacion, uint kdeKnninit, bool conPCA) { //si varioAlpha es false es porque estoy variando el kdeKnn
    vector<resultados>* estadistica;
    int i;
    if(varioAlpha){
        i=alpha;
    }else{
        i= kdeKnninit;
    }
    for (vector<pair<vector<resultados >,double> >::iterator it = estadisticas.begin() ; it != estadisticas.end(); ++it) {
        vector<resultados >& estadistica = it->first;
        string accuracy = to_string(it->second);
        ofstream salida;
        string valor_parametro = int2stringConCantidadDigitos(4, i);
        if (conPCA){
            if(varioAlpha){
                salida = getFlujo(nombreArchivo + "K_" + to_string(kdeKnninit) + "Alpha_" + valor_parametro);
            } else{
                salida = getFlujo(nombreArchivo + "Alpha_" + to_string(alpha) + "K_" + valor_parametro);
            }
        } else{
            salida = getFlujo(nombreArchivo + "K_" + valor_parametro);
        }
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
            i-=variacion; //le resto lo que fui variando el alpha
        }else{
            i+=variacion; //le sumo lo que fui variando el kDeKnn
        }
    }
}
//------------------------- Escritura de las estadisticas -------------------------//

vector<pair<vector<resultados >,double> > kFold (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k, uint kdeKnn, uint alpha, bool conPCA, bool varioAlpha) {
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
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona); //to store the random numbers
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector< vector<int>> foldsXpersona (n);
    mt19937 g(static_cast<uint32_t>(time(0)));
    for (uint i = 0; i<n;++i){
        shuffle(folds.begin(), folds.end(),g);
        foldsXpersona[i] = folds;
    }
	// la idea es que voy a tener muestras balanceadas, entonces para cada persona voy a tener la misma cantidad de imagenes en test y en train
		// como cada persona tiene 10 imagenes, el k puede ser 1, 2, 5 o 10, k = 1 no tiene mucho sentido*/
    vector<pair<vector<resultados >,double> > res;
    for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
        vector<vector<double> > trainXTemp;
        vector<vector<double> > testYTemp;
        vector<uint> labelsXTemp;
        vector<uint> labelsYTemp;
        for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
            for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
                uint temp = (j*imagenesPorPersona)+foldsXpersona[j][u];
                if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
                    testYTemp.push_back(trainX[temp]);
                    labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
                } else{
                    trainXTemp.push_back(trainX[temp]);
                    labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
                }

            }
        }

        if (conPCA){
            if(varioAlpha){
                vector<vector<double>> V = PCATecho(trainXTemp,alpha);
                uint size_V = V[0].size();
                vector<pair<vector<resultados >,double> > resVariandoAlphaParaUnFold;
                for(uint h = 0; h <= alpha-1; h+=20) {//esto sirve para iterar el alpha (voy borrando columnas de la matriz V dependiendo del h)
                    for (uint o = 0; o < V.size(); o++){ //borro las columnas de V que necesito borrar para variar el alpha
                        V[o].erase(V[o].begin()+size_V-h, V[o].end());
                    }
                    vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
                    vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
                    //for(uint y = 0; y < 321; ++y) { //este for seria para variar el kDeKnn
                    vector<resultados > resultadosTemp (0);
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
                escribirEstadisiticas("./Resultados/ResultadosVariandoAlpha"+to_string(kdeKnn)+"/ResultadosVariandoAlpha", resVariandoAlphaParaUnFold,i,alpha,true,20,kdeKnn,true); //primer true es que estoy variando el alpha, 20 es la variacion del alpha
//el segundo bool es que estoy usando PCA

                V = PCATecho(trainXTemp,61);
                size_V = V[0].size();
                vector<pair<vector<resultados >,double> > resVariandoAlphaParaUnFoldFina;
                for(uint h = 0; h < 61; ++h) {//puse 41 para que alpha varie desde 41 a 1
                    for (uint o = 0; o < V.size(); o++){ //borro las columnas de V que necesito borrar para variar el alpha
                        V[o].erase(V[o].begin()+size_V-h, V[o].end());
                    }
                    vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
                    vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
                    //for(uint y = 0; y < 321; ++y) { //este for seria para variar el kDeKnn
                    vector<resultados > resultadosTemp (0);
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
                    }
                    double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);

                    resVariandoAlphaParaUnFoldFina.push_back(make_pair(resultadosTemp,accuracyTemp));
                    //} //aca terminaria el for que varia el kDeKnn
                }
                escribirEstadisiticas("./Resultados/ResultadosVariandoAlphaFina"+to_string(kdeKnn)+"/ResultadosVariandoAlphaFina", resVariandoAlphaParaUnFoldFina,i,61,true,1,kdeKnn,true); //41 porque el alpha varia de 41 a 1
            }else{//vario k
                vector<vector<double>> V = PCATecho(trainXTemp,alpha);
                uint size_V = V[0].size();
                vector<pair<vector<resultados >,double> > resVariandoKParaUnFold;
                vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
                vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
                for(uint y = 1; y < kdeKnn; y+=20) { //este for seria para variar el kDeKnn
                    vector<resultados > resultadosTemp (0);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,y);
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
                    resVariandoKParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp));
                } //aca terminaria el for que varia el kDeKnn
                escribirEstadisiticas("./Resultados/ResultadosVariandoKConPCA"+to_string(alpha)+"/ResultadosVariandoKConPCA", resVariandoKParaUnFold,i,alpha,false,20,1,true); //false es que estoy variando el k, 20 es la variacion del alpha
//el segundo bool es que estoy usando PCA
                vector<pair<vector<resultados >,double> > resVariandoKParaUnFoldFina;
                for(uint y = 1; y <= 41; y++) { //este for seria para variar el kDeKnn
                    vector<resultados > resultadosTemp (0);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,y);
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
                    resVariandoKParaUnFoldFina.push_back(make_pair(resultadosTemp,accuracyTemp));
                } //aca terminaria el for que varia el kDeKnn
                escribirEstadisiticas("./Resultados/ResultadosVariandoKConPCAFina"+to_string(alpha)+"/ResultadosVariandoKConPCAFina", resVariandoKParaUnFoldFina,i,alpha,false,1,1,true); //false es que estoy variando el k, 20 es la variacion del alpha
//el segundo bool es que estoy usando PCA


            }
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
                double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
                resVariandoKSinPCAParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp)); //los resultados entran dependiendo del k, los de k=1 van primeros y los de k =321 van ultimos
            } //aca terminaria el for que varia el kDeKnn

            escribirEstadisiticas("./Resultados/ResultadosVariandoKSinPCA/ResultadosVariandoKSinPCA", resVariandoKSinPCAParaUnFold,i,0,false,20,1,false); //false es que estoy variando elkdeKnn, pongo 0 porque no importa en este caso ya que no estoy variando el alpha, 20 es porque vario el k de a 20, el 1 es porque el k arranca en 1 en este caso

            vector<pair<vector<resultados >,double> > resVariandoKSinPCAParaUnFoldFina;
            for(uint y = 1; y <= 41; ++y) { //este for seria para variar el kDeKnn
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

                double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
                resVariandoKSinPCAParaUnFoldFina.push_back(make_pair(resultadosTemp,accuracyTemp));
            } //aca terminaria el for que varia el kDeKnn

            escribirEstadisiticas("./Resultados/ResultadosVariandoKSinPCAFina/ResultadosVariandoKSinPCAFina", resVariandoKSinPCAParaUnFoldFina,i,0,false,1,1,false);//el primer 1 es porque vario el k de a 1 y el segundo 1 porque esta iteracion arranca con k = 1
        }
    }
    return res;
}

void kFoldDark (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k) {
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
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona); //to store the random numbers
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector< vector<int>> foldsXpersona (n);
    mt19937 g(static_cast<uint32_t>(time(0)));
    for (uint i = 0; i<n;++i){
        shuffle(folds.begin(), folds.end(),g);
        foldsXpersona[i] = folds;
    }
	// la idea es que voy a tener muestras balanceadas, entonces para cada persona voy a tener la misma cantidad de imagenes en test y en train
		// como cada persona tiene 10 imagenes, el k puede ser 1, 2, 5 o 10, k = 1 no tiene mucho sentido*/
    vector<pair<vector<resultados >,double> > res;
    for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
        vector<vector<double> > trainXTemp;
        vector<vector<double> > testYTemp;
        vector<uint> labelsXTemp;
        vector<uint> labelsYTemp;
        for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
            for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
                uint temp = (j*imagenesPorPersona)+foldsXpersona[j][u];
                if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
                    testYTemp.push_back(trainX[temp]);
                    labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
                } else{
                    trainXTemp.push_back(trainX[temp]);
                    labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
                }

            }
        }

        vector<vector<double>> V = PCATecho(trainXTemp,31);
        vector<pair<vector<resultados >,double> > resVariandoDarkParaUnFold;
        vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
        for(uint y = 1; y <= 10; ++y) { //este for seria para variar lo que oscurece

            vector<vector<double> > testYTemp2 = testYTemp;
            multMatEsc(testYTemp2,y*0.10);
            testYTemp2 = multMat(testYTemp2,V);
            vector<resultados > resultadosTemp (0);
            vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,1);
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
            resVariandoDarkParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp));
        } //aca terminaria el for que varia el kDeKnn



        escribirEstadisiticas("./Resultados/ResultadosDark/ResultadosDark", resVariandoDarkParaUnFold,i,1,false,20,1,true);
//si usas esta funcion el nombre de los archivos va a ser cualquier cosa pero bue jaja
        //de ultima hayque hacer otra funcion para que escriba bien el nombre de los archivos
    }
}

void kFoldDarkSat (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k) {
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
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona); //to store the random numbers
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector< vector<int>> foldsXpersona (n);
    mt19937 g(static_cast<uint32_t>(time(0)));
    for (uint i = 0; i<n;++i){
        shuffle(folds.begin(), folds.end(),g);
        foldsXpersona[i] = folds;
    }
	// la idea es que voy a tener muestras balanceadas, entonces para cada persona voy a tener la misma cantidad de imagenes en test y en train
		// como cada persona tiene 10 imagenes, el k puede ser 1, 2, 5 o 10, k = 1 no tiene mucho sentido*/
    vector<pair<vector<resultados >,double> > res;
    for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
        vector<vector<double> > trainXTemp;
        vector<vector<double> > testYTemp;
        vector<uint> labelsXTemp;
        vector<uint> labelsYTemp;
        for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
            for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
                uint temp = (j*imagenesPorPersona)+foldsXpersona[j][u];
                if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
                    testYTemp.push_back(trainX[temp]);
                    labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
                } else{
                    trainXTemp.push_back(trainX[temp]);
                    labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
                }

            }
        }

        vector<vector<double>> V = PCATecho(trainXTemp,31);
        vector<pair<vector<resultados >,double> > resVariandoDarkParaUnFold;
        vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
        for(uint y = 0; y <= 9; ++y) { //este for seria para variar lo que oscurece, primera iteracion lo dejo como esta y la ultima le resto 180 a todos

            vector<vector<double> > testYTemp2 = testYTemp;
            for(uint o = 0; o < testYTemp2.size(); ++o){
                for(uint h = 0; h < testYTemp2[0].size(); ++h){
                    testYTemp2[o][h]-=20*y;
                    if (testYTemp2[o][h]< 0){
                        testYTemp2[o][h]=0;
                    }
                }
            }
            testYTemp2 = multMat(testYTemp2,V);
            vector<resultados > resultadosTemp (0);
            vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,1);
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
            resVariandoDarkParaUnFold.push_back(make_pair(resultadosTemp,accuracyTemp));
        } //aca terminaria el for que varia el kDeKnn



        escribirEstadisiticas("./Resultados/ResultadosDarkSat/ResultadosDarkSat", resVariandoDarkParaUnFold,i,10,true,20,1,true);
//si usas esta funcion el nombre de los archivos va a ser cualquier cosa pero bue jaja
        //de ultima hayque hacer otra funcion para que escriba bien el nombre de los archivos
    }
}

void medirTiempos (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k, uint kdeKnn, uint alpha, bool conPCA, bool varioAlpha) {
    //codigo para calcular la cantidad de imagenes por persona suponiendo que las muestras son balanceadas y la cantidad de clases
    uint imagenesPorPersona = 0;
    uint cantidadDeClases = 0;
    uint personasTemp = 0;
    for( uint i = 0; i < labelsX.size(); i++){
        if (labelsX[i] != personasTemp){
            personasTemp = labelsX[i];
            cantidadDeClases++;
        }
        if(personasTemp == 1){
            imagenesPorPersona++;
        }
    }
//***********************************************************************//
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona);
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector<vector<double> > trainXTemp;
    vector<vector<double> > testYTemp;
    vector<uint> labelsXTemp;
    vector<uint> labelsYTemp;
    for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
        for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
            uint temp = (j*imagenesPorPersona)+u;
            if (u < imagenesPorPersona/k){ //si estoy en el fold que quiero
                testYTemp.push_back(trainX[temp]);
                labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
            } else{
                trainXTemp.push_back(trainX[temp]);
                labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
            }

        }
    }
    if(conPCA){
        if(varioAlpha){
            vector<vector<unsigned long> > vectorTiemposYAlpha (17);
            vector<vector<double>> V = PCATecho(trainXTemp,alpha);
            uint size_V = V[0].size();
            for(uint h = 0; h <= alpha-1; h+=20) {//esto sirve para iterar el alpha (voy borrando columnas de la matriz V dependiendo del h)
                for (uint o = 0; o < V.size(); o++){ //borro las columnas de V que necesito borrar para variar el alpha
                    V[o].erase(V[o].begin()+size_V-h, V[o].end());
                }
                vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
                vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
                for (int i = 0; i < 20; i++) {
                    unsigned long start, end;
                    unsigned long delta = 0;
                    RDTSC_START(start);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,kdeKnn);
                    RDTSC_STOP(end);
                    delta = end - start;//cada delta es el tiempo que tarda en calcular el PCA+ aplicar el cambio de base + calcular el knn para todos los elementos de testY
                    vectorTiemposYAlpha[(size_V-h-1)/20].push_back(delta);
                }
            }
            escribirTiempos("./Resultados/TiemposVariandoAlpha"+to_string(kdeKnn)+"/TiemposVariandoAlpha", vectorTiemposYAlpha,true,true,20,kdeKnn,1);

            vector<vector<unsigned long> > vectorTiemposYAlphaFina (61);
            V = PCATecho(trainXTemp,61);
            size_V = V[0].size();
            for(uint h = 0; h <= 60; ++h) {//esto sirve para iterar el alpha (voy borrando columnas de la matriz V dependiendo del h)
                for (uint o = 0; o < V.size(); o++){ //borro las columnas de V que necesito borrar para variar el alpha
                    V[o].erase(V[o].begin()+size_V-h, V[o].end());
                }
                vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
                vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
                for (int i = 0; i < 20; i++) {
                    unsigned long start, end;
                    unsigned long delta = 0;
                    RDTSC_START(start);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,kdeKnn);
                    RDTSC_STOP(end);
                    delta = end - start;//cada delta es el tiempo que tarda en calcular el PCA+ aplicar el cambio de base + calcular el knn para todos los elementos de testY
                    vectorTiemposYAlphaFina[size_V-h-1].push_back(delta);
                }
            }
            escribirTiempos("./Resultados/TiemposVariandoAlphaFina"+to_string(kdeKnn)+"/TiemposVariandoAlphaFina", vectorTiemposYAlphaFina,true,true,1,kdeKnn,1);

        }else{
            vector<vector<unsigned long> > vectorTiemposYK (17);
            vector<vector<double>> V = PCATecho(trainXTemp,alpha);
            vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
            vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
            for(uint y = 1; y <= kdeKnn; y+=20){
                for (int i = 0; i < 20; i++) {
                    unsigned long start, end;
                    unsigned long delta = 0;
                    RDTSC_START(start);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,y);
                    RDTSC_STOP(end);
                    delta = end - start;//cada delta es el tiempo que tarda en calcular el PCA+ aplicar el cambio de base + calcular el knn para todos los elementos de testY
                    vectorTiemposYK[(y-1)/20].push_back(delta);
                }
            }
            escribirTiempos("./Resultados/TiemposVariandoKConPCA"+to_string(alpha)+"/TiemposVariandoKConPCA", vectorTiemposYK,true,false,20,1,alpha);


            vector<vector<unsigned long> > vectorTiemposYKFina (41);
            V = PCATecho(trainXTemp,alpha);
            trainXTemp2 = multMat(trainXTemp,V);
            testYTemp2 = multMat(testYTemp,V);
            for(uint y = 1; y <= 41; ++y){
                for (int i = 0; i < 20; i++) {
                    unsigned long start, end;
                    unsigned long delta = 0;
                    RDTSC_START(start);
                    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,y);
                    RDTSC_STOP(end);
                    delta = end - start;//cada delta es el tiempo que tarda en calcular el PCA+ aplicar el cambio de base + calcular el knn para todos los elementos de testY
                    vectorTiemposYKFina[y-1].push_back(delta);
                }
            }
            escribirTiempos("./Resultados/TiemposVariandoKConPCAFina"+to_string(alpha)+"/TiemposVariandoKConPCAFina", vectorTiemposYKFina,true,false,1,1,alpha);
        }

    }else{
        vector<vector<unsigned long> > vectorTiemposYK (17);
        for(uint y = 1; y <= kdeKnn; y+=20) { //este for seria para variar el kDeKnn
            for (int i = 0; i < 20; i++) {
                unsigned long start, end;
                unsigned long delta = 0;
                RDTSC_START(start);
                vector<uint> vectordeKnns = vectorDeKnns(trainXTemp,labelsXTemp,testYTemp,y);
                RDTSC_STOP(end);
                delta = end - start;//cada delta es el tiempo que tarda en calcular el knn para todos los elementos de testY
                vectorTiemposYK[(y-1)/20].push_back(delta);
            }

        }
        escribirTiempos("./Resultados/TiemposVariandoKSinPCA/TiemposVariandoKSinPCA", vectorTiemposYK,false,false,20,1,0); //el primer 20 es la variacion, el 1 es el k inicial, el 0 es el alpha que aca no importa

        vector<vector<unsigned long> > vectorTiemposYKFina (41);
        for(uint y = 1; y <= 41; ++y) { //este for seria para variar el kDeKnn
            for (int i = 0; i < 20; i++) {
                unsigned long start, end;
                unsigned long delta = 0;
                RDTSC_START(start);
                vector<uint> vectordeKnns = vectorDeKnns(trainXTemp,labelsXTemp,testYTemp,y);
                RDTSC_STOP(end);
                delta = end - start;//cada delta es el tiempo que tarda en calcular el knn para todos los elementos de testY
                vectorTiemposYKFina[y-1].push_back(delta);
            }

        }
        escribirTiempos("./Resultados/TiemposVariandoKSinPCAFina/TiemposVariandoKSinPCAFina", vectorTiemposYKFina,false,false,1,1,0);
    }
}

void escribirTiemposDataset(string nombreArchivo, vector<vector<unsigned long> > &tiempos, int variacion) { //si varioAlpha es false es porque estoy variando el kdeKnn
    //vector<unsigned long>* tiempo;
    int i = 41;

    for (vector<vector<unsigned long> >::iterator it = tiempos.begin() ; it != tiempos.end(); ++it) {
        vector<unsigned long>& tiempo = *it;
        ofstream salida;
        string valor_parametro = int2stringConCantidadDigitos(4, i);
        salida = getFlujo(nombreArchivo + "TamDataset_" + valor_parametro);
        for (vector<unsigned long>::iterator it2 = tiempo.begin() ; it2 != tiempo.end(); ++it2) {
            salida << *it2 << endl;
        }
        salida.close();
        i+=variacion; //le sumo lo que fui variando
    }
}

void escribirDistorted (string nombreArchivo, pair<vector<resultados >,double> &estadisticas) {

    string accuracy = to_string(estadisticas.second);
    ofstream salida;
    salida = getFlujo(nombreArchivo + "Distorted");
    string precision = "";
    string recall = "";
    string f1="";
    for (vector<resultados>::iterator it = estadisticas.first.begin() ; it != estadisticas.first.end(); ++it) {
        precision += to_string(it->precision) + "\t";
        recall += to_string(it->recall) + "\t";
        f1 += to_string(it->f1) + "\t";
    }
    salida << accuracy << endl;
    salida << precision << endl;
    salida << recall << endl;
    salida << f1 << endl;
    salida.close();
}

void testDistorted (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k) {
    uint imagenesPorPersona = 0;
    uint cantidadDeClases = 0;
    uint personasTemp = 0;
    for( uint i = 0; i < labelsX.size(); i++){
        if (labelsX[i] != personasTemp){
            personasTemp = labelsX[i];
            cantidadDeClases++;
        }
        if(personasTemp == 1){
            imagenesPorPersona++;
        }
    }
//***********************************************************************//
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona);
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector<vector<double> > trainXTemp;
    vector<vector<double> > testYTemp;
    vector<uint> labelsXTemp;
    vector<uint> labelsYTemp;
    for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
        for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
            uint temp = (j*imagenesPorPersona)+u;
            if (u < imagenesPorPersona/k){ //si estoy en el fold que quiero
                testYTemp.push_back(trainX[temp]);
                labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
            } else{
                trainXTemp.push_back(trainX[temp]);
                labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
            }

        }
    }

    vector<vector<double>> V = PCATecho(trainXTemp,31);
    vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
    vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
    vector<resultados > resultadosTemp (0);
    vector<uint> vectordeKnns = vectorDeKnns(trainXTemp,labelsXTemp,testYTemp,1);
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

    double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
    pair<vector<resultados >,double> resDistorted = make_pair(resultadosTemp,accuracyTemp);


    escribirDistorted("./Resultados/ResultadosDistorted/ResultadosDistorted", resDistorted);


}

void medirTiemposDataset (const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k) {
    //codigo para calcular la cantidad de imagenes por persona suponiendo que las muestras son balanceadas y la cantidad de clases
    uint imagenesPorPersona = 0;
    uint cantidadDeClases = 0;
    uint personasTemp = 0;
    for( uint i = 0; i < labelsX.size(); i++){
        if (labelsX[i] != personasTemp){
            personasTemp = labelsX[i];
            cantidadDeClases++;
        }
        if(personasTemp == 1){
            imagenesPorPersona++;
        }
    }
//***********************************************************************//
    uint n = trainX.size()/imagenesPorPersona; //n es la cantidad de personas
    vector<int> folds(imagenesPorPersona);
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector<vector<double> > trainXTemp;
    vector<vector<double> > testYTemp;
    vector<uint> labelsXTemp;
    vector<uint> labelsYTemp;
    for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
        for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
            uint temp = (j*imagenesPorPersona)+u;
            if (u < imagenesPorPersona/k){ //si estoy en el fold que quiero
                testYTemp.push_back(trainX[temp]);
                labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
            } else{
                trainXTemp.push_back(trainX[temp]);
                labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
            }

        }
    }

//*********************Codigo para achicar data set de training ***************//
    uint imagenesPorPersonaEnFold = (k-1)*imagenesPorPersona/k;
    uint imagenesPorPersonaEnFoldInit = imagenesPorPersonaEnFold;
    vector<vector<unsigned long> > vectorTiemposYDataset (imagenesPorPersonaEnFoldInit);
    for (uint j = 0; j < imagenesPorPersonaEnFoldInit; ++j){//itero sobre la cantidad de imagenes por persona en el fold inicial, la idea es dejar siempre las muestras balanceadas, pero al final tener una de cada uno

        vector<vector<double>> V = PCATecho(trainXTemp,31);
        vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
        vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
        for (int i = 0; i < 20; i++) {
            unsigned long start, end;
            unsigned long delta = 0;
            RDTSC_START(start);
            vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,1);
            RDTSC_STOP(end);
            delta = end - start;//cada delta es el tiempo que tarda en calcular el PCA+ aplicar el cambio de base + calcular el knn para todos los elementos de testY
            vectorTiemposYDataset[imagenesPorPersonaEnFoldInit-1-j].push_back(delta);
        }

        for (uint u = 0; u < n; ++u){//itero sobre la cantidad de personas
            trainXTemp.erase(trainXTemp.begin()+((n-1-u)*imagenesPorPersonaEnFold)); //borro la ultima imagen de cada uno en cada paso
            labelsXTemp.erase(labelsXTemp.begin()+((n-1-u)*imagenesPorPersonaEnFold));
        }
        --imagenesPorPersonaEnFold;

    }
//*********************Codigo para achicar data set de training END***************//
    escribirTiemposDataset("./Resultados/TiemposVariandoTamDataset/TiemposVariandoTamDataset", vectorTiemposYDataset,n);
}

void escribirMatrizDeConfusion(string nombreArchivo, vector< vector<uint> > &matrizConfusion) {
    ofstream salida(nombreArchivo, ios_base::out);//getFlujo(nombreArchivo);
    for (vector< vector<uint> >::iterator it = matrizConfusion.begin() ; it != matrizConfusion.end(); ++it) {
        vector<uint>& fila = *it;
        string fila_str = "";
        for (vector<uint>::iterator it = fila.begin() ; it != fila.end(); ++it) {
            fila_str += to_string(*it) + "\t";
        }
        //cout << fila_str << endl;
        salida << fila_str << endl;
    }
    salida.close();
}


vector<vector< vector<uint> > > kfoldmatConf(const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k){
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
    uint n = trainX.size()/imagenesPorPersona;
    vector<int> folds(imagenesPorPersona);
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector< vector<int>> foldsXpersona (n);
    mt19937 g(static_cast<uint32_t>(time(0)));
    for (uint i = 0; i<n;++i){
        shuffle(folds.begin(), folds.end(),g);
        foldsXpersona[i] = folds;
    }
    vector<vector< vector<uint> > > matrices (0);
    for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
        vector<vector<double> > trainXTemp;
        vector<vector<double> > testYTemp;
        vector<uint> labelsXTemp;
        vector<uint> labelsYTemp;
        for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
            for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
                uint temp = (j*imagenesPorPersona)+foldsXpersona[j][u];
                if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
                    testYTemp.push_back(trainX[temp]);
                    labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
                } else{
                    trainXTemp.push_back(trainX[temp]);
                    labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
                }

            }
        }
        vector<vector<double>> V = PCATecho(trainXTemp,31);
        vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
        vector<vector<double> > testYTemp2 = testYTemp;
        multMatEsc(testYTemp2,0.5);
        testYTemp2 = multMat(testYTemp,V);
        vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,1);

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

//**********************Codigo para la matriz de confusion END**************// //Usar cuando ya tengamos parametros elegidos
        matrices.push_back(matrizConfusion);
    }
    return matrices;
}



void escribirDatosTamMatriz(string nombreArchivo, vector<pair<vector<resultados >,double> > &estadisticas, uint kDekfold, int variacion) {
    vector<resultados>* estadistica;
    int i = 328;
    for (vector<pair<vector<resultados >,double> >::iterator it = estadisticas.begin() ; it != estadisticas.end(); ++it) {
        vector<resultados >& estadistica = it->first;
        string accuracy = to_string(it->second);
        ofstream salida;
        string valor_parametro = int2stringConCantidadDigitos(4, i);
        salida = getFlujo(nombreArchivo + "Tamaño_" + valor_parametro);
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
        salida.close();
        i-=variacion; //le resto lo que fui variando el alpha

    }
}

void kfoldTamDataset(const vector<vector<double> >& trainX, const vector<clase_t>& labelsX, uint k){
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
    uint n = trainX.size()/imagenesPorPersona;
    vector<int> folds(imagenesPorPersona);
    for(uint i = 0; i < imagenesPorPersona; ++i){
        folds[i] = i;
    }
    vector< vector<int>> foldsXpersona (n);
    mt19937 g(static_cast<uint32_t>(time(0)));
    for (uint i = 0; i<n;++i){
        shuffle(folds.begin(), folds.end(),g);
        foldsXpersona[i] = folds;
    }
    for(uint i = 0; i<k; i++){ //itero sobre la cantidad de folds
        vector<vector<double> > trainXTemp;
        vector<vector<double> > testYTemp;
        vector<uint> labelsXTemp;
        vector<uint> labelsYTemp;
        for (uint j = 0; j < n; j++){ //itero sobre la cantidad de personas
            for (uint u = 0; u < imagenesPorPersona; u++) {//itero sobre la cantidad de imagenes por persona
                uint temp = (j*imagenesPorPersona)+foldsXpersona[j][u];
                if (u >= i*imagenesPorPersona/k && u < (i+1)*imagenesPorPersona/k){ //si estoy en el fold que quiero
                    testYTemp.push_back(trainX[temp]);
                    labelsYTemp.push_back(labelsX[temp]); //agrego el elemento a test
                } else{
                    trainXTemp.push_back(trainX[temp]);
                    labelsXTemp.push_back(labelsX[temp]); //agrego el elemento a train
                }

            }
        }
        vector<pair<vector<resultados >,double> > resVariandoDatasetUnFold;
//*********************Codigo para achicar data set de training ***************//
        uint imagenesPorPersonaEnFold = (k-1)*imagenesPorPersona/k;
        uint imagenesPorPersonaEnFoldInit = imagenesPorPersonaEnFold;
        for (uint j = 0; j < imagenesPorPersonaEnFoldInit; ++j){//itero sobre la cantidad de imagenes por persona en el fold inicial, la idea es dejar siempre las muestras balanceadas, pero al final tener una de cada uno

            vector<vector<double>> V = PCATecho(trainXTemp,31);
            vector<vector<double> > trainXTemp2 = multMat(trainXTemp,V);
            vector<vector<double> > testYTemp2 = multMat(testYTemp,V);
            vector<uint> vectordeKnns = vectorDeKnns(trainXTemp2,labelsXTemp,testYTemp2,1);
            vector<resultados > resultadosTemp (0);
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

            double accuracyTemp = accuracy(labelsYTemp,vectordeKnns);
            resVariandoDatasetUnFold.push_back(make_pair(resultadosTemp,accuracyTemp));
            for (uint u = 0; u < n; ++u){//itero sobre la cantidad de personas
                trainXTemp.erase(trainXTemp.begin()+((n-1-u)*imagenesPorPersonaEnFold)); //borro la ultima imagen de cada uno en cada paso
                labelsXTemp.erase(labelsXTemp.begin()+((n-1-u)*imagenesPorPersonaEnFold));
            }
            --imagenesPorPersonaEnFold;

        }
//*********************Codigo para achicar data set de training END***************//
        escribirDatosTamMatriz("./Resultados/ResultadosVariandoTamDataset/ResultadosVariandoTamDataset", resVariandoDatasetUnFold,i,n);//41 es la variacion del tamaño del dataset

    }
}
