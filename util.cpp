#include <string>
#include <iostream>
#include "ppmloader.h"
#include "dirent.h"
#include <vector>
#include <stdio.h>
#include <map>

#include "util.h"


using namespace std;

typedef vector<string> stringvec;

typedef unsigned char uchar;

using namespace std;

int ancho = 0;
int alto = 0;

void listarDirectorio(const string& directorio, stringvec& v)
{
    string nomArch;
    DIR* dirp = opendir(directorio.c_str());
    struct dirent * dp;
    if (dirp == NULL) {
        throw runtime_error("no se encontro directorio " + directorio + "!");
    }
    while ((dp = readdir(dirp)) != NULL) {
        string nomArch = dp->d_name;
        if (nomArch.compare(".") != 0 && nomArch.compare("..") != 0)
            v.push_back(directorio + "/" + nomArch);
    }
    closedir(dirp);
}

vector<double> datosAVector(uchar &datos, int tamano) {
    vector<double> ret (0);
    for (size_t i = 0; i < tamano; ++i) {
        double bt = static_cast<double>(*((unsigned char *)&datos + i));
        ret.push_back(bt);
    }
    return ret;
}

void vectorADatos(vector<double> vec, int cantidadCopiada, uchar* datos) {
    size_t i = 0;
    for (vector<double>::iterator it = vec.begin() ; it != vec.end(); ++it) {
        double byte = *it;
        *(datos + cantidadCopiada + i) = static_cast<unsigned char>(byte);
        i++;
    }
}

void cargarDataSet(stringvec listaImagenes, int tamanoDeReferencia, vector<vector<double>> *dataSet) {
    int ancho = 0, alto = 0;
    uchar* datos = NULL;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
    uint i = 0;

    for (stringvec::iterator it = listaImagenes.begin() ; it != listaImagenes.end(); ++it) {
        string nombreArchivo = *it;
        bool ret = LoadPPMFile(&datos, &ancho, &alto, &pt, nombreArchivo.c_str());
        if (!ret || ancho == 0|| alto == 0|| pt!=PPM_LOADER_PIXEL_TYPE_GRAY_8B) {
            throw runtime_error("no se puede cargar el archivo" + nombreArchivo + "!");
        }
        if (tamanoDeReferencia != ancho*alto)
            throw runtime_error("La imagen " + nombreArchivo + " difiere en tamano o profundidad de color y canales.");

        vector<double> vec = datosAVector(*datos, tamanoDeReferencia);
        (*dataSet).push_back(vec);
    }
}


void listarImagenes(string pathDataSet, stringvec& listadoImagenes) {
    stringvec listadoDirectorios;
    listarDirectorio(pathDataSet, listadoDirectorios);
    for (vector<string>::iterator it = listadoDirectorios.begin() ; it != listadoDirectorios.end(); ++it) {
        string path = *it;
        listarDirectorio(path, listadoImagenes);
    }
}


int getTamanoImagenes(string pathImagen, int *ancho, int *alto) {
    uchar *datos;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
    bool ret = LoadPPMFile(&datos, ancho, alto, &pt, pathImagen.c_str());
    if (!ret || ancho == 0|| alto == 0|| pt!=PPM_LOADER_PIXEL_TYPE_GRAY_8B) {
        throw runtime_error("no se puede cargar el archivo" + pathImagen + "!");
    }
    return (*ancho)*(*alto);
}


void convertirMatrizAImagen(string pathImagen, int cantidadDeImagenes, vector<vector<double>>* dataSet) {
    uchar* data = new uchar[cantidadDeImagenes*ancho*alto];
    int copiado = 0;
    for (int i=0; i<cantidadDeImagenes; i++) {
        vectorADatos((*dataSet)[i], copiado, data);
        copiado += ancho*alto;
    }

    SavePPMFile(pathImagen.c_str(), data, ancho, alto*cantidadDeImagenes, PPM_LOADER_PIXEL_TYPE_GRAY_8B, "salida junta");
}


const vector<string> explode(const string& s, const char& c)
{
    string buff{""};
    vector<string> v;

    for(auto n:s)
    {
        if(n != c) buff+=n; else
        if(n == c && buff != "") { v.push_back(buff); buff = ""; }
    }
    if(buff != "") v.push_back(buff);

    return v;
}

void getEtiquetas(stringvec *listaImagenes, vector<uint>* etiquetas) {
    for (uint i = 0; i < listaImagenes->size(); i++) {
        (*etiquetas)[i] = i;

    }

}

void cargarDataSetEnMatriz(string pathAlDataSet, vector<vector<double>>* dataSet, vector<uint>* labelsX) {
    stringvec listaImagenes;
	  stringvec listaDirectorios;

    listarDirectorio(pathAlDataSet, listaDirectorios);

	  listarImagenes(pathAlDataSet, listaImagenes);

    //copy(listaImagenes.begin(), listaImagenes.end(), ostream_iterator<string>(cout, "\n"));
    getEtiquetas(&listaDirectorios, labelsX);

    /*string s = (etiquetas->find(44))->second;
    cout << s << endl;*/

    int tamanoDeReferencia = getTamanoImagenes(listaImagenes[0], &ancho, &alto);
    //cout << tamanoDeReferencia << endl;
    cargarDataSet(listaImagenes, tamanoDeReferencia, dataSet);

    //codigo que escribe 20 de las imagenes cargagas en el vector para corroborar que las imagenes se
    // pasaron correctamente a la matriz

    string pathImagen = "./nomsal.pgm";
}



//------------------------ Parseo de la entrada -------------------------------//

bool contiene(char *argv[], const string *cadena) {
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[5];
    return param1.compare(*cadena) || param2.compare(*cadena) || param3.compare(*cadena) || param4.compare(*cadena);
}

string obtener(char *argv[], const string *cadena) {
    string ret;
    string param1 = argv[1], param2 = argv[3], param3 = argv[5], param4 = argv[5];
    if (param1.compare(*cadena)) ret = argv[2];
    if (param2.compare(*cadena)) ret = argv[4];
    if (param3.compare(*cadena)) ret = argv[6];
    if (param4.compare(*cadena)) ret = argv[7];
    return ret;
}

bool obtenerParametros(int argc, char * argv[], string *metodo, string *trainSet, string *testSet, string *classif) {
    bool ret = false;
    const string param1 = "-m", param2 = "-i", param3 = "-q", param4 = "-o";

    if (argc == 9 && contiene(argv, &param1) && contiene(argv, &param2) && contiene(argv, &param3) && contiene(argv, &param4)) {
        *metodo = obtener(argv, &param1);
        *trainSet = obtener(argv, &param2);
        *testSet = obtener(argv, &param3);
        *classif = obtener(argv, &param4);
        ret = (metodo != NULL && trainSet != NULL && testSet != NULL && classif != NULL);
    }
    return ret;
}
