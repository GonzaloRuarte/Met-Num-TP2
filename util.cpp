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


void listarDirectorio(const string& directorio, stringvec& v)
{
    string nomArch;
    DIR* dirp = opendir(directorio.c_str());
    struct dirent * dp;
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
    int i = 0;

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


void verificarMatrizAImagen(string pathImagen, int cantidadDeImagenes, int alto, int ancho, vector<vector<double>> dataSet) {
    uchar* data = new uchar[cantidadDeImagenes*ancho*alto];
    int copiado = 0;
    for (int i=0; i<cantidadDeImagenes; i++) {
        vectorADatos(dataSet[i], copiado, data);
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

void getEtiquetas(stringvec *listaImagenes, map<int, string> etiquetas) {
    const char separador = '/';
    for (int i = 0; i < listaImagenes->size(); i++) {
        etiquetas[i] = explode((*listaImagenes)[i], separador)[2];

    }

}

vector<vector<double>>& cargarDataSetEnMatriz(string pathAlDataSet) {
    int ancho = 0;
    int alto = 0;
    vector<vector<double>> dataSet(0);
    stringvec listaImagenes;

    listarImagenes(pathAlDataSet, listaImagenes);
    //copy(listaImagenes.begin(), listaImagenes.end(), ostream_iterator<string>(cout, "\n"));
    map<int, string> *etiquetas = new map<int, string>;
    getEtiquetas(&listaImagenes, *etiquetas);

    string s = (etiquetas->find(44))->second;
    cout << s << endl;

    int tamanoDeReferencia = getTamanoImagenes(listaImagenes[0], &ancho, &alto);
    //cout << tamanoDeReferencia << endl;
    cargarDataSet(listaImagenes, tamanoDeReferencia, &dataSet);


    return dataSet;
}
