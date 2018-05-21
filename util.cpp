#include <string>
#include <iostream>
#include <fstream>
#include "ppmloader.h"
#include "dirent.h"
#include <vector>
#include <stdio.h>
#include <map>
#include <sys/stat.h>
#include <sstream>

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
            throw runtime_error("no se puede cargar el archivo " + nombreArchivo + "!");
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

void getEtiquetas(stringvec *listaDirectorios, vector<uint>* etiquetas) {
		uint i = 0;
		for (vector<string>::iterator it = listaDirectorios->begin() ; it != listaDirectorios->end(); ++it) {
			i++;
			string path = *it;
			stringvec listadoImagenes;
        		listarDirectorio(path, listadoImagenes);
			for (vector<string>::iterator it = listadoImagenes.begin() ; it != listadoImagenes.end(); ++it) {
        			(*etiquetas).push_back(i);
			}
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
/*
 * nombreArchivo:   el nombre del test, pero sin el '.in' ni el '.out' el metodo se va a encargar de cargar los datos de ambos.
 * dataSet:         la matriz que contendra las imagenes de entrada.
 * labels:          los labels de las imagenes cargadas en dataSet.
 * autovalores:     los 15 autovalores de mayor magnitud de la matriz de covarianza Mx (v1, v2, ..., v15) ordenados decrecientemente.
 * */
void cargarTest(string nombreArchivo, vector<vector<double>> *dataSet, vector<uint> *labels, vector<double> *autovalores) {
    fstream entrada(nombreArchivo + ".in", ios_base::in);

    vector<string>* listaImagenes = new vector<string>(0);

    string lectura;
    bool path = true;
    while(entrada >> lectura) {
        lectura = explode(lectura, ',').at(0);
        if (path) {
            listaImagenes->push_back("./" + lectura);
            path = false;
        } else {
            labels->push_back(stoi(lectura));
            path = true;
        }
    }
    entrada.close();

    int tamanoDeReferencia = getTamanoImagenes(listaImagenes->at(0), &ancho, &alto);
    cargarDataSet(*listaImagenes, tamanoDeReferencia, dataSet);

    fstream salida(nombreArchivo + ".expected", ios_base::in);
    while(salida >> lectura) {
        autovalores->push_back(stod(lectura));
    }
    salida.close();

}

/*
 * nombreArchivo:   el nombre del archivo
 * dataSet:         la matriz que contendra las imagenes de entrada.
 * labels:          los labels de las imagenes cargadas en dataSet.
 * */
void cargarSet(string nombreArchivo, vector<vector<double>> *dataSet, vector<uint> *labels) {
    fstream entrada(nombreArchivo, ios_base::in);

    vector<string>* listaImagenes = new vector<string>(0);

    string lectura;
    bool path = true;
    while(entrada >> lectura) {
        lectura = explode(lectura, ',').at(0);
        if (path) {
            listaImagenes->push_back("./" + lectura);
            path = false;
        } else {
            labels->push_back(stoi(lectura));
            path = true;
        }
    }
    entrada.close();

    int tamanoDeReferencia = getTamanoImagenes(listaImagenes->at(0), &ancho, &alto);
    cargarDataSet(*listaImagenes, tamanoDeReferencia, dataSet);
}

/*
 * nombreArchivo:   el nombre del archivo
 * dataSet:         la matriz que contendra las imagenes de entrada.
 * labels:          los labels de las imagenes cargadas en dataSet.
 * */
void cargarSet(string nombreArchivo, vector<vector<double>> *dataSet) {
    fstream entrada(nombreArchivo, ios_base::in);

    vector<string>* listaImagenes = new vector<string>(0);

    string lectura;
    while(entrada >> lectura) {
        listaImagenes->push_back("./" + lectura);
    }
    entrada.close();

    int tamanoDeReferencia = getTamanoImagenes(listaImagenes->at(0), &ancho, &alto);
    cargarDataSet(*listaImagenes, tamanoDeReferencia, dataSet);
}

void guardarClasificacion(string nombreArchivo, vector<uint> &clasificacion) {
    ofstream salida(nombreArchivo, ios_base::out);
    for (vector<uint>::iterator it = clasificacion.begin() ; it != clasificacion.end(); ++it) {
        uint clase = *it;
        string clase_str = "";
        clase_str += to_string(*it) + "\t";
        //cout << fila_str << endl;
        salida << clase_str << endl;
    }
    salida.close();
}

bool existeArchivo(const string& nombreArchivo) {
    struct stat buf;
    bool ret = false;
    if (stat(nombreArchivo.c_str(), &buf) != -1) {
        ret = true;
    }
    return ret;
}

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void crearCarpetaSiNoExiste(string path) {
    if (!is_file_exist(path.c_str())) {
        mkdir(path.c_str(), S_IRWXU);
    }
}

void crearUltimaCarpetaSiNoExiste(string nombreArchivo) {
    vector<string> v{explode(nombreArchivo, '/')};
    string carpeta = "";
    for (int i=0; i<v.size()-1; i++)
        carpeta += v[i] + "/";
    crearCarpetaSiNoExiste(carpeta);
}

/*
 * Esta funcion es para que si el archivo ya existe lo abra y devuelve el flujo para agregarle datos
 * y si no existe lo crea y devuelve el flujo agregarle datos.
 */
ofstream getFlujo(const string& nombreArchivo) {
    crearUltimaCarpetaSiNoExiste(nombreArchivo);
    if (existeArchivo(nombreArchivo)) {
        ofstream ret(nombreArchivo, ios_base::app);
        return ret;
    } else {
        ofstream ret2(nombreArchivo, ios_base::out);
        return ret2;
    }
}

string int2stringConCantidadDigitos(int cantidadDigitos, int i){
    stringstream flujo;
    flujo.fill  ('0');
    flujo.width (cantidadDigitos);
    flujo << i;
    return(flujo.str());
}



