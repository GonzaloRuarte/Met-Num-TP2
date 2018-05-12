from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from itertools import izip_longest

# Cada archivo contiene todos los folds en los que se realizo la medicion para ese valor de parametro
# La lectura toma todas las mediciones y las promedia

def leer_estadisticas_de_archivo(nombre_archivo):
    estadistica = {}
    estadistica['accuracy'] = 0.0
    estadistica['precision'] = []
    estadistica['recall'] = []
    estadistica['f1'] = []
    with open(nombre_archivo, 'r') as f:
        contenido = f.readlines()
        cantidad_de_folds = len(contenido) / 4
        for i in range(0, len(contenido), 4):
            estadistica['accuracy']  = estadistica['accuracy'] + float(contenido[i])            
            estadistica['precision'] = [sum(n) for n in izip_longest(estadistica['precision'], [float(x) for x in contenido[i+1].split()], fillvalue=0.0)]
            estadistica['recall']    = [sum(n) for n in izip_longest(estadistica['recall'], [float(x) for x in contenido[i+2].split()], fillvalue=0.0)]
            estadistica['f1']        = [sum(n) for n in izip_longest(estadistica['f1'], [float(x) for x in contenido[i+3].split()], fillvalue=0.0)]
    
    estadistica['accuracy']  = estadistica['accuracy'] / cantidad_de_folds
    estadistica['precision'] = [n / cantidad_de_folds for n in  estadistica['precision']]
    estadistica['recall'] =    [n / cantidad_de_folds for n in  estadistica['recall']]
    estadistica['f1'] =        [n / cantidad_de_folds for n in  estadistica['f1']]
    return estadistica

def leer_estadisticas(ruta):
    estadisticas = []
    nombre_de_archivos = listdir(ruta)
    nombre_de_archivos.sort()
    for nombre_archivo in nombre_de_archivos:
        estadisticas.append(leer_estadisticas_de_archivo(ruta + nombre_archivo))
    return estadisticas    
    
def get_recall(clase, estadisticas):
    ret = []
    for estadistica in estadisticas: 
        ret.append(estadistica['recall'][clase])
    return ret



def get_precision(clase, estadisticas):
    ret = []
    for estadistica in estadisticas:
        ret.append(estadistica['precision'][clase])
    return ret

def get_promedio(estadisticas, tipo):
    ret = []
    for estadistica in estadisticas: 
        suma = 0
        for n in estadistica[tipo]:
            suma = suma + n 
        promedio = suma / len(estadistica[tipo])
        ret.append(promedio)
    return ret

def graficar_estadisticas(presicion, recall, titulo, nombre_archivo):
    plt.clf()
    plt.step(recall, precision, color='b', alpha=0.2, where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title(titulo)
    fig = plt.gcf()
    plt.show()
    plt.draw()
    fig.savefig(nombre_archivo, dpi=600)

estadisticas = leer_estadisticas('/home/christian/Resultados/')
recall = get_promedio(estadisticas, 'recall')
precision = get_promedio(estadisticas, 'precision')
    
graficar_estadisticas(precision, recall, 
                      'curva Precision-Recall variando el alfa multiclase premediada', 
                      '/home/christian/graficosResultados/variandoElAlfa.png')


