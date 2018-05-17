# -*- coding: utf-8 -*-
import sys
from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
import math
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from itertools import izip_longest
import collections
import random

# Cada archivo contiene todos los folds en los que se realizo la medicion para ese valor de parametro
# La lectura toma todas las mediciones y las promedia

print sys.stdout.encoding

def leer_estadisticas_de_archivo(nombre_archivo):
    estadistica = {}
    estadistica['accuracy'] = 0.0
    estadistica['precision'] = []
    estadistica['recall'] = []
    estadistica['f1'] = []
    i=0;
    with open(nombre_archivo, 'r') as f:
        contenido = f.readlines()
        cantidad_de_folds = len(contenido) / 4
        for i in range(0, len(contenido), 4):
            #print contenido[i+1].split()[0]
            estadistica['accuracy']  = estadistica['accuracy'] + float(contenido[i])            
            estadistica['precision'] = [sum(n) for n in izip_longest(estadistica['precision'], 
                                                    [float(x) for x in contenido[i+1].split()], fillvalue=0.0)]
            estadistica['recall']    = [sum(n) for n in izip_longest(estadistica['recall'], 
                                                    [float(x) for x in contenido[i+2].split()], fillvalue=0.0)]
            estadistica['f1']        = [sum(n) for n in izip_longest(estadistica['f1'], 
                                                    [float(x) for x in contenido[i+3].split()], fillvalue=0.0)]
    
    estadistica['accuracy']  = estadistica['accuracy'] / cantidad_de_folds
    estadistica['precision'] = [n / cantidad_de_folds for n in  estadistica['precision']]
    estadistica['recall'] =    [n / cantidad_de_folds for n in  estadistica['recall']]
    estadistica['f1'] =        [n / cantidad_de_folds for n in  estadistica['f1']]
    return estadistica

def leer_estadisticas(ruta):
    estadisticas = []
    nombre_de_archivos = listdir(ruta)
    nombre_de_archivos.sort()
    etiquetas = []
    for nombre_archivo in nombre_de_archivos:
        etiquetas.append(int(nombre_archivo.split('_')[1]))
        estadisticas.append(leer_estadisticas_de_archivo(ruta + nombre_archivo))
    return {'estadisticas': estadisticas, 'etiquetas': etiquetas}
    
def obtener(clase, estadisticas, tipo):
    ret = []
    for estadistica in estadisticas: 
        ret.append(estadistica[tipo][clase])
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

def getPaleta(nombre, cantidad_de_colores):
    colores = []
    mapa_de_color = plt.cm.get_cmap(nombre)
    for i in np.arange(cantidad_de_colores):
        j = random.uniform(0, 1)
        while j in colores:
            j = random.uniform(0, 1)
        colores.append(mapa_de_color(j))
    return colores

# Obtenemos el maximo valor maximiza el parametro medido para cada clase.
def get_valores_que_maximiza_parametro_por_clase(estadisticas, tipo):
    cantidad_de_clases = len(estadisticas[0][tipo])
    valores_del_parametro = [0] * cantidad_de_clases
    valores_maximos = [0.0] * cantidad_de_clases
    valor_del_parametro = 0;
    for estadistica in estadisticas:
        estadistica_actual = estadistica[tipo]
        j = 0;
        for n in estadistica_actual:
            if (valores_maximos[j] < n):
                valores_del_parametro[j] = valor_del_parametro
                valores_maximos[j] = n
                
            j = j + 1
        valor_del_parametro = valor_del_parametro + 1

    return [valores_del_parametro, valores_maximos]
    
def graficar_41_clases(directorio_estadisticas, tipo, titulo, nombre_archivo):
    plt.clf()
    plt.figure(figsize=(10,6))
    
    estadisticas = leer_estadisticas(directorio_estadisticas)
    ests = estadisticas['estadisticas']
    
    colores = getPaleta('Dark2', 41)

    with plt.style.context('default'):
        for i in np.arange(40):
            j = random.sample(np.arange(10000),  1)[0]
            plt.plot(estadisticas['etiquetas'], obtener(i, ests, tipo))

    plt.title(titulo)
    fig = plt.gcf()
    ax = plt.gca()
    ax.set_xticks(estadisticas['etiquetas'])
    ax.grid(color='gray', linestyle='dashed', linewidth=1)
    plt.show()
    plt.draw()
    fig.savefig(nombre_archivo, dpi=600)
    
    
def graficar_estadisticas(precision, recall, f1, etiquetasX, titulo, nombre_archivo):
    plt.clf()
    plt.figure(figsize=(10,6))
    
    plt.plot(etiquetasX, precision, 'r', label = 'precision') 
    plt.plot(etiquetasX, recall, 'b', label = 'recall')  
    plt.plot(etiquetasX, f1, 'g', label = 'f1') 
    plt.legend()
    
    plt.xlim([1, 321])
    plt.title(titulo)
    fig = plt.gcf()
    ax = plt.gca()
    ax.set_xticks(etiquetasX)
    ax.grid(color='gray', linestyle='dashed', linewidth=1)
    plt.show()
    plt.draw()
    fig.savefig(nombre_archivo, dpi=600)


    
    
def graficar(directorio_estadisticas, titulo_grafico, nombre_archivo_imagen):
    estadisticas = leer_estadisticas(directorio_estadisticas)
#    get_valores_que_maximiza_parametro_por_clase(estadisticas['estadisticas'], 'recall') 
#    get_valores_que_maximiza_parametro_por_clase(estadisticas['estadisticas'], 'precision') 

    recall = get_promedio(estadisticas['estadisticas'], 'recall')
    precision = get_promedio(estadisticas['estadisticas'], 'precision')   
    f1 = get_promedio(estadisticas['estadisticas'], 'f1')   
    clase = 2
    graficar_estadisticas(precision, recall, f1, estadisticas['etiquetas'], 
                          titulo_grafico, nombre_archivo_imagen)    


    
def graficar_estadisticas_barras(valores_maximizados, etiquetas, titulo, nombre_archivo, indice):
    plt.clf()
    n = len(valores_maximizados[0])
    
    cmapBlues=plt.cm.Blues
    #cmapDivering=plt.cm.viridis
    colores = getPaleta('Dark2', n)    

    plt.bar(np.arange(n), valores_maximizados[1], color=colores)  # Dibujamos el gráfico de barras
    #plt.ylim(550,650)  # Limitamos los valores del eje y al range definido [450, 550]
    plt.title(titulo)  # Colocamos el título
    #cmap=plt.cm.Blues
    etqs = []
    i = 0
    for k in valores_maximizados[0]:
        etqs.append(etiquetas[k] + ' (clase ' + str(i + indice + 1) + ')')
        i = i + 1
    plt.xticks(np.arange(n), etqs, rotation = 90)  # Colocamos las etiquetas del eje x
    fig = plt.gcf()
    plt.show()
    plt.draw()   
    fig.subplots_adjust(wspace=0.1, top=0.9, right=0.9, left=0.1, bottom=0.30)
    fig.savefig(nombre_archivo, dpi=600)

def dividir(valores, indice_inicial, indice_final):
    return [valores[0][indice_inicial:indice_final], valores[1][indice_inicial:indice_final]]
    

def graficar_barras(directorio_estadisticas, tipo, titulo_grafico, nombre_archivo_imagen):
    estadisticas = leer_estadisticas(directorio_estadisticas)
    ests = estadisticas['estadisticas']
    etiquetas = estadisticas['etiquetas']
    valores_maximizados = get_valores_que_maximiza_parametro_por_clase(ests, tipo)
    valores_maximizados1 = dividir(valores_maximizados, 0, 13)
    valores_maximizados2 = dividir(valores_maximizados, 14, 27)
    valores_maximizados3 = dividir(valores_maximizados, 27, 40)
    graficar_estadisticas_barras(valores_maximizados1, etiquetas, titulo_grafico, nombre_archivo_imagen + '_1.png', 0)
    graficar_estadisticas_barras(valores_maximizados2, etiquetas, titulo_grafico, nombre_archivo_imagen + '_2.png', 14)
    graficar_estadisticas_barras(valores_maximizados3, etiquetas, titulo_grafico, nombre_archivo_imagen + '_3.png', 27)
    
#graficar_barras('/home/christian/Resultados/', 'recall',
#                      'valores del parametro k que maximizan el recall',
#                      '/home/christian/graficosResultados/barras_k_recall')

#graficar_barras('/home/christian/Resultados/', 'precision',
#                      'valores del parametro k que maximizan el precision',
#                      '/home/christian/graficosResultados/barras_k_precision')

graficar('/home/christian/Resultados/', 
                      'precision recall y f1 promediados, variando el parametro k',
                      '/home/christian/graficosResultados/precision_recall_promediados_k.png')
       
graficar_41_clases('/home/christian/Resultados/', 
                   'recall', 'recall de las 41 clases variando el parametro k', 
                   '/home/christian/graficosResultados/recall_41_clases_k.png')
graficar_41_clases('/home/christian/Resultados/', 'precision', 
                   'precision de las 41 clases variando el parametro k', 
                   '/home/christian/graficosResultados/precision_41_clases_k.png')
graficar_41_clases('/home/christian/Resultados/', 'f1', 
                   'f1 de las 41 clases variando el parametro k', 
                   '/home/christian/graficosResultados/f1_41_clases_k.png')
#graficar('/home/christian/Resultados/', 
#         'curva Precision-Recall variando el alfa multiclase premediada',
#         '/home/christian/graficosResultados/variandoElAlfa.png')

