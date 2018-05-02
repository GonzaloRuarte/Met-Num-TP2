function [ matriz_imagenes ] = importar_imagenes_chicas( tp2_folder )
%IMPORTAR_IMAGENES_CHICAS importa las imagenes a una matriz
%   Le das la carpeta donde estan y las mete cada una en una fila de la
%   matriz.
%   En mi caso tp2_folder es '~/Documents/MetNum/TP2/archivos_tp2/',
%   "~/" en matlab, hace que empiece la busqueda desde la carpeta home del 
%   sistema operativo ubuntu.
    X = importdata([tp2_folder, 'ImagenesCarasRed/s1/1.pgm']);
    X = reshape(X',1,[]);
    matriz_imagenes = X;
    for im_n = 2:10
        image_number = char(string(im_n)); %combierto el numero im_n en un string de numeros, y despues en un arreglo de chars.
        X = importdata([tp2_folder, 'ImagenesCarasRed/s1/', image_number, '.pgm']);
        X = reshape(X',1,[]);
        matriz_imagenes = [matriz_imagenes; X];
    end
    for suj_n = 2:41
        subject_number = char(string(suj_n));
        for im_n = 1:10
            image_number = char(string(im_n)); %combierto el numero im_n en un string de numeros, y despues en un arreglo de chars.
            X = importdata([tp2_folder, 'ImagenesCarasRed/s', subject_number, '/', image_number, '.pgm']);
            X = reshape(X',1,[]);
            matriz_imagenes = [matriz_imagenes; X];
        end
    end

end

