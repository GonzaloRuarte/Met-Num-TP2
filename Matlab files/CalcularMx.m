function [ Mx ] = CalcularMx( tp2_folder )
%CALCULARMX Le das la carpeta de las imagenes(en mi caso 
%'~/Documents/MetNum/TP2/archivos_tp2/', o 
%'/home/leandro/Documents/MetNum/TP2/archivos_tp2/'), y te calcula la matriz Mx
%del TP2
%   Usa las funciones importar imagenes, matriz Semivarianza y
%   VectorMedias.
    imgs = double(importar_imagenes(tp2_folder));
    Vm = VectorMedias(imgs);
    X = MatrizSemivarianza(imgs, Vm);
    [n, m] = size(X);
    raiz = sqrt(n-1);
    X = X / raiz;
    Mx = X' * X;
end

