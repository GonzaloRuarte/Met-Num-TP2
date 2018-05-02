function [ img ] = Autovector_a_Autoimagen( autovector, anchoim, altoim )
%AUTOVECTOR_A_AUTOIMAGEN revierte la transformacion de imagen a vector
%hecha en la funcion "importar_imagenes"
%   Usa Reshape de forma poco intuitiva pero esta bien (mirar en matlab wiki o
%   probar uno mismo) y finalmente lo transpone para que quede de las
%   mismas dimesiones que la imagen original.
    img = reshape(autovector, anchoim, altoim); %parece al revez, pero esto revierte el reshape hecho en importar imagenes
    img = img'; %ahora la hacemos una matriz de altoim filas x anchoim columnas, quedando correctamente.
end

