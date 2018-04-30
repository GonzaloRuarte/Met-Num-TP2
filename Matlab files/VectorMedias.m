function [ vector_media ] = VectorMedias( matriz_muestra )
%VECTORMEDIAS hace el promedio de cada columna y lo guarda en un vector.
%   Detailed explanation goes here
    [alto, ancho] = size(matriz_muestra);
    vector_media = zeros(1, ancho);
    for i = 1:ancho
        for j = 1:alto
            vector_media(i) = vector_media(i) + matriz_muestra(j,i);
        end
    end
    vector_media = vector_media / alto;
    vector_media = vector_media';
end

