function [ matriz_semivarianza ] = MatrizSemivarianza( matriz_muestra, vector_media )
%MATRIZSEMIVARIANZA le resta a cada fila de matriz_muestra el vector_media.
%   Después de restar a cada fila, a ese resultado lo divide como escalar
%   la raiz cuadrada de n-1.
    [n, m] = size(matriz_muestra);
    matriz_semivarianza = zeros(n,m);
    for i = 1:n
        for j = 1:m
            matriz_semivarianza(i,j) = (matriz_muestra(i,j) - vector_media(j));
        end
    end
end

