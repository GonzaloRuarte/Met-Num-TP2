function [ autovalor, autovector ] = MetodoPotencia( Mx, accuracy )
%METODOPOTENCIA Calcula el mayor autovalor y mayor autovector
%   El algoritmo termina cuando la norma2 de la diferencia entre
%   Mx*v-lambda*v, su norma2 es menor a accuracy. Lo ideal es que su norma
%   2 sea 0 (<=> esa resta da vector 0, osea son iguales)
    pair = size(Mx);
    n = pair(2);
    autovector = rand(n,1); %Creo a un autovector del ancho de la matriz.
    autovector = autovector / norm(autovector); % Lo normalizo.
    autovalor = (autovector' * Mx * autovector)/ (autovector' * autovector);
    while (norm(Mx * autovector - autovalor*autovector) > accuracy)
        autovector = Mx * autovector;
        autovector = autovector / norm(autovector);
        autovalor = (autovector' * Mx * autovector)/ (autovector' * autovector);
    end
end

