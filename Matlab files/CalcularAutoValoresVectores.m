function [ V, D ] = CalcularAutoValoresVectores( Mx, accuracy, autocaras )
%CALCULARAUTOVALORESVECTORES Calcula V y D del TP
%   Utiliza el metodo de la potencia y calcula a V y D con la accuracy
%   dada, pero solo calcula la cantidad de autocaras pedidas.
    n = size(Mx);
    n = n(2);
    V = zeros(n,autocaras);
    D = zeros(n,n);
    MatrizIntermedia = Mx;
    for i = 1 : autocaras
        tic
        [autovalor, autovector] = MetodoPotencia(MatrizIntermedia, accuracy);
        MatrizIntermedia = Deflacion(MatrizIntermedia, autovalor, autovector);
        V(:, i) = autovector(:);
        D(i,i) = autovalor;
        i
        toc
    end
end

