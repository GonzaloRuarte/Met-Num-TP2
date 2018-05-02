function [ Mx2 ] = Deflacion( Mx, autovalor, autovector )
%DEFLACION Saca el autovector dado, para poder seguir usando Método de
%Potencia
%   Devuelve A-lambda*v*v'
    Mx2 = Mx - (autovalor*(autovector*autovector'));

end

