function [ alphas, tiempos ] = importarDatosFina( valork )
%IMPORTARDATOSFINA Summary of this function goes here
%   Detailed explanation goes here
    tiempos = [];
    alphas = [];
    for i = 1 : 9
        tiempos = [tiempos; importdata(['TiemposVariandoAlphaFina', valork, '/TiemposVariandoAlphaFinaK_', valork, 'Alpha_000', char(string(i))])];
        alphas = [alphas; ones(20,1) + (i-1)];
    end
    for i = 10 : 61
        tiempos = [tiempos; importdata(['TiemposVariandoAlphaFina', valork, '/TiemposVariandoAlphaFinaK_', valork, 'Alpha_00', char(string(i))]) ];
        alphas = [alphas; ones(20,1) + (i-1)];
    end
    
end

