function [ alphas, tiempos ] = importarDatos( valork )
%IMPORTDATA Summary of this function goes here
%   Detailed explanation goes here
    tiempos = [];
    alphas = [];
    for i = 0 : 4
        tiempos = [tiempos; importdata(['TiemposVariandoAlpha', valork, '/TiemposVariandoAlphaK_', valork, 'Alpha_00', char(string(2*i)), '1'])];
        alphas = [alphas; ones(20,1) + (20*i)];
    end
    for i = 5 : 16
        tiempos = [tiempos; importdata(['TiemposVariandoAlpha', valork, '/TiemposVariandoAlphaK_', valork, 'Alpha_0', char(string(2*i)), '1']) ];
        alphas = [alphas; ones(20,1) + (20*i)];
    end
end

