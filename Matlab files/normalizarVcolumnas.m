function [ V_norm ] = normalizarVcolumnas( V )
%NORMALIZARV Summary of this function goes here
%   Detailed explanation goes here
    [n,m] = size(V);
    for i = 1:m
        V_norm(:,i) = V(:,i) / norm(V(:,i),2);
    end

end

