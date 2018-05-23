function [ x_res, y_res ] = podarOutliers( x, y )
%PODAR les das puntos en x , con valores de x iguales, y poda entre
%los de y a los que se van muy lejos, para los iguales en x.
%   Detailed explanation goes here
    [n,m] = size(y);
    if (m>n)
        x = x';
        y = y';
    end
    n = max(n,m);
    
    n = n/20;
    x_res=[];
    y_res=[];
    x_temp = [];
    y_temp = [];
    for i = 0:n-1
        x_temp = x(i*20+1:(i+1)*20);
        y_temp = y(i*20+1:(i+1)*20);
        y_temp = sort(y_temp);
        while( iqr(y_temp)*1.5 < abs(max(y_temp) - mean(y_temp)) )
            [tam,tam2] = size(y_temp);
            tam = max(tam,tam2);
            x_temp = x_temp(1:tam-1);
            y_temp = y_temp(1:tam-1);
        end
        x_res = [x_res;x_temp];
        y_res = [y_res;y_temp];
    end

end

