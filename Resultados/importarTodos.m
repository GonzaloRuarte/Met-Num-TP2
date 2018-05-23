function [ alphasM, tiemposM, alphasfinasM, tiemposfinasM ] = importarTodos( )
%IMPORTARTODOS Summary of this function goes here
%   Detailed explanation goes here
    tiemposM = [];
    alphasM = [];
    for i = 1:4
        if (i == 1)
            text = '1';
        elseif(i==2)
            text='20';
        elseif(i==3)
            text = '40';
        else
            text = '100';
        end
        [alphasM(i,:),tiemposM(i,:)] = importarDatos(text);
    end

    for i = 1:4
        if (i == 1)
            text = '1';
        elseif(i==2)
            text = '20';
        elseif(i==3)
            text = '40';
        else
            text = '100';
        end
        [alphasfinasM(i,:),tiemposfinasM(i,:)] = importarDatosFina(text);
    end
end

