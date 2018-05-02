function [ vector_autoimagenes ] = Crear_n_Autoimagenes( tp2_folder, accuracy, n, sonImagenesChicas )
%CREAR_N_AUTOIMAGENES Importa las imagenes y construye los primeros n
%autovectores. Nota importante, para acceder a la imagen i se usa
%vector_autoimagenes(:,:,i).
%   Le das la carpeta donde están las imagenes y calcula Mx
%   En mi caso tp2_folder es '~/Documents/MetNum/TP2/archivos_tp2/' o 
%   '/home/leandro/Documents/MetNum/TP2/archivos_tp2/',
%   "~/" en matlab, hace que empiece la busqueda desde la carpeta home del 
%   sistema operativo ubuntu.
%   Realiza el método de la potencia solo para los n primeros autovectores
%   pedidos, y finalmente convierte aquellos autovectores, en matrices del
%   tamaño de las imagenes originales, y coloca esas matrices en un vector.
%   sonImagenesChicas es un bool, si es true, utiliza las imagenes chicas
%   en la carpeta dada, sino usa las grandes.
    if(sonImagenesChicas)
        Mx = CalcularMxChicas(tp2_folder);
    else
        Mx = CalcularMx(tp2_folder);
    end
    V = CalcularAutoValoresVectores(Mx,accuracy,n);
    matrix_width = size(Mx);
    matrix_width = matrix_width(2);
    if(matrix_width==10304)
        anchoim = 92;
        altoim = 112;
    else
        anchoim = 23;
        altoim = 28;
    end
    vector_autoimagenes = zeros(altoim, anchoim, n);
    for i = 1:n
        img = Autovector_a_Autoimagen( V(:,i), anchoim, altoim);
        vector_autoimagenes(:,:,i) = img; %guardo la img i en la iesima posicion
    end
end

