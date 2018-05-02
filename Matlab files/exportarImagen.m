function [  ] = exportarImagen( autoimagenes, out_folder )
%EXPORTARIMAGEN agarra la hipermatriz de imagenes y exporta las imagenes a
%la carpeta dada como parametro out_folder, las imagenes van a ser "i.jpg"
%con i el numero del autovector.
%   Lo que dice arriba.
    n = size(autoimagenes);
    n = n(3);
    for i = 1 : n
        image_number = char(string(i));
        im_intermedio = autoimagenes(:,:,i);
        im_intermedio = im_intermedio - min(min(im_intermedio));
        im_intermedio = im_intermedio * (255 / max(max(im_intermedio)));
        uint_im = uint8(im_intermedio);
        imwrite(uint_im,[out_folder, image_number, '.jpg']);
    end
end

