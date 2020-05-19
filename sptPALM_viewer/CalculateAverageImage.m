function AvIm = CalculateAverageImage

[ImageName, ImageDirectory] = uigetfile('*.tif', 'Select a movie associated to the MTT data you have load');
ImageFullName = strcat(ImageDirectory, ImageName);
ImInfo = imfinfo(ImageFullName);
NImages = length(ImInfo);

if NImages>100
    Dn = floor(NImages/100);
    Ntotal = 100;
else
    Dn = 1;
    Ntotal = NImages;
end


AvIm = uint32(zeros(ImInfo(1).Height, ImInfo(1).Width));

if NImages > 1
    fprintf('Calculating the average image ...     ')
    for nimage = 1 : Dn : NImages
        fprintf('\b\b\b\b%03i%%', round(100*nimage/NImages))
        im = imread(ImageFullName, 'Index', nimage);
        im = uint32(im);
        
        AvIm = AvIm + im;
    end
    AvIm = uint16(round(AvIm/Ntotal));
    fprintf('\r\n');
    
else
    AvIm = imread(ImageFullName);
end