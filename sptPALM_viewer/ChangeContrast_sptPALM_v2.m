function ChangeContrast_sptPALM_v2(h)

Up = get(h.Slider_UpperContrast, 'Value');
Low = get(h.Slider_LowerContrast, 'Value');

N = round(get(h.Slider_SelectFrame, 'Value'));
im = imread(h.MovieDisplayFullName, 'Index', N);

if isfield(h, 'MovieDisplayROI')
    im = imcrop(im, h.MovieDisplayROI);
end

Im_Corrected = imadjust(im,[Low Up],[]);

axes(h.MainAxes)
imshow(Im_Corrected)
axis image
colormap('Gray')