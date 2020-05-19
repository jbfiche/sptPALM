%% In order to compare the theoretical results from the modelization (Hirsch
%% et al. 2013. Plos One), several acquisitions were performed with our
%% emCCD (ANDOR 897). The idea is to compare the experimental values for the
%% noise with the simulation and check whether the two are similar.
%% ================================================================

clear all
close all
clc

% First test was performed in background conditions, with an electronic
% gain of 0
% ---------

% cd('/mnt/PALM_dataserv/DATA/JB/2020/DATA_Alexandre/Test_emCCD/')
cd('/home/fiche/Desktop/Alexandre_data/DATA_Alexandre/Test_emCCD/')

imName = 'gain_0_50ms_background.tif';
Nframes = imfinfo(imName);
Nframes = length(Nframes)-1;

E = zeros(Nframes,1);
Pixel_values_distribution = cell(Nframes,1);

for nframe = 1 : Nframes
    
    im = imread(imName, nframe);
    im = reshape(im, [numel(im),1]);
    im = double(im);
    E(nframe) = mean(im);
    pixel = unique(im);
    values_distribution = zeros(size(pixel,1),2);
    
    if nframe == 1
        Min = min(pixel);
        Max = max(pixel);
    end
    if Min>min(pixel)
        Min = min(pixel);
    end
    if Max<max(pixel)
        Max = max(pixel);
    end
    
    for nvalues = 1 : size(pixel,1)
        values_distribution(nvalues,:) = [pixel(nvalues), sum(im(:)==pixel(nvalues))];
    end
    Pixel_values_distribution{nframe} = values_distribution;
end

Offset = round(mean(E));
NPixelValues = Max - Min + 1;
Pixel_values_distribution_final = zeros(NPixelValues,2);
Pixel_values_distribution_final(:,1) = Min : 1 : Max;

for nframe = 1 : Nframes
    
    Pixel_values = Pixel_values_distribution{nframe};
    for nvalues = 1 : size(Pixel_values,1)
        Idx = Pixel_values(nvalues,1)-Min+1;
        Pixel_values_distribution_final(Idx,2) = Pixel_values_distribution_final(Idx,2) + Pixel_values(nvalues,2);
    end
end

cdf_0 = Pixel_values_distribution_final;
for nvalues = 1 : size(Pixel_values_distribution_final,1)
    cdf_0(nvalues,2) = sum(Pixel_values_distribution_final(1:nvalues,2));
end
cdf_0(:,1) = cdf_0(:,1) - Offset;
cdf_0(:,2) = cdf_0(:,2)/cdf_0(end,2);


figure(1)
plot(cdf_0(:,1), cdf_0(:,2), '-ob')
axis tight
axis square
xlabel('A/D counts')
ylabel('cdf')
legend({'gain=0'})

% Second test was performed in background conditions, with an electronic
% gain of 100
% -----------

cd('/mnt/PALM_dataserv/DATA/JB/2020/DATA_Alexandre/Test_emCCD/')
imName = 'gain_100_50ms_background.tif';
Nframes = imfinfo(imName);
Nframes = length(Nframes)-1;

Pixel_values_distribution = cell(Nframes,1);

for nframe = 1 : Nframes
    
    im = imread(imName, nframe);
    im = reshape(im, [numel(im),1]);
    im = double(im);
    pixel = unique(im);
    values_distribution = zeros(size(pixel,1),2);
    
    if nframe == 1
        Min = min(pixel);
        Max = max(pixel);
    end
    if Min>min(pixel)
        Min = min(pixel);
    end
    if Max<max(pixel)
        Max = max(pixel);
    end
    
    for nvalues = 1 : size(pixel,1)
        values_distribution(nvalues,:) = [pixel(nvalues), sum(im(:)==pixel(nvalues))];
    end
    Pixel_values_distribution{nframe} = values_distribution;
end

NPixelValues = Max - Min + 1;
Pixel_values_distribution_final = zeros(NPixelValues,2);
Pixel_values_distribution_final(:,1) = Min : 1 : Max;

for nframe = 1 : Nframes
    
    Pixel_values = Pixel_values_distribution{nframe};
    for nvalues = 1 : size(Pixel_values,1)
        Idx = Pixel_values(nvalues,1)-Min+1;
        Pixel_values_distribution_final(Idx,2) = Pixel_values_distribution_final(Idx,2) + Pixel_values(nvalues,2);
    end
end

cdf_1 = Pixel_values_distribution_final;
for nvalues = 1 : size(Pixel_values_distribution_final,1)
    cdf_1(nvalues,2) = sum(Pixel_values_distribution_final(1:nvalues,2));
end
cdf_1(:,1) = cdf_1(:,1) - Offset;
cdf_1(:,2) = cdf_1(:,2)/cdf_1(end,2);


figure(2)
cla
hold off
plot(cdf_0(:,1), cdf_0(:,2), '-ob')
hold on
plot(cdf_1(:,1), cdf_1(:,2), '-og')
axis tight
axis square
xlabel('A/D counts')
ylabel('cdf')
legend({'gain=0', 'gain=100'})

% Second test was performed in background conditions, with an electronic
% gain of 200
% -----------

cd('/mnt/PALM_dataserv/DATA/JB/2020/DATA_Alexandre/Test_emCCD/')
imName = 'gain_200_50ms_background.tif';
Nframes = imfinfo(imName);
Nframes = length(Nframes)-1;

Pixel_values_distribution = cell(Nframes,1);

for nframe = 1 : Nframes
    
    im = imread(imName, nframe);
    im = reshape(im, [numel(im),1]);
    im = double(im);
    pixel = unique(im);
    values_distribution = zeros(size(pixel,1),2);
    
    if nframe == 1
        Min = min(pixel);axis tight
axis square
xlabel('A/D counts')
ylabel('cdf')
legend({'gain=0', 'gain=100', 'gain=200'})
        Max = max(pixel);
    end
    if Min>min(pixel)
        Min = min(pixel);
    end
    if Max<max(pixel)
        Max = max(pixel);
    end
    
    for nvalues = 1 : size(pixel,1)
        values_distribution(nvalues,:) = [pixel(nvalues), sum(im(:)==pixel(nvalues))];
    end
    Pixel_values_distribution{nframe} = values_distribution;
end

NPixelValues = Max - Min + 1;
Pixel_values_distribution_final = zeros(NPixelValues,2);
Pixel_values_distribution_final(:,1) = Min : 1 : Max;

for nframe = 1 : Nframes
    
    Pixel_values = Pixel_values_distribution{nframe};
    for nvalues = 1 : size(Pixel_values,1)
        Idx = Pixel_values(nvalues,1)-Min+1;
        Pixel_values_distribution_final(Idx,2) = Pixel_values_distribution_final(Idx,2) + Pixel_values(nvalues,2);
    end
end

cdf_2 = Pixel_values_distribution_final;
for nvalues = 1 : size(Pixel_values_distribution_final,1)
    cdf_2(nvalues,2) = sum(Pixel_values_distribution_final(1:nvalues,2));
end
cdf_2(:,1) = cdf_2(:,1) - Offset;
cdf_2(:,2) = cdf_2(:,2)/cdf_2(end,2);

figure(3)
cla
hold off
plot(cdf_0(:,1), cdf_0(:,2), '-ob')
hold on
plot(cdf_1(:,1), cdf_1(:,2), '-og')
plot(cdf_2(:,1), cdf_2(:,2), '-or')
axis tight
axis square
xlabel('A/D counts')
ylabel('cdf')
legend({'gain=0', 'gain=100', 'gain=200'})