clc
clear all
close all

% Define the camera parameters
% ----------------------------

i = 0.5; % Intensity in photons
q = 0.9; % Quantum yield
c = 0; % Charges (in electrons) that are created during the transfer from the detector to the EM register
g = 200; % Electronic gain of the emCCD detector
f = 11.97; % A/D factor or detector sensitivity
r = 42; % Readout noise

% Define the intensity (count) probability function
% -------------------------------------------------

l = i*q + c;

P_pos = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic)^2/(2*r^2)) + 2/g * ncx2pdf(2*l,4,2*f*Nic/g);
P_neg = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic)^2/(2*r^2));

% Calculate the distribution according to the parameters defined above
% --------------------------------------------------------------------

Min_count = -1000;
Max_count = 5000;
Count = Min_count : 1 : Max_count;
pdf = zeros(1,Max_count-Min_count);
pdf = double(pdf);

for Nic = Min_count : 1 : Max_count
    
    n = Nic - Min_count + 1;
    if Nic>0
        pdf(n) = P_pos(Nic);
    else
        pdf(n) = P_neg(Nic);
    end
end

figure(1)
plot(Count, pdf, '-b')

% Infer the cumulative distribution
% ---------------------------------

cdf = zeros(1,Max_count-Min_count);
cdf = double(cdf);

for n = 1 : size(pdf,2)
    if n == 1
        cdf(n) = pdf(n);
    else
        cdf(n) = cdf(n-1) + pdf(n);
    end
end

cdf = cdf/max(cdf);

figure(2)
plot(Count, cdf, '-b')

% Infer the inverse cumulative distribution
% -----------------------------------------

[~, Idx] = max(pdf);
min_step = ( cdf(Idx+1) - cdf(Idx-1) ) /2;

Count = transpose(Count);
cdf = transpose(cdf);
icdf_range = min_step : min_step/2 : 1;
icdf_range = transpose(icdf_range);
icdf = zeros(numel(icdf_range),2);

for n = 1 : numel(icdf_range)
    
    idx = dsearchn(cdf,icdf_range(n));
    icdf(n,:) = [icdf_range(n), Count(idx)];
end

figure(3)
plot(icdf(:,1), icdf(:,2), '-b')
    