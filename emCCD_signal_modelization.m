clc
clear all
close all


% Define the camera parameters
% ----------------------------

i = 0; % Intensity in photons
q = 0.9; % Quantum yield
c = 1; % Charges (in electrons) that are created during the transfer from the detector to the EM register
g = 100; % Electronic gain of the emCCD detector
f = 63; % A/D factor or detector sensitivity
r = 72; % Readout noise

% Define the intensity (count) probability function
% -------------------------------------------------

l = i*q + c;

P_pos = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic)^2/(2*r^2)) + 2/g * ncx2pdf(2*l,4,2*f*Nic/g);
P_neg = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic)^2/(2*r^2));

% Calculate the distribution according to the parameters defined above
% --------------------------------------------------------------------

Min_count = -1000;
Max_count = 14000;
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

% Recalculate the range of count values according to the values of the pdf
% ------------------------------------------------------------------------

pdf = pdf / sum(pdf);
Idx_Min = find(pdf>0.0001, 1, 'first');
Idx_Max = find(pdf>0.0001, 1, 'last');

Count = Count(Idx_Min:Idx_Max);
pdf = pdf(Idx_Min:Idx_Max);
NCount = size(Count,2);

figure(1)
subplot(1,3,1)
plot(Count, pdf, '-b')
axis square
xlabel('A/D count')
ylabel('pdf')

% Infer the cumulative distribution
% ---------------------------------

cdf = zeros(1,NCount);
cdf = double(cdf);

for n = 1 : NCount
    if n == 1
        cdf(n) = pdf(n);
    else
        cdf(n) = cdf(n-1) + pdf(n);
    end
end

cdf = cdf/max(cdf);

subplot(1,3,2)
plot(Count, cdf, '-ob')
axis square
xlabel('A/D count')
ylabel('cdf')

% Infer the inverse cumulative distribution
% -----------------------------------------

fitobject = fit(cdf',Count','linearinterp');

X_fit = 0 : 0.001 : 1;
Y_fit = fitobject(X_fit);

subplot(1,3,3)
plot(cdf, Count, '-ob')
hold on
plot(X_fit, Y_fit, '-r')
axis square