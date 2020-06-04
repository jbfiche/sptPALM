%*****************************
%
% emCCD_signal_modelization.m
%
% ****************************
%
% JB Fiche
% Mar, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function is calculating the inverse cumulative distribution
% of the A/D counts for an emCCD camera assuming the model described in
% Hirsch et al. Plos ONE, 2013. The output is an array containing the
% interpolation of the inverse distribution for a given number of incident
% photons (from 0 to a maximum value).
% -------------------------------------------------------------------------
% Specific:
% -------------------------------------------------------------------------
% To fix:
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function emCCD_noise_distribution = emCCD_signal_modelization(h)

% Define the camera parameters
% ----------------------------

MeanPhotons = h.SimulationParameters.MeanPhotons; % Average number of incident photons
q = h.SimulationParameters.QY; % quantum yield
g = h.SimulationParameters.Gain; % gain of the camera
f = h.SimulationParameters.CCDsensitivity; % CCD sensitivity (e/AD counts)
r = h.SimulationParameters.ReadoutNoise; % readout noise (e)
c = 5; % Spurious charges (values estimated experimentally)

verbose = 0;

emCCD_noise_distribution = cell(MeanPhotons,1);

fprintf('\r\n')
fprintf('---------------------------- \r\n')
fprintf('Calculating the noise simulation for the emCCD ...     ')

for n_photons = 0 : MeanPhotons
    
    fprintf('\b\b\b\b%03i%%', round(100*n_photons/MeanPhotons))
    
    % Define the intensity (count) probability function
    % -------------------------------------------------
    
    l = n_photons*q + c;
    
    P_pos = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic).^2/(2*r^2)) + 2/g * ncx2pdf(2*l,4,2*f*Nic/g);
    P_neg = @(Nic) 1/(sqrt(2*pi)*r) * exp(-l - (f*Nic).^2/(2*r^2));
    
    % Calculate the distribution according to the parameters defined above
    % --------------------------------------------------------------------
    
    Min_count = -1000;
    Max_count = 14000;
    Count_neg = Min_count : 1 : 0;
    pdf_neg = P_neg(Count_neg);
    
    Count_pos = 1 : 1 : Max_count;
    pdf_pos = P_pos(Count_pos);
    
    Count = cat(2, Count_neg, Count_pos);
    pdf = cat(2, pdf_neg, pdf_pos);
    
    % Recalculate the range of count values according to the values of the pdf
    % ------------------------------------------------------------------------
    
    pdf = pdf / sum(pdf);
    Idx_Min = find(pdf>0.0001, 1, 'first');
    Idx_Max = find(pdf>0.0001, 1, 'last');
    
    Count = Count(Idx_Min:Idx_Max);
    pdf = pdf(Idx_Min:Idx_Max);
    NCount = size(Count,2);
    
    if verbose
        figure(2)
        subplot(1,3,1)
        plot(Count, pdf, '-b')
        axis square
        xlabel('A/D count')
        ylabel('pdf')
    end
    
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
    
    if verbose
        subplot(1,3,2)
        plot(Count, cdf, '-ob')
        axis square
        xlabel('A/D count')
        ylabel('cdf')
    end
    
    % Infer the inverse cumulative distribution
    % -----------------------------------------
    
    fitobject = fit(cdf',Count','linearinterp');
    
    if verbose
        X_fit = 0 : 0.001 : 1;
        Y_fit = fitobject(X_fit);
        subplot(1,3,3)
        plot(cdf, Count, '-ob')
        hold on
        plot(X_fit, Y_fit, '-r')
        axis square
    end
    
    % Save the results of the simulations
    % -----------------------------------
    
    emCCD_noise_distribution{n_photons+1} = fitobject;
end

fprintf('\r\n')