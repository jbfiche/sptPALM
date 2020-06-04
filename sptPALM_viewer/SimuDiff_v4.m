%*****************************
%
% SimuDiff_v4.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function is used to create simulated Brownian trajectories
% and test the analysis. The idea is to have a tool to test the different
% analysis methods and see how reproducible/ accurate the results can be.
% -------------------------------------------------------------------------
% Specific: The trajectories are only Brownian.
% -------------------------------------------------------------------------
% To fix:
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function SimuDiff_v4(h)

%% Definition of the parameters used for the simulation
%% ====================================================

PixelSize = h.SimulationParameters.PixelSize; % in um
AcquisitionTime = h.SimulationParameters.AcquisitionTime; % in s

NProteinsTot_Init = h.SimulationParameters.NProteinsTot_Init; % Define the total pool of proteins that can be activated and imaged
MeanPhotons = h.SimulationParameters.MeanPhotons; % Average number of incident photons
MeanProbActivation = h.SimulationParameters.MeanProbActivation; % The probability is defined as the mean number of protein activated per frame
MeanPhotoBleachingTime = h.SimulationParameters.MeanPhotoBleachingTime; % in s
Ton = h.SimulationParameters.Ton; % in s
Toff1 = h.SimulationParameters.Toff1; % in s
Toff2 = h.SimulationParameters.Toff2; % in s

MaxBlink = h.SimulationParameters.MaxBlink; % in frames
MinTrajLength = h.SimulationParameters.MinTrajLength; %in frames

ImageSize = h.SimulationParameters.ImageSize; % px
Offset = h.SimulationParameters.Offset; % intensity offset for the image
GaussStampSize = h.SimulationParameters.GaussStampSize; % parameters defining the half-size of the window where the 2D gaussian will be calculated
LimitResolution = h.SimulationParameters.LimitResolution; % in um

QY = h.SimulationParameters.QY; % quantum yield
Gain = h.SimulationParameters.Gain; % gain of the camera
CCDsensitivity = h.SimulationParameters.CCDsensitivity; % CCD sensitivity (e/AD counts)
ReadoutNoise = h.SimulationParameters.ReadoutNoise; % readout noise (e)

Diff1 = str2double(get(h.Simulation_Diff1, 'String')); % um^2/s
Diff2 = str2double(get(h.Simulation_Diff2, 'String')); % um^2/s
PopulationRatio = str2double(get(h.Simulation_Fraction, 'String'))/100; % Define the ratio of population 1 with respect to overall population of proteins available
NFrames = str2double(get(h.Simulation_NFrames, 'String'));

if isnan(Diff1)
    hwarn = warndlg('You need to indicate a non-zero value for the coefficient of diffusion #1');
    uiwait(hwarn)
    return
end

if isnan(Diff2)
    Diff2 = 0;
    PopulationRatio = 1;
end

Trajectories = {};
LifeTime_all = [];
TrajLength = [];
SingleStepLength = [];
DetectionList = [];

%% Define the folder where the results of the analysis and the movies will
%% be saved
%% ========

FolderName = uigetdir;
cd(FolderName)

%% Start the simulation by calculating all the trajectories according to
%% the parameters defined for the emitters.
%% ========================================

clc
fprintf('Calculating the trajectories ...     ')

NProteinsTot = NProteinsTot_Init;

Ton = Ton/AcquisitionTime;
Toff1 = Toff1/AcquisitionTime;
Toff2 = Toff2/AcquisitionTime;

for nFrame = 1 : NFrames
    
    fprintf('\b\b\b\b%03i%%', round(100*nFrame/NFrames))
    
    % Define how many FPs are photo-activated for this frame
    % ------------------------------------------------------
    
    N = random('Binomial', NProteinsTot, MeanProbActivation/NProteinsTot);
    NProteinsTot = NProteinsTot - N;
    
    % For each newly activated proteins, it defines its first detection
    % position as well as for how many frames the FPs will be detected.
    % -----------------------------------------------------------------
    
    for nProt = 1 : N
        
        X = random('Uniform', 1, ImageSize);
        Y = random('Uniform', 1, ImageSize);
        LifeTime = round(random('Exponential', round(MeanPhotoBleachingTime/AcquisitionTime)));
        LifeTime_all = cat(1, LifeTime_all, LifeTime);
        
        % Define whether the protein belongs to population 1 or 2
        % -------------------------------------------------------
        
        Pop = random('Binomial', 1, PopulationRatio);
        if Pop == 1
            Diff = Diff1;
        else
            Diff = Diff2;
        end
        
        % From it, it is possible to define all the positions of the
        % proteins during its movement assuming a Brownian motion and an
        % average diffusion coefficient.
        % The movement being brownian, the distance r is described by a
        % normal law N(0,sqrt(4Dt)). The x/y coordinates are then described
        % by a normal law of average 0 and standard deviation 2Dt. 
        %
        % Note that the values r,x,y are expressed in um.
        % -----------------------------------------------
        
        Traj = zeros(LifeTime,4);
        Traj(1,:) = [nFrame, X*PixelSize, Y*PixelSize, 0];
        
        for dt = 1 : LifeTime-1
            
            x = random('Normal', 0, sqrt(2*Diff*AcquisitionTime));
            y = random('Normal', 0, sqrt(2*Diff*AcquisitionTime));
            
            Traj(dt+1,:) = [nFrame+dt, ...
                Traj(dt,2)+x, ...
                Traj(dt,3)+y, ...
                abs(sqrt(x^2+y^2))];
        end
        
        % Using the values of ton, toff and the LifeTime, the ON/OFF
        % properties of the molecule are calculated
        % -----------------------------------------
        
        EmissionState = zeros(LifeTime,1);
        T = 1;
        Emitting = 1;
        
        while T < LifeTime
            if Emitting
                
                ton = round(random('Exponential', round(Ton)));
                if T+ton<LifeTime
                    EmissionState(T:T+ton) = 1;
                    SingleStepLength = cat(1, SingleStepLength, Traj(T+1:T+ton,4)); % The first point is not taken into account since the distance is at least calculated between the "n" and "n+1" points
                else
                    EmissionState(T:LifeTime) = 1;
                    SingleStepLength = cat(1, SingleStepLength, Traj(T+1:LifeTime,4));
                end
                T = T+ton+1;
                Emitting = 0;
            else
                
                % When the molecule is not emitting, we choose between the
                % two values of the off time (short = blink) or long (dark
                % state).
                
                if random('Binomial', 1, 0.5)<0.5
                    toff = round(random('Exponential', round(Toff1)));
                else
                    toff = round(random('Exponential', round(Toff2)));
                end
                T = T+toff+1;
                Emitting = 1;
            end
        end
        
        % From the array "EmissionState", the lists of detections that
        % will be used to build the simulation movie is created. Note that
        % all the detections detected after the maximum number of frames
        % are removed.
        % ------------
        
        Traj = Traj(EmissionState==1,:);
        Traj = Traj(Traj(:,1)<=NFrames,:);
        DetectionList = cat(1, DetectionList, Traj);
        
        % According to the values of "MaxBlink" and "MinTrajLength", the
        % trajectory "Traj" is splitted into several sub-trajectories that
        % are saved in the cell "Trajectories".
        % --------------------------------------
        
        Frame_SubTraj = [];
        First_Frame = 1;
        Npoint = 1;
        
        for nT = 2 : size(Traj,1)
            
            dt = Traj(nT,1)-Traj(nT-1,1);
            if dt > MaxBlink+1 && nT < size(Traj,1) && Npoint>=MinTrajLength
                
                Frame_SubTraj = cat(1, Frame_SubTraj, [First_Frame, nT-1, Npoint]);
                First_Frame = nT;
                Npoint = 1;
                
            elseif dt > MaxBlink+1 && nT < size(Traj,1) && Npoint<MinTrajLength
                
                First_Frame = nT;
                Npoint = 1;
                
            elseif nT == size(Traj,1) && Npoint+1>=MinTrajLength
                
                Frame_SubTraj = cat(1, Frame_SubTraj, [First_Frame, nT, Npoint+1]);
                
            else
                
                Npoint = Npoint+1;
            end
        end
        
        for nSubTraj = 1 : size(Frame_SubTraj,1)
            TrajLength = cat(1, TrajLength, Frame_SubTraj(nSubTraj,3));
            subtraj = Traj(Frame_SubTraj(nSubTraj,1):Frame_SubTraj(nSubTraj,2),:);
            Trajectories = cat(1, Trajectories, subtraj');
        end
    end
end

fprintf('\r\n')

%% Open the figure where the results will be saved
%% ===============================================

hPlot = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen

%% Plot all the trajectories - the color-code is a function of the
%% detection/activation time
%% =========================

fprintf('Plotting the trajectories ...     ');

Color = jet;
NTraj = size(Trajectories, 1);
ntraj_color = ceil(NTraj/size(Color,1));

for ntraj = 1 : NTraj
    
    fprintf('\b\b\b\b%03i%%', round(100*ntraj/NTraj))
    
    X = Trajectories{ntraj}(2,:)/PixelSize;
    Y = Trajectories{ntraj}(3,:)/PixelSize;
    line(Y, X, 'Color', Color(ceil(ntraj/ntraj_color),:),'LineWidth',1)
    
end

fprintf('\r')

axis square
axis off
box off
rectangle('Position', [1, 1, ImageSize, ImageSize])

saveas(hPlot, 'Trajectories.png');

%% Analyze the simulated trajectories using a modified version of the
%% function written for sptPALM_viewer and plot the results
%% ========================================================

MaxDisplayTime = 0.5; % (s) The maximum display time for the MSD curve
p = 4; % Only the p first points from the MSD curve are used for the calculation of the diffusion coefficient
MinNPoint = 0.75; % (%) minimum number of points for the trajectory
MinNPointMSD = 3; % Minimum number of point for an estimation of the MSD at a specific lag time

fileID = Simulated_Trajectory_analysis_v3(hPlot, ...
    FolderName, ...
    MaxBlink, ...
    MinNPoint, ...
    MinTrajLength, ...
    1000*AcquisitionTime, ...
    MaxDisplayTime, ...
    p, ...
    MinNPointMSD, ...
    Trajectories, ...
    PixelSize);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Diffusion coefficient of population 1(um^2/s)');
fprintf(fileID, '\n\n\n %4.5f\n', Diff1);
fprintf(fileID, '\r\n %s', 'Diffusion coefficient of population 2(um^2/s)');
fprintf(fileID, '\n\n\n %4.5f\n', Diff2);
fprintf(fileID, '\r\n %s', 'Ratio of population 2 with respect to the overall number of proteins available');
fprintf(fileID, '\n\n\n %4.5f\n', PopulationRatio);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Image size (px)');
fprintf(fileID, '\n\n\n %4.2f\n', ImageSize);
fprintf(fileID, '\r\n %s', 'Number of frames');
fprintf(fileID, '\n\n\n %4.2f\n', NFrames);
fprintf(fileID, '\r\n %s', 'Acquisition time (s)');
fprintf(fileID, '\n\n\n %4.2f\n', AcquisitionTime);
fprintf(fileID, '\r\n %s', 'Gain');
fprintf(fileID, '\n\n\n %4.2f\n', Gain);
fprintf(fileID, '\r\n %s', 'Intensity offset');
fprintf(fileID, '\n\n\n %4.2f\n', Offset);
fprintf(fileID, '\r\n %s', 'emCCD quantum yield QY');
fprintf(fileID, '\n\n\n %4.2f\n', QY);
fprintf(fileID, '\r\n %s', 'emCCD sensitivity (e/AD counts)');
fprintf(fileID, '\n\n\n %4.2f\n', CCDsensitivity);
fprintf(fileID, '\r\n %s', 'emCCD readout noise (e)');
fprintf(fileID, '\n\n\n %4.2f\n', ReadoutNoise);
fprintf(fileID, '\r\n %s', 'Limit of resolution (in um)');
fprintf(fileID, '\n\n\n %4.2f\n', LimitResolution);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Number of proteins available initially in the pool');
fprintf(fileID, '\n\n\n %4.2f\n', NProteinsTot_Init);
fprintf(fileID, '\r\n %s', 'Number of proteins effectively activated');
fprintf(fileID, '\n\n\n %4.2f\n', NProteinsTot_Init-NProteinsTot);
fprintf(fileID, '\r\n %s', 'Average number for the photo-activation probability');
fprintf(fileID, '\n\n\n %4.2f\n', MeanProbActivation);
fprintf(fileID, '\r\n %s', 'Average number of emitted photons per molecule');
fprintf(fileID, '\n\n\n %4.2f\n', MeanPhotons);
fprintf(fileID, '\r\n %s', 'Average life time of emitters before photo-bleaching (s)');
fprintf(fileID, '\n\n\n %4.2f\n', MeanPhotoBleachingTime);
fprintf(fileID, '\r\n %s', 'Ton (average emission time in frames)');
fprintf(fileID, '\n\n\n %4.2f\n', Ton);
fprintf(fileID, '\r\n %s', 'Toff_1 (average off time in frames - short one)');
fprintf(fileID, '\n\n\n %4.2f\n', Toff1);
fprintf(fileID, '\r\n %s', 'Toff_2 (average off time in frames - long one)');
fprintf(fileID, '\n\n\n %4.2f\n', Toff2);

fclose(fileID);

% Plot the distribution of trajectories duration
% ----------------------------------------------

hold off
cla

Bin = min(LifeTime_all) : (max(LifeTime_all)-min(LifeTime_all))/50 : max(LifeTime_all);
histogram(LifeTime_all, Bin, 'Normalization', 'probability');

ax = gca;
ax.FontSize = h.FontSize;
axis square
box on
xlabel('PhotoBleaching time (frames)')
ylabel('Frequency')
title('Distribution of photo-bleaching time for the experiment')

saveas(hPlot, 'Trajectories_Length.png');

% Plot empirical single step length distribution for all the trajectories
% -----------------------------------------------------------------------

hold off
cla

[f1,x1] = ecdf(SingleStepLength);
plot(x1,f1, '-r', 'LineWidth', 1)

ax = gca;
ax.FontSize = h.FontSize;
axis square
axis([0 max(x1) 0 1])
box on
xlabel('Step length (um)')
ylabel('Cumulative distribution')
title('Cumulative distribution of the step length')

saveas(hPlot, 'Cumulative_distribution.png');
close(hPlot)

%% Save the list of all the positions in a txt file
%% ================================================

save('Detection_list.txt', 'DetectionList', '-ascii');

%% Create the simulated movie
%% ==========================

if get(h.Simulation_Verbose, 'Value') == 1
    set(h.sptPALM_DisplayMovie, 'Visible', 'on');
end

% Define the folder where the movies are going to be saved
% --------------------------------------------------------

cd(FolderName)
LookForSimulationFolder = dir('Simulated_Movies');
if ~isempty(LookForSimulationFolder)
    rmdir('Simulated_Movies', 's');
end

mkdir(FolderName, 'Simulated_Movies');
cd('Simulated_Movies')

% Calculate the noise simulation for the emCCD detector according to the
% parameters defined by the user
% ------------------------------

emCCD_noise_distribution = emCCD_signal_modelization(h);

% Convert the positions saved in the variable "DetectionList" in pixel
% --------------------------------------------------------------------

DetectionList(:,2:3) = DetectionList(:,2:3)/PixelSize;

% Initialize the template for a single emitter
% --------------------------------------------

Gaussian_2D = @(A,x,y,x0,y0,s) A/(s*sqrt(2*pi)) * exp(-(x-x0).^2/s^2).*exp(-(y-y0).^2/s^2);
StampSize = (2*GaussStampSize+1)^2;

% According to the simulated trajectories, all the detections were saved in
% the array called "DetectionList". For each frame, the noise of the emCCD
% detector is simulated based on a Poisson distribution (Noise1) and the
% intensity detected for each emitter is simulated using the function
% "Gaussian_2D". Note that the intensity of each pixel is calculated using
% a 2D gaussian with a std equal to the diffraction limit but, for each
% pixel, a noise is also added to simulate the noise associated to the
% electronic amplification (Normal law).
% -------------------------------------

NMovies = ceil(NFrames/1000);

for nmovie = 1 : NMovies
    
    fprintf('Calculating the movies # %i ...     ', nmovie)
    MovieName = strcat('Simulation_', num2str(nmovie-1), '.tif');
    
    for nframe = 1 : 1000
        
        fprintf('\b\b\b\b%03i%%', round(100*nframe/1000))
        Frame = 1000*(nmovie-1)+nframe;
        
        if Frame <= NFrames
            
            % Initialize the image
            % --------------------
            
            I_photon = zeros(ImageSize*ImageSize,1);
            I_counts = zeros(ImageSize*ImageSize,1);
            
            % Calculate the position of each emitter on the image and the
            % number of photons for each pixel
            % --------------------------------
            
            EmittersList = DetectionList(DetectionList(:,1)==Frame, 2:3);
            
            for nemitter = 1 : size(EmittersList,1)
                x = EmittersList(nemitter,1);
                y = EmittersList(nemitter,2);
                
                [X,Y] = meshgrid(x-GaussStampSize:1:x+GaussStampSize , y-GaussStampSize:1:y+GaussStampSize);
                X = reshape(round(X), [StampSize,1]);
                Y = reshape(round(Y), [StampSize,1]);
                Idx_within_im_border = find(X>0 & Y>0 & X<ImageSize & Y<ImageSize);
                
                if ~isempty(Idx_within_im_border)
                    Photons = Gaussian_2D(random('Poisson', MeanPhotons), X(Idx_within_im_border), Y(Idx_within_im_border), x, y, LimitResolution/PixelSize);
                    LinearIdx = sub2ind([ImageSize,ImageSize], X(Idx_within_im_border), Y(Idx_within_im_border));
                    I_photon(LinearIdx) = I_photon(LinearIdx) + Photons;
                end
            end
            
            I_photon = round(I_photon);
            
            % Using the model developped by Hirsch et al. for emCCD, calculate
            % the A/D counts for each pixel
            % -----------------------------
            
            Unique_Intensity = unique(I_photon);
            
            for n_value = 1 : length(Unique_Intensity)
                
                nphoton = Unique_Intensity(n_value);
                Idx = find(I_photon == nphoton);
                N = random('Uniform', 0, 1, [length(Idx),1]);
                I_counts(Idx) = emCCD_noise_distribution{nphoton+1}(N);
            end
            
            I_counts = round(I_counts);
            
            % Create the final simulated image and save it
            % --------------------------------------------
            
            I_counts = uint16(Offset + I_counts);
            im = reshape(I_counts, [ImageSize, ImageSize]);
            
            ImRegistered = 0;
            while ~ImRegistered
                try
                    if Frame == 1
                        imwrite(im, MovieName,'WriteMode','overwrite')
                        ImRegistered = 1;
                    else
                        imwrite(im, MovieName,'WriteMode','append')
                        ImRegistered = 1;
                    end
                    
                catch 
                    ImRegistered = 0;
                end
            end
            
            if get(h.Simulation_Verbose, 'Value') == 1
                axes(h.MainAxes)
                imagesc(im)
                axis image
                colorbar
                colormap('Gray')
                title(sprintf('Frame number %i', Frame));
            end
        end
    end
    fprintf('\r')
end