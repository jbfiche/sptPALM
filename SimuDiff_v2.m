function SimuDiff_v2(h)

%% Definition of the parameters used for the simulation
%% ====================================================

PixelSize = h.SimulationParameters.PixelSize; % in µm
AcquisitionTime = h.SimulationParameters.AcquisitionTime; % in s

NProteinsTot_Init = h.SimulationParameters.NProteinsTot_Init; % Define the total pool of proteins that can be activated and imaged
MeanProbActivation = h.SimulationParameters.MeanProbActivation; % The probability is defined as the mean number of protein activated per frame
MeanPhotoBleachingTime = h.SimulationParameters.MeanPhotoBleachingTime; % in s
Ton = h.SimulationParameters.Ton; % in s
Toff1 = h.SimulationParameters.Toff1; % in s
Toff2 = h.SimulationParameters.Toff2; % in s

MaxBlink = h.SimulationParameters.MaxBlink; % in frames
MinTrajLength = h.SimulationParameters.MinTrajLength; %in frames

ImageSize = h.SimulationParameters.ImageSize; % px
Offset = h.SimulationParameters.Offset; % intensity offset for the image
Noise1 = h.SimulationParameters.Noise1; % parameter lambda for the Poisson distribution describing the noise of the EMCCD camera
Noise2 = h.SimulationParameters.Noise2; % parameter describing the normal noise for each pixel and linked to the photon detection
GaussStampSize = h.SimulationParameters.GaussStampSize; % parameters defining the half-size of the window where the 2D gaussian will be calculated
LimitResolution = h.SimulationParameters.LimitResolution; % in µm
MeanSingleEventIntensity = h.SimulationParameters.MeanSingleEventIntensity; % Mean intensity of single emitters

Diff1 = str2double(get(h.Simulation_Diff1, 'String')); % µm²/s
Diff2 = str2double(get(h.Simulation_Diff2, 'String')); % µm²/s
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

verbose = get(h.Simulation_Verbose, 'Value');

%% Define the folder where the results of the analysis and the movies will
%% be saved
%% ========

FolderName = uigetdir;
cd(FolderName)

%% Start the simulation by calculating all the trajectories according to
%% parameters defined in the first part
%% ====================================

hwb = waitbar(0, 'Calculating the trajectories...');
NProteinsTot = NProteinsTot_Init;

Ton = Ton/AcquisitionTime;
Toff1 = Toff1/AcquisitionTime;
Toff2 = Toff2/AcquisitionTime;

for nFrame = 1 : NFrames
    
    waitbar(nFrame/NFrames);
    
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
        % proteins during its movement assuming a certain type of movement
        % and a typical diffusion coefficient
        % ------------------------------------
        
        Traj = zeros(LifeTime,4);
        Traj(1,:) = [nFrame, X, Y, 0];
        
        for dt = 1 : LifeTime-1
            %             d = abs(sqrt(4*Diff*AcquisitionTime)*random('Normal', 1, 0.2));
            d = sqrt(4*Diff*AcquisitionTime)/PixelSize;
            Theta = random('Uniform', 0, 360);
            Traj(dt+1,:) = [nFrame+dt, Traj(dt,2)+d*cosd(Theta), Traj(dt,3)+d*sind(Theta), d];
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
        % will be used to build the simulation movie is created
        % ---------------------------------------------------------
        
        Traj = Traj(EmissionState==1,:);
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
            Trajectories = cat(1, Trajectories, Traj(Frame_SubTraj(nSubTraj,1):Frame_SubTraj(nSubTraj,2),:));
        end
        
    end
end

close(hwb)

%% Open the figure where the results will be saved
%% ===============================================

hPlot = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen

%% Plot all the trajectories - the color-code is a function of the
%% detection/activation time
%% =========================

hold off
cla
hold on

Color = jet;
NTraj = size(Trajectories, 1);
ntraj_color = ceil(NTraj/size(Color,1));

for ntraj = 1 : NTraj
    
    X = Trajectories{ntraj}(:,2);
    Y = Trajectories{ntraj}(:,3);
    line(Y, X, 'Color', Color(ceil(ntraj/ntraj_color),:),'LineWidth',1)
    
end

axis square
axis off
box off
rectangle('Position', [1, 1, ImageSize, ImageSize])

saveas(gcf, 'Trajectories.fig', 'fig')
export_fig(hPlot, 'Trajectories.png');

%% Analyze the simulated trajectories using a modified version of the
%% function written for sptPALM_CBS_v2 and plot the results
%% ========================================================

MaxDisplayTime = 0.5; % (s) The maximum display time for the MSD curve
p = 4; % Only the p first points from the MSD curve are used for the calculation of the diffusion coefficient
MinNPoint = 0.75; % (%) minimum number of points for the trajectory
MinNPointMSD = 3; % Minimum number of point for an estimation of the MSD at a specific lag time
MaxStepLength = 10; % (frame) Maximum distance between two consecutive detections in the same trajectory
DiffCalculationMethod = 1; % 1=calculation based on the fit of the first "p" points - 2=calculate the weighted average

fileID = Simulated_Trajectory_analysis_v1(hPlot, FolderName, MaxBlink, MinNPoint, MinTrajLength, 1000*AcquisitionTime, PixelSize, MaxDisplayTime, p, MinNPointMSD, MaxStepLength, DiffCalculationMethod, Trajectories);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Diffusion coefficient of population 1(µm²/s)');
fprintf(fileID, '\n\n\n %4.5f\n', Diff1);
fprintf(fileID, '\r\n %s', 'Diffusion coefficient of population 2(µm²/s)');
fprintf(fileID, '\n\n\n %4.5f\n', Diff2);
fprintf(fileID, '\r\n %s', 'Ration of population 2 with respect to the overall number of proteins available');
fprintf(fileID, '\n\n\n %4.5f\n', PopulationRatio);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Image size (px)');
fprintf(fileID, '\n\n\n %4.2f\n', ImageSize);
fprintf(fileID, '\r\n %s', 'Number of frames');
fprintf(fileID, '\n\n\n %4.2f\n', NFrames);
fprintf(fileID, '\r\n %s', 'Intensity offset');
fprintf(fileID, '\n\n\n %4.2f\n', Offset);
fprintf(fileID, '\r\n %s', 'Mean EMCCD dark noise intensity (defined by a Poisson law)');
fprintf(fileID, '\n\n\n %4.2f\n', Noise1);
fprintf(fileID, '\r\n %s', 'Mean EMCCD amplification noise (defined by a Normal law)');
fprintf(fileID, '\n\n\n %4.2f\n', Noise2);
fprintf(fileID, '\r\n %s', 'Limit of resolution (in µm)');
fprintf(fileID, '\n\n\n %4.2f\n', LimitResolution);
fprintf(fileID, '\r\n %s', 'Mean intensity of the single fluorescent events');
fprintf(fileID, '\n\n\n %4.2f\n', MeanSingleEventIntensity);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Number of proteins available initially in the pool');
fprintf(fileID, '\n\n\n %4.2f\n', NProteinsTot_Init);
fprintf(fileID, '\r\n %s', 'Number of proteins effectively activated');
fprintf(fileID, '\n\n\n %4.2f\n', NProteinsTot_Init-NProteinsTot);
fprintf(fileID, '\r\n %s', 'Average number for the photo-activation probability');
fprintf(fileID, '\n\n\n %4.2f\n', MeanProbActivation);
fprintf(fileID, '\r\n %s', 'Average life time of emitters before photo-bleaching (s)');
fprintf(fileID, '\n\n\n %4.2f\n', MeanPhotoBleachingTime);
fprintf(fileID, '\r\n %s', 'Ton (average emission time in frames)');
fprintf(fileID, '\n\n\n %4.2f\n', Ton);
fprintf(fileID, '\r\n %s', 'Toff_1 (average off time in frames - short one)');
fprintf(fileID, '\n\n\n %4.2f\n', Toff1);
fprintf(fileID, '\r\n %s', 'Toff_2 (average off time in frames - long one)');
fprintf(fileID, '\n\n\n %4.2f\n', Toff2);

fclose(fileID);

% Plot the distribution of trajectories length
% --------------------------------------------

hold off
cla

Bin = min(LifeTime_all) : (max(LifeTime_all)-min(LifeTime_all))/50 : max(LifeTime_all);
[Counts, Bin] = hist(LifeTime_all, Bin);
plot(Bin, Counts/sum(Counts), '-ob')

ax = gca;
ax.FontSize = h.FontSize;
axis square
box on
xlabel('PhotoBleaching time (frames)')
ylabel('Frequency')
title('Distribution of photo-bleaching time for the experiment')

export_fig(hPlot, 'Trajectories_Length.png');

% Plot empirical single step length distribution for all the trajectories
% -----------------------------------------------------------------------

hold off
cla

[f1,x1] = ecdf(SingleStepLength*PixelSize);
plot(x1,f1, '-r', 'LineWidth', 1)

ax = gca;
ax.FontSize = h.FontSize;
axis square
axis([0 max(x1) 0 1])
box on
xlabel('Step length (nm)')
ylabel('Cumulative distribution')
title('Cumulative distribution of the step length')

export_fig(hPlot, 'Cumulative_distribution.png');

close(hPlot)

%% Create the simulated movie
%% ==========================

if verbose
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

% Initialize the template for the image according to the aquisition
% parameters
% ----------

im0 = uint16(Offset*ones(ImageSize*ImageSize,1));
[~,Idx] = sort(DetectionList(:,1));
DetectionList = DetectionList(Idx,:);

Gaussian_2D = @(x,y,x0,y0,s,A) A*exp(-(x-x0).^2/s^2).*exp(-(y-y0).^2/s^2)*random('Normal', 1, Noise2);
[Xmesh,Ymesh] = meshgrid(-GaussStampSize:1:GaussStampSize);
Xmesh = reshape(Xmesh, [(2*GaussStampSize+1)^2,1]);
Ymesh = reshape(Ymesh, [(2*GaussStampSize+1)^2,1]);

hwb = waitbar(0, 'Calculating the movies ...');
NMovies = ceil(NFrames/1000);

% According to the simulated trajectories, all the detection were saved in
% the array called "DetectionList". For each frame, the noise of the emCCD
% detector is simulated based on a Poisson distribution (Noise1) and the
% intensity detected for each emitter is simulated using the function
% "Gaussian_2D". Note that the intensity of each pixel is calculated using
% a 2D gaussian with a std equal to the diffraction limit but, for each
% pixel, a noise is also added to simulate the noise associated to the
% electronic amplification (Normal law).
% -------------------------------------

for nmovie = 1 : NMovies
    
    MovieName = strcat('AC', num2str(nmovie-1), '.tif');
    
    for nframe = 1000*(nmovie-1)+1 : 1000*nmovie
        
        waitbar(nframe/NFrames);
        
        im = im0 + uint16(random('Poisson', Noise1, [ImageSize*ImageSize,1]));
        im = reshape(im, [ImageSize, ImageSize]);
        
        EmittersList = DetectionList(DetectionList(:,1)==nframe, 2:3);
        for nemitter = 1 : size(EmittersList,1)
            
            if EmittersList(nemitter,1)>0 && EmittersList(nemitter,1)<ImageSize && EmittersList(nemitter,2)>0 && EmittersList(nemitter,2)<ImageSize
                
                Xpixel_mesh = round(EmittersList(nemitter,1)) + Xmesh;
                Ypixel_mesh = round(EmittersList(nemitter,2)) + Ymesh;
                A = random('Poisson', MeanSingleEventIntensity);
                Intensity = uint16(Gaussian_2D(Xpixel_mesh, Ypixel_mesh, EmittersList(nemitter,1), EmittersList(nemitter,2), LimitResolution/PixelSize, A));
                
                for npixel = 1 : size(Xpixel_mesh,1)
                    
                    if Xpixel_mesh(npixel)>0 && Xpixel_mesh(npixel)<ImageSize && Ypixel_mesh(npixel)>0 && Ypixel_mesh(npixel)<ImageSize
                        im(Xpixel_mesh(npixel), Ypixel_mesh(npixel)) = Intensity(npixel) + im(Xpixel_mesh(npixel), Ypixel_mesh(npixel));
                    end
                end
            end
        end
        
        ImRegistered = 0;
        while ~ImRegistered
            try
                if nframe == 1
                    imwrite(im, MovieName,'WriteMode','overwrite')
                    ImRegistered = 1;
                else
                    imwrite(im, MovieName,'WriteMode','append')
                    ImRegistered = 1;
                end
                
            catch exception
                ImRegistered = 0;
            end
        end
        
        if verbose
            axes(h.MainAxes)
            imagesc(im)
            axis image
            colorbar
            colormap('Gray')
            title(sprintf('Frame number %i', nframe));   
        end
    end
end

close(hwb)
if verbose
    set(h.sptPALM_DisplayMovie, 'Visible', 'off');
end

