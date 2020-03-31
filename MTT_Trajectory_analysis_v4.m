%*****************************
%
% MTT_Trajectory_analysis_v4.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function is used to calculate the MSD and instant diffusion
% coefficient of each trajectories.
% -------------------------------------------------------------------------
% Specific:
% -------------------------------------------------------------------------
% To fix:
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function h = MTT_Trajectory_analysis_v4(h)

%% Parameter for the analysis
%% ==========================

MaxBlinks = str2double(get(h.MaxBlinks, 'String')); % Maximum number of frames two successive detections can be separated of
MinNPoint = str2double(get(h.MinNumberPoints, 'String')); % Minimum number of frames for each trajectory
MinTrajLength_MSDCalculation = str2double(get(h.MinTrajLength, 'String')); % Minium number of frames for the trajectories used for the calculation of the MSD
AcquisitionTime = str2double(get(h.AcquisitionTime, 'String')); % in ms
PixelSize = str2double(get(h.PixelSize, 'String')); %  in um
MaxDisplayTime = str2double(get(h.MaxDisplayTime, 'String')); % The MSD curve will be displayed only from zero to this value in s
p = str2double(get(h.NumberPointsMSDFit, 'String')); % Use the first "p" points of the MSD to estimate the Dapp
MinNPointMSD = str2double(get(h.MinimumNumberPointsMSD, 'String')); % Minimum number of points used to calculate each values of the MSD
MaxStepLength = str2double(get(h.MaxStepLength, 'String')); % Minimum number of points used to calculate each values of the MSD

FontSize = h.FontSize;
DiffCalculationMethod = get(h.DiffusionCalculationMethod, 'Value');
Reconstructed_Traj = h.Reconstructed_Traj;

Plot_Traj = get(h.Plot_trajectories, 'Value');

%% Make sure the display window is visible
%% =======================================

ax = h.MainAxes;
set(h.sptPALM_DisplayMovie, 'Visible', 'on');

%% Remove the results from the previous analysis from structure "h" and
%% reinitialize the structure "Results"
%% ====================================

h = sptPALM_initialize(h, 'Reset_h');

if ~isempty(dir(h.ResultsFileName))
    delete(h.ResultsFileName)
end

Results = matfile(h.ResultsFileName, 'Writable', true);
Results.FileToAnalyse = h.FileToAnalyse;
Results.DirectoryName = h.DirectoryName;
Results.SingleStep_Length = h.SingleStep_Length;
Results.Length_Traj = h.Length_Traj;
Results.PixelSize = PixelSize;
Results.AcquisitionTime = AcquisitionTime;
Results.Reconstructed_Traj = Reconstructed_Traj;

%% Filter the trajectories according to the parameters selected
%% ============================================================

[Reconstructed_Traj_Filtered, NTraj_Filter] = Filter_Trajectories(Reconstructed_Traj, MaxBlinks, MinTrajLength_MSDCalculation, MinNPoint, MaxStepLength);
set(h.NTrajectoriesFiltered, 'String', num2str(NTraj_Filter)); % Display the # of trajectories selected after applying the filters
h.Reconstructed_Traj_Filtered = Reconstructed_Traj_Filtered;

%% Look whether we want to define an ROI. If a ROI is defined, return
%% only the trajectories detected within the ROI
%% =============================================

[h, Area, Reconstructed_Traj_ROI, NTraj_ROI] = Select_Trajectories_ROI(h, Results, ax, Reconstructed_Traj_Filtered, PixelSize);
set(h.NTrajectoriesROI, 'String', num2str(NTraj_ROI));  % Display the # of trajectories detected within the ROI
h.Reconstructed_Traj_ROI = Reconstructed_Traj_ROI;

%% Analyze the trajectories and calculate the MSD for each trajectory
%% (added as the fifth row)
%% ========================

% if h.Parallel_computing.Value == 0
%     [MSD_all,MSD_weight,Reconstructed_Traj_MSD] = MSD_calculation(Reconstructed_Traj_ROI,NTraj_ROI,MinNPointMSD,p);
% else
%     [MSD_all,MSD_weight,Reconstructed_Traj_MSD] = MSD_calculation_parallel_computing(Reconstructed_Traj_ROI,NTraj_ROI,MinNPointMSD,p);
% end

[MSD_all,MSD_weight,Reconstructed_Traj_MSD] = MSD_calculation(Reconstructed_Traj_ROI,NTraj_ROI,MinNPointMSD,p);

NTraj_MSD = size(MSD_all,1);
set(h.NTrajectoriesMSD, 'String', num2str(NTraj_MSD)); % Display the # of trajectories accepted for MSD calculation
h.Reconstructed_Traj_MSD = Reconstructed_Traj_MSD;

%% Return the track density and display it on the GUI
%% ==================================================

if isfield(h, 'AvIm') && exist('Area','var')~=0
    Density = round(1000*NTraj_ROI/(sum(Area)*PixelSize^2))/1000;
    set(h.TrackDensity, 'String', num2str(Density)); % Display the density of tracks per �m�
    h.Density = Density;
elseif isfield(h, 'AvIm') && exist('Area', 'var')==0
    Density = round(1000*NTraj_ROI/(size(h.AvIm,1)*size(h.AvIm,2)*PixelSize^2))/1000;
    set(h.TrackDensity, 'String', num2str(Density)); % Display the density of tracks per �m�
    h.Density = Density;
else
    Density = NaN;
    set(h.TrackDensity, 'String', ''); % Display the density of tracks per �m�
end

%% The apparent diffusion coefficient is calculated for each single trajectory
%% using the first "p" points of the MSD curves
%% ============================================

% For each point, the apparent coefficient diffusion is calculated
% using either a linear fit on the first p points of the average of the MSD
% for a lagtime of 1 frame
% -------------------------

% if h.Parallel_computing.Value == 0
%     [MSD_all,Reconstructed_Traj_MSD_accepted,Dapp] = Diff_calculation(DiffCalculationMethod,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD);
% else
%     [MSD_all,Reconstructed_Traj_MSD_accepted,Dapp] = Diff_calculation_parallel_computing(DiffCalculationMethod,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD);
% end

[MSD_all,Reconstructed_Traj_MSD_accepted,Dapp] = Diff_calculation(DiffCalculationMethod,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD);

NTraj_Diff = size(Dapp,1);
set(h.NTrajectoriesDiff, 'String', num2str(NTraj_Diff)); % Display the # of trajectories selected for D calculation
h.Dapp = Dapp;

%% Plot the distribution of apparent diffusion coefficient
%% =======================================================

LogDapp = log10(Dapp);
[NbrGaussianFit, D_mean, varargout] = FitGaussianDistribution_v2(LogDapp, MSD_all, FontSize, ax, round(MaxDisplayTime*1000/AcquisitionTime), Reconstructed_Traj_MSD_accepted, DiffCalculationMethod);
MSD_FIT = varargout{1};

T = cat(2, LogDapp(:,1), Dapp(:,1));
T = array2table(T, 'VariableNames', {'logD', 'D'});
writetable(T, 'Saved_Diffusion_Coeff.txt');

%% Plot the MSD curve
%% ==================

Lag = 1 : 1 : size(MSD_FIT,1);

axes(ax)
hold off
cla

if NbrGaussianFit == 2
    
    MSD_1 = MSD_FIT(:,1:3);
    MSD_2 = MSD_FIT(:,4:6);
    
    hold on
    errorbar(Lag*AcquisitionTime/1000, MSD_1(:,1), MSD_1(:,2), '-o', 'Color', [1 0.5 0])
    errorbar(Lag*AcquisitionTime/1000, MSD_2(:,1), MSD_2(:,2), '-o', 'Color', [0 0.4 1])
    errorbar(Lag*AcquisitionTime/1000, MSD_1(:,3), MSD_1(:,2), '-s', 'Color', [1 0.5 0], 'LineWidth', 2)
    errorbar(Lag*AcquisitionTime/1000, MSD_2(:,3), MSD_2(:,2), '-s', 'Color', [0 0.4 1], 'LineWidth', 2)
    
    axis square
    ax.FontSize = FontSize;
    box on
    xlabel('Time(s)')
    ylabel('MSD (um²)')
    
    Title = sprintf('log(D_1) = %.2f um²/s -- log(D_2) = %.2f um²/s', D_mean(1,1), D_mean(1,2));
    title(Title);
    legend('D1 average', 'D2 average', 'D1 median', 'D2 median', 'Location', 'northwest');
    
    t = Lag*AcquisitionTime/1000;
    T = cat(2,t', MSD_1(:,1), MSD_1(:,3), MSD_1(:,2),...
        MSD_2(:,1), MSD_2(:,3), MSD_2(:,2));
    T = array2table(T, 'VariableNames', {'Time_s', 'MSD1_Average', 'MSD1_Median', 'MSD1_Error_Bar', ...
        'MSD2_Average', 'MSD2_Median', 'MSD2_Error_Bar'});
    writetable(T, 'Saved_MSD.txt');
    
else
    
    MSD = MSD_FIT;
    
    errorbar(Lag*AcquisitionTime/1000, MSD(:,1),  MSD(:,2), '-o', 'Color', [0 0.4 1])
    hold on
    errorbar(Lag*AcquisitionTime/1000, MSD(:,3),  MSD(:,2), '-s', 'Color', [0 0.4 1], 'LineWidth', 2)
    [fitobject,~] = fit(Lag(1:4)'*AcquisitionTime/1000, MSD(1:4,1), 'poly1');
    
    Fit = fitobject.p1 * Lag'*AcquisitionTime/1000 + fitobject.p2;
    plot(Lag*AcquisitionTime/1000, Fit, '-r')
    axis square
    axis tight
    box on
    ax.FontSize = FontSize;
    xlabel('Time(s)')
    ylabel('MSD (um²/s)')
    legend('Mean values', 'Median values', 'Location', 'northwest')
    
    t = Lag*AcquisitionTime/1000;
    T = cat(2,t', MSD(:,1), MSD(:,3), MSD(:,2));
    T = array2table(T, 'VariableNames', {'Time_s', 'MSD1_Average', 'MSD1_Median', 'MSD1_Error_Bar'});
    writetable(T, 'Saved_MSD.txt');
    
end

saveas(ax, 'MSD_Curves.png');

%% If the option was selected, the trajectories are plotted on a single graph
%% ==========================================================================

if Plot_Traj
    
    hPlot = figure;
    set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
    scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
    set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen
    
    if isfield(h, 'AvIm')
        imagesc(h.AvIm)
        colormap('gray')
    else
        fig = gcf;
        fig.Color = [0.6 0.6 0.6];
    end
    
    axis image
    axis off
    legend off
    title('')
    
    Color = jet;
    ntraj_color = ceil(NTraj_ROI/size(Color,1));
    
    for ntraj = 1 : NTraj_ROI
        
        if isfield(h, 'AvIm')
            X = Reconstructed_Traj_ROI{ntraj}(2,:)/PixelSize;
            Y = Reconstructed_Traj_ROI{ntraj}(3,:)/PixelSize;
        else
            X = Reconstructed_Traj_ROI{ntraj}(2,:);
            Y = Reconstructed_Traj_ROI{ntraj}(3,:);
        end
        line(Y, X, 'Color', Color(ceil(ntraj/ntraj_color),:),'LineWidth',1)
        
    end
    
    if isfield(h, 'AvIm')
        
        ScaleBar = [5, 5, 1/PixelSize, 1];
        rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
        
    else
        AxisLimits = axis;
        Box = [AxisLimits(1), AxisLimits(3), AxisLimits(2)-AxisLimits(1), AxisLimits(4)-AxisLimits(3)];
        rectangle('Position', Box, 'EdgeColor', [0 0 0], 'LineWidth', 2)
        
        ScaleBar = [AxisLimits(1)+1, AxisLimits(3)+1, 1, 0.1];
        rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
    end
    
    saveas(hPlot, 'Trajectories.png');
end

%% Save the parameters in the file called MTT_sptPALM_analysis.mat
%% ===============================================================

Results.MaxBlinks = MaxBlinks;
Results.MinNPoint = MinNPoint;
Results.MinTrajLength_MSDCalculation = MinTrajLength_MSDCalculation;
Results.MaxDisplayTime = MaxDisplayTime;
Results.p = p;
Results.MaxStepLength = MaxStepLength;
Results.MinNPointMSD = MinNPointMSD;

Results.Reconstructed_Traj_Filtered = Reconstructed_Traj_Filtered;

Results.Reconstructed_Traj_ROI = Reconstructed_Traj_ROI;
h.Reconstructed_Traj_ROI = Reconstructed_Traj_ROI;

Results.MSD_all = MSD_all;
Results.MSD_weight = MSD_weight;
Results.Reconstructed_Traj_MSD = Reconstructed_Traj_MSD;
h.Reconstructed_Traj_MSD = Reconstructed_Traj_MSD;

Results.Dapp = Dapp;
Results.Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD_accepted;
h.Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD_accepted;

if NbrGaussianFit==1
    h.Reconstructed_Traj_Diff = varargout{2};
    h.FittedDiffDistribution = varargout{3};
    h.DiffDistribution = varargout{4};
    Results.Reconstructed_Traj_Diff = varargout{2};
    Results.FittedDiffDistribution = varargout{3};
    Results.DiffDistribution = varargout{4};
else
    h.Reconstructed_Traj_DiffPop1 = varargout{2};
    h.Reconstructed_Traj_DiffPop2 = varargout{3};
    h.Fraction = varargout{4};
    h.FittedDiffDistribution = varargout{5};
    h.DiffDistribution = varargout{6};
    Results.Reconstructed_Traj_DiffPop1 = varargout{2};
    Results.Reconstructed_Traj_DiffPop2 = varargout{3};
    Results.Fraction = varargout{4};
    Results.FittedDiffDistribution = varargout{5};
    Results.DiffDistribution = varargout{6};
end

h.D_mean = D_mean;
h.MSD_FIT = MSD_FIT;
h.NbrGaussianFit = NbrGaussianFit;

Results.NbrGaussianFit = NbrGaussianFit;
Results.D_mean = D_mean;
Results.MSD_FIT = MSD_FIT;
Results.Density = Density;

%% The figure are saved and the parameters as well in a txt file
%% ==============================================================

if isunix
    fileID = fopen(strcat(h.DirectoryName, '/', 'Parameters_analysis.txt'),'w');
else
    fileID = fopen(strcat(h.DirectoryName, '\', 'Parameters_analysis.txt'),'w');
end

fprintf(fileID, 'The following parameters are only the one used for the MATLAB analysis');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');

fprintf(fileID, '\r\n %s', 'Maximum number of blinks');
fprintf(fileID, '\n\n\n %4.2f\n', MaxBlinks);
fprintf(fileID, '\r\n %s', 'Minimum length trajectory');
fprintf(fileID, '\n\n\n %4.2f\n', MinNPoint);
fprintf(fileID, '\r\n %s', 'Minimum length trajectory for the MSD calculation');
fprintf(fileID, '\n\n\n %4.2f\n', MinTrajLength_MSDCalculation);
fprintf(fileID, '\r\n %s', 'Maximum display time for the MSD');
fprintf(fileID, '\n\n\n %4.2f\n', MaxDisplayTime);
fprintf(fileID, '\r\n %s', 'Number of points used to calculate the apparent D');
fprintf(fileID, '\n\n\n %4.2f\n', p);
fprintf(fileID, '\r\n %s', 'Minimum number of distance values used to calculate each point on the MSD curve');
fprintf(fileID, '\n\n\n %4.2f\n', MinNPointMSD);
fprintf(fileID, '\r\n %s', 'Maximum accepted distance separating two consecutive detections');
fprintf(fileID, '\n\n\n %4.2f\n', MaxStepLength);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Number of movies analyzed ');
fprintf(fileID, '\n\n\n %s\n', get(h.NMovies, 'String'));
fprintf(fileID, '\r\n %s', 'Acquisition time (ms)');
fprintf(fileID, '\n\n\n %4.2f\n', AcquisitionTime);
fprintf(fileID, '\r\n %s', 'Pixel size (um)');
fprintf(fileID, '\n\n\n %4.2f\n', PixelSize);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Initial number of trajectories ');
fprintf(fileID, '\n\n\n %4.2f\n', size(Reconstructed_Traj,1));
fprintf(fileID, '\r\n %s', 'Number of trajectories after filtering for the blinking and duration');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_Filter);
fprintf(fileID, '\r\n %s', 'Number of trajectories after ROI selection');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_ROI);
fprintf(fileID, '\r\n %s', 'Number of trajectories validated for MSD calculation');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_MSD);
fprintf(fileID, '\r\n %s', 'Number of trajectories validated for the Dapp calculation');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_Diff);
fprintf(fileID, '\r\n %s', 'Density of tracks detected (/um²)');
fprintf(fileID, '\n\n\n %4.2f\n', Density);

fclose(fileID);
set(h.SaveForTesseler, 'Enable', 'on')
set(h.PlotPreviousAnalysis, 'Enable', 'on')
set(h.DataTypePlot, 'Enable', 'on')
