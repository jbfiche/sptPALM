function h = MTT_Trajectory_analysis_v4(h)

%% Parameter for the analysis
%% ==========================

MaxBlinks = str2num(get(h.MaxBlinks, 'String')); % Maximum number of frames two successive detections can be separated of
MinNPoint = str2num(get(h.MinNumberPoints, 'String')); % Minimum number of frames for each trajectory
MinTrajLength_MSDCalculation = str2num(get(h.MinTrajLength, 'String')); % Minium number of frames for the trajectories used for the calculation of the MSD
AcquisitionTime = str2num(get(h.AcquisitionTime, 'String')); % in ms
PixelSize = str2num(get(h.PixelSize, 'String')); %  in µm
MaxDisplayTime = str2num(get(h.MaxDisplayTime, 'String')); % The MSD curve will be displayed only from zero to this value in s
p = str2num(get(h.NumberPointsMSDFit, 'String')); % Use the first "p" points of the MSD to estimate the Dapp
MinNPointMSD = str2num(get(h.MinimumNumberPointsMSD, 'String')); % Minimum number of points used to calculate each values of the MSD
MaxStepLength = str2num(get(h.MaxStepLength, 'String')); % Minimum number of points used to calculate each values of the MSD

FontSize = h.FontSize;
DiffCalculationMethod = get(h.DiffusionCalculationMethod, 'Value');
Reconstructed_Traj = h.Reconstructed_Traj;

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

%% Open the figure where the results will be saved
%% ===============================================

hPlot = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen

%% Filter the trajectories according to the parameters selected
%% ============================================================

[Reconstructed_Traj_Filtered, NTraj_Filter] = Filter_Trajectories(Reconstructed_Traj, MaxBlinks, MinTrajLength_MSDCalculation, MinNPoint, MaxStepLength);
set(h.NTrajectoriesFiltered, 'String', num2str(NTraj_Filter)); % Display the # of trajectories selected after applying the filters
h.Reconstructed_Traj_Filtered = Reconstructed_Traj_Filtered;

%% Look whether we want to define an ROI or not. If a ROI is defined, return
%% only the trajectories detected within the ROI
%% =============================================

[h, Area, Reconstructed_Traj_ROI, NTraj_ROI] = Select_Trajectories_ROI(h, Results, hPlot, Reconstructed_Traj_Filtered, PixelSize);
set(h.NTrajectoriesROI, 'String', num2str(NTraj_ROI));  % Display the # of trajectories detected within the ROI
h.Reconstructed_Traj_ROI = Reconstructed_Traj_ROI;

%% Analyze the trajectories and calculate the MSD for each trajectory
%% (added as the fifth row)
%% ========================

hwaitbar = waitbar(0,'Calculating the MSD');
MSD_all = {};
MSD_weight = {};
Idx_TrajAcceptedMSD = [];
NLongest = 0;

for  ntraj = 1 : NTraj_ROI
    
    MSD = [];
    Weight = [];
    
    waitbar(ntraj/NTraj_ROI);
    Traj = Reconstructed_Traj_ROI{ntraj};
    LagMax = (Traj(1,end) - Traj(1,1)); % Calculate the maximum lag time for this trajectory
    
    % Calculate the MSD. The lagtime goes from "1" to "LagTime-(MinNPointMSD-1)"
    % since we want, for each value of the MSD, an average over at least 
    % "MinNPointMSD" different values.
    % For each lagtime values, the MSD is kept only if there are at least
    % "MinNPointMSD" distances used for the calculation. Else, the value is
    % discarted.
    % --------- 
    
    for lag = 1 : LagMax-(MinNPointMSD-1)
        
        D_all = [];
%         D = 0;
        n = 1;
        nMSD = 0;
        MaxLagPoint = Traj(1, end) - lag; % Return the value of time after which there is no point left in the trajectory that could be separated by a lagtime equals to "lag"
        MaxNPoint = find(Traj(1,:)>MaxLagPoint,1); % Return a list of points that could not be used for the MSD calculation for lagtime = lag;
        
        while n < MaxNPoint && MaxNPoint >= MinNPointMSD
            
            Idx = find(Traj(1,:)==Traj(1,n)+lag, 1); % Check that two events have been detected with the rigth lag time
            if ~isempty(Idx)
                
                Xi = Traj(2,n);
                Xj = Traj(2,Idx);
                Yi = Traj(3,n);
                Yj = Traj(3,Idx);
                d = (Xi - Xj)^2 + (Yi - Yj)^2;
%                 D = D + d;
                D_all = cat(1, D_all, d);
                nMSD = nMSD+1;
            end
            
            n = n+1;
        end
        
        % If the number of segments for the calculation of the MSD is lower
        % than MinNPointMSD, the calculation is stopped
        % ---------------------------------------------
        
        if nMSD > MinNPointMSD
            MSD(lag) = mean(D_all);
            Weight(lag) = std(D_all);
%             Weight(lag) = nMSD;
        else
            break
        end
    end
    
    % Several trajectories (with blinks) can still gives MSD that have less
    % that p points, which can result in bug/error for the calculation of 
    % the apparent diffusion coefficient. In order to avoid this issue, a 
    % last check is performed here and all the MSD array with less than p 
    % points are discarded.
    % ---------------------
    
    if size(MSD,2)>=p
        MSD_all = cat(1, MSD_all, MSD);
        MSD_weight = cat(1, MSD_weight, Weight);
        Idx_TrajAcceptedMSD = cat(1, Idx_TrajAcceptedMSD, ntraj);
    end
    
%     % For the purpose of specific analysis, the longest tracks are
%     % searched below
%     % --------------
%     
%     if  length(MSD)>45
%         
%         NLongest = NLongest+1;
%         
%         figure(5)
%         ax = gca;
%         
%         T = 1 : 1 : length(MSD);
%         errorbar(T*0.02,MSD,Weight, '-ob')
%         axis square
%         
%         figure(6)
%         X = Traj(2,:);
%         Y = Traj(3,:);
%         plot(X,Y,'-ok','MarkerSize', 10, 'LineWidth', 0.5, 'MarkerFaceColor', 'k')
%         axis equal
%         
%         disp('')
%     end
    
end

close(hwaitbar);
Reconstructed_Traj_MSD = Reconstructed_Traj_ROI(Idx_TrajAcceptedMSD);
NTraj_MSD = size(MSD_all,1);
set(h.NTrajectoriesMSD, 'String', num2str(NTraj_MSD)); % Display the # of trajectories accepted for MSD calculation
h.Reconstructed_Traj_MSD = Reconstructed_Traj_MSD;

%% Return the track density and display it on the GUI
%% ==================================================

if isfield(h, 'AvIm') && exist('Area', 'var')~=0
    Density = round(1000*NTraj_ROI/(sum(Area)*PixelSize^2))/1000;
    set(h.TrackDensity, 'String', num2str(Density)); % Display the density of tracks per µm²
    h.Density = Density;
elseif isfield(h, 'AvIm') && exist('Area', 'var')==0
    Density = round(1000*NTraj_ROI/(size(h.AvIm,1)*size(h.AvIm,2)*PixelSize^2))/1000;
    set(h.TrackDensity, 'String', num2str(Density)); % Display the density of tracks per µm²
    h.Density = Density;
else
    Density = NaN;
    set(h.TrackDensity, 'String', ''); % Display the density of tracks per µm²
end

%% The apparent diffusion coefficient is calculated for each single trajectory
%% using the first "p" points of the MSD curves
%% ============================================

% For each point, the apparent coefficient diffusion is calculated
% using either a linear fit on the first p points of the average of the MSD
% for a lagtime of 1 frame
% -------------------------

Idx_MSD_accepted = [];
hwaitbar = waitbar(0, 'Calculating the apparent diffusion coefficient ...');

Lag = 1 : 1 : p;
Lag = Lag*AcquisitionTime/1000;
Dapp = [];

for nMSD = 1 : size(MSD_all,1)
    
    waitbar(nMSD /size(MSD_all,1));
    MSD = MSD_all{nMSD}(1:p);
    Weight = MSD_weight{nMSD}(1:p);
    Nnan = sum(isnan(MSD));
    Nzeros = sum(MSD(:)>0);
    
    if DiffCalculationMethod == 2 && size(MSD,2)>=p
        
        % Calculation of Dapp using the fit
        % ---------------------------------
              
        if Nzeros == p && Nnan == 0
            
            %                 figure(1)
            %                 cla
            %                 hold on
            %                 plot(Lag', MSD_all(nMSD, 1 : Minp+(p-1))', 'o')
            %                 plot(0, 0, 'o')
            
            [fitobject,gof] = fit(Lag', MSD', 'Poly1', 'Weights', Weight', 'Upper', [Inf, min(MSD)]);
            if fitobject.p1>0
                Dapp(end+1,1) = fitobject.p1/4; % Return the coefficient of diffusion
                %                     Dapp(end,2) = fitobject.p2; % Return the dynamic localization uncertainty (4 s^2)
                %                     D(end+1,1) = fitobject.p1; % Return the coefficient of diffusion
                %                     D(end,2) = sqrt(fitobject.p2); % Return the dynamic localization uncertainty (4 s^2)
                %                         fitobject = ezfit(Lag', MSD_all(nMSD, 1 : Minp+(p-1) )', 'poly1'); % Replaced the function fit by ezfit to avoid licence problem
                %                         if fitobject.m(2)>0
                %                             D(end+1,1) = fitobject.m(2)/4; % Return the coefficient of diffusion
                %                             D(end,2) = sqrt(fitobject.m(1)/4); % Return the dynamic localization uncertainty (4 s^2)
                Idx_MSD_accepted(end+1,1) = nMSD; % Save the index of the MSD values that are used for the calculation of the apparent diffusion coefficient
            end
        end
        
    else
        
        % Calculation of Dapp using the average
        % -------------------------------------

        if Nzeros == p && Nnan == 0
                Dapp(end+1,1) = sum(MSD.*Weight./(4*Lag))/sum(Weight);
                Idx_MSD_accepted(end+1,1) = nMSD;
        end
        
    end
end

MSD_all = MSD_all(Idx_MSD_accepted); % Keep only the points in MSD_all that were used to calculate the distribution of apparent diffusion coefficient
Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD(Idx_MSD_accepted); % Keep only the trajectories that associated to the accepted MSD 

NTraj_Diff = size(Dapp,1);
set(h.NTrajectoriesDiff, 'String', num2str(NTraj_Diff)); % Display the # of trajectories selected for D calculation
h.Dapp = Dapp;
close(hwaitbar)

%% Plot the distribution of apparent diffusion coefficient
%% =======================================================

figure(hPlot)
ax = gca;
hold off
cla

LogDapp = log10(Dapp);
Bin = min(LogDapp(:,1)) : (max(LogDapp(:,1))-min(LogDapp(:,1)))/50 : max(LogDapp(:,1));
[N,Bin] = hist(LogDapp(:,1), Bin);
N = 100*N/sum(N);

bar(Bin, N, 'FaceColor', [0 0.4 1])
ax.FontSize = FontSize;
axis square
box on
xlabel('Log of apparent diffusion coefficient (µm²/s)')
ylabel('Fraction of molecule (%)')

T = cat(2, LogDapp(:,1), Dapp(:,1));
T = array2table(T, 'VariableNames', {'logD', 'D'});
writetable(T, 'Saved_Diffusion_Coeff.txt');

[NbrLorentzianFit, D_mean, varargout] = FitLorentzianDistribution_v2(N, Bin, LogDapp, MSD_all, FontSize, hPlot, round(MaxDisplayTime*1000/AcquisitionTime), Reconstructed_Traj_MSD_accepted, DiffCalculationMethod);
MSD_FIT = varargout{1};

%% Plot the MSD curve
%% ==================

Lag = 1 : 1 : size(MSD_FIT,1);

figure(hPlot)
ax = gca;
hold off
cla

if NbrLorentzianFit == 2
    
    MSD_1 = MSD_FIT(:,1:3);
    MSD_2 = MSD_FIT(:,4:6);
    
    errorbar(Lag*AcquisitionTime/1000, MSD_1(:,1), MSD_1(:,2), '-o', 'Color', [1 0.5 0])
    hold on
    errorbar(Lag*AcquisitionTime/1000, MSD_2(:,1), MSD_2(:,2), '-o', 'Color', [0 0.4 1])
%     [fitobject1,gof1] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_1(Idx(1:4),1), 'poly1');
%     [fitobject2,gof2] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_2(Idx(1:4),1), 'poly1');
    
    errorbar(Lag*AcquisitionTime/1000, MSD_1(:,3), MSD_1(:,2), '-s', 'Color', [1 0.5 0], 'LineWidth', 2)
    errorbar(Lag*AcquisitionTime/1000, MSD_2(:,3), MSD_2(:,2), '-s', 'Color', [0 0.4 1], 'LineWidth', 2)
%     [fitobject1,gof1] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_1(Idx(1:4),1), 'poly1');
%     [fitobject2,gof2] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_2(Idx(1:4),1), 'poly1');
    
    axis square
    ax.FontSize = FontSize;
    box on
    xlabel('Time(s)')
    ylabel('MSD (µm²)')
    
    Title = sprintf('log(D_1) = %.2f µm²/s -- log(D_2) = %.2f µm²/s', D_mean(1,1), D_mean(1,2));
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
    ylabel('MSD (µm²/s)')
    legend('Mean values', 'Median values', 'Location', 'northwest')
    
    t = Lag*AcquisitionTime/1000;
    T = cat(2,t', MSD(:,1), MSD(:,3), MSD(:,2));
    T = array2table(T, 'VariableNames', {'Time_s', 'MSD1_Average', 'MSD1_Median', 'MSD1_Error_Bar'});
    writetable(T, 'Saved_MSD.txt');
    
end

export_fig(hPlot, 'MSD_Curves.png');
% export_fig(hPlot, 'MSD_Curves.pdf', '-pdf');

%% The trajectories are plotted on a single graph
%% ==============================================

figure(hPlot)
hold off
cla

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

export_fig(hPlot, 'Trajectories.png');
% export_fig(hPlot, 'Trajectories.pdf', '-pdf');

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

if NbrLorentzianFit==1
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
h.NbrLorentzianFit = NbrLorentzianFit;

Results.NbrLorentzianFit = NbrLorentzianFit;
Results.D_mean = D_mean;
Results.MSD_FIT = MSD_FIT;
Results.Density = Density;

%% The figure are saved and the parameters as well in a txt file
%% ==============================================================

fileID = fopen(strcat(h.DirectoryName, '\', 'Parameters_analysis.txt'),'w');
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
fprintf(fileID, '\r\n %s', 'Pixel size (µm)');
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
fprintf(fileID, '\r\n %s', 'Density of tracks detected (/µm²)');
fprintf(fileID, '\n\n\n %4.2f\n', Density);

fclose(fileID);
set(h.SaveForTesseler, 'Enable', 'on')
set(h.PlotPreviousAnalysis, 'Enable', 'on')
set(h.DataTypePlot, 'Enable', 'on')
close(hPlot)
