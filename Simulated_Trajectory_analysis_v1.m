function fileID = Simulated_Trajectory_analysis_v1(hPlot, FolderName, MaxBlink, MinNPoint, MinTrajLength_MSDCalculation, AcquisitionTime, PixelSize, MaxDisplayTime, p, MinNPointMSD, MaxStepLength, DiffCalculationMethod, Reconstructed_Traj)

%% Parameter for the analysis
%% ==========================

clc
FontSize = 15;

%% Filter the trajectories according to the parameters selected
%% ============================================================

[Reconstructed_Traj_ROI, NTraj_ROI] = Filter_Simulated_Trajectories(Reconstructed_Traj, MaxBlink, MinTrajLength_MSDCalculation, MinNPoint, MaxStepLength);

fprintf('\r\n');
fprintf('\r\n');
fprintf('------------------------------------------------------\r\n');
fprintf('%i trajectories have been selected for the calculation', NTraj_ROI)
fprintf('\r\n');

%% Analyze the trajectories and calculate the MSD
%% for each trajectory (added as the fifth row)
%% ================================================

hwaitbar = waitbar(0,'Calculating the MSD');
MSD_all = {};
MSD_weight = {};
Idx_TrajAcceptedMSD = [];

for  ntraj = 1 : NTraj_ROI
    
    MSD = [];
    Weight = [];
    
    waitbar(ntraj/NTraj_ROI);
    Traj = transpose(Reconstructed_Traj_ROI{ntraj});
    LagMax = (Traj(1, end) - Traj(1,1)); % Calculate the maximum lag time for this trajectory
    
    % Calculate the MSD. The lagtime goes from "1" to "LagTime-p" since we
    % want, for each value of the MSD, an average over at least "p" different
    % values.
    % For each lagtime values, the MSD is kept only if there are at least
    % "MinNPointMSD" distances used for the calculation. Else, the value is
    % discarted.
    % --------- 
    
    for lag = 1 : LagMax-(MinNPointMSD-1)
        
        D = 0;
        n = 1;
        nMSD = 0;
        MaxLagPoint = Traj(1, end) - lag; % Return the value of time after which there is no point left in the trajectory that could be separated by a lagtime equals to "lag"
        MaxNPoint = find(Traj(1,:)>MaxLagPoint,1); % Return a list of points that could not be used for the MSD calculation for lagtime = lag;
        
        while n < MaxNPoint && (MaxNPoint-1) >= MinNPointMSD
            
            Idx = find(Traj(1,:)==Traj(1,n)+lag, 1); % Check that two events have been detected with the rigth lag time
            if ~isempty(Idx)
                
                Xi = Traj(2,n)*PixelSize;
                Xj = Traj(2,Idx)*PixelSize;
                Yi = Traj(3,n)*PixelSize;
                Yj = Traj(3,Idx)*PixelSize;
                d = (Xi - Xj)^2 + (Yi - Yj)^2;
                D = D + d;
                nMSD = nMSD+1;
            end
            
            n = n+1;
        end
        
        if nMSD > MinNPointMSD
            MSD(lag) = 1/nMSD*D;
            Weight(lag) = nMSD;
        end
    end
    
    % Despite the selection, several trajectories (with blinks) can
    % still gives MSD that have less that p points. In order to avoid
    % this issue, a last check is performed here and all the MSD array
    % with less than p points are discarded.
    % -------------------------------------
    
    if size(MSD,2)>p
        MSD_all = cat(1, MSD_all, MSD);
        MSD_weight = cat(1, MSD_weight, Weight);
        Idx_TrajAcceptedMSD = cat(1, Idx_TrajAcceptedMSD, ntraj);
    end
end

% Save the results
% ----------------

Reconstructed_Traj_MSD = Reconstructed_Traj_ROI(Idx_TrajAcceptedMSD);
Results.MSD_all = MSD_all;
Results.MSD_weight = MSD_weight;
Results.Reconstructed_Traj_MSD = Reconstructed_Traj_MSD;
close(hwaitbar);

NTraj_MSD = size(MSD_all,1);

fprintf('%i trajectories have been selected for the calculation', NTraj_MSD)
fprintf('\r\n');

%% The apparent diffusion coefficient is calculated for each single trajectory
%% using the first "p" points of the MSD curves
%% ============================================

% For each point, the apparent coefficient diffusion is calculated
% using either a linear fit on the first p points of the average of the MSD
% for a lagtime of 1 frame
% -------------------------

hwaitbar = waitbar(0, 'Calculating the apparent diffusion coefficient ...');

Lag = 1 : 1 : p;
Lag = Lag*AcquisitionTime/1000;
Dapp_Method1 = [];
Dapp_Method2 = [];
Idx_MSD_accepted_Method1 = [];
Idx_MSD_accepted_Method2 = [];

for nMSD = 1 : size(MSD_all,1)
    
    waitbar(nMSD /size(MSD_all,1));
    MSD = MSD_all{nMSD}(1:p);
    Weight = MSD_weight{nMSD}(1:p);
    Nnan = sum(isnan(MSD));
    Nzeros = sum(MSD(:)>0);
    
    if size(MSD,2)>=p
        
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
                Dapp_Method1(end+1,1) = fitobject.p1/4; % Return the coefficient of diffusion
                %                     Dapp(end,2) = fitobject.p2; % Return the dynamic localization uncertainty (4 s^2)
                %                     D(end+1,1) = fitobject.p1; % Return the coefficient of diffusion
                %                     D(end,2) = sqrt(fitobject.p2); % Return the dynamic localization uncertainty (4 s^2)
                %                         fitobject = ezfit(Lag', MSD_all(nMSD, 1 : Minp+(p-1) )', 'poly1'); % Replaced the function fit by ezfit to avoid licence problem
                %                         if fitobject.m(2)>0
                %                             D(end+1,1) = fitobject.m(2)/4; % Return the coefficient of diffusion
                %                             D(end,2) = sqrt(fitobject.m(1)/4); % Return the dynamic localization uncertainty (4 s^2)
                Idx_MSD_accepted_Method1(end+1,1) = nMSD; % Save the index of the MSD values that are used for the calculation of the apparent diffusion coefficient
            end
        end
        
        
        % Calculation of Dapp using the average
        % -------------------------------------

        if Nzeros == p && Nnan == 0
                Dapp_Method2(end+1,1) = sum(MSD.*Weight./(4*Lag))/sum(Weight);
                Idx_MSD_accepted_Method2(end+1,1) = nMSD;
        end
        
    end
end

MSD_all_Method1 = MSD_all(Idx_MSD_accepted_Method1);
MSD_all_Method2 = MSD_all(Idx_MSD_accepted_Method2); % Keep only the points in MSD_all that were used to calculate the distribution of apparent diffusion coefficient
Reconstructed_Traj_MSD_accepted_Method1 = Reconstructed_Traj_MSD(Idx_MSD_accepted_Method1);
Reconstructed_Traj_MSD_accepted_Method2 = Reconstructed_Traj_MSD(Idx_MSD_accepted_Method2); % Keep only the trajectories that associated to the accepted MSD 

NTraj_Diff_Method1 = size(Dapp_Method1,1);
NTraj_Diff_Method2 = size(Dapp_Method2,1);

fprintf('%i trajectories have been selected for the calculation (method based on the fit)', NTraj_Diff_Method1)
fprintf('%i trajectories have been selected for the calculation (method based on the weighted average)', NTraj_Diff_Method2)
fprintf('\r\n');
close(hwaitbar)

%% Plot the distribution of apparent diffusion coefficient for two methods
%% =======================================================================

for Method = 1 : 2
    
    if Method == 1
        Dapp = Dapp_Method1;
        MSD_all = MSD_all_Method1;
        Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD_accepted_Method1;
    else
        Dapp = Dapp_Method2;
        MSD_all = MSD_all_Method2;
        Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD_accepted_Method2;
    end
    
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
    
    [NbrLorentzianFit, D_mean, varargout] = Simulation_FitLorentzianDistribution_v1(N, Bin, LogDapp, MSD_all, FontSize, hPlot, round(MaxDisplayTime*1000/AcquisitionTime), Reconstructed_Traj_MSD_accepted);
    
    if Method == 1
        MSD_FIT_Method1 = varargout{1};
        export_fig(hPlot, 'Diffusion_distribution_Fit_Method.png');
    else
        MSD_FIT_Method2 = varargout{1};
        export_fig(hPlot, 'Diffusion_distribution_Weighted_Average_Method.png');
    end
end

%% Plot the MSD curve for the method based on the fits
%% ==================================================s

for Method = 1 : 2
    
    if Method == 1
        MSD_FIT = MSD_FIT_Method1;
    else
        MSD_FIT = MSD_FIT_Method2;
    end
    
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
        Title = sprintf('log(D) = %.2f µm²/s', D_mean);
        title(Title);
        legend('Mean values', 'Median values', 'Location', 'northwest')
        
    end
    
    if Method == 1
        export_fig(hPlot, 'MSD_Curves_Fit_Method.png');
    else
        export_fig(hPlot, 'MSD_Curves_Weighted_Average_Method.png');
    end
end

%% The figure are saved and the parameters as well in a txt file
%% ==============================================================

fileID = fopen(strcat(FolderName, '\', 'Parameters_analysis.txt'),'w');
fprintf(fileID, 'The following parameters were used for the simulation');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');

fprintf(fileID, '\r\n %s', 'Maximum number of blinks');
fprintf(fileID, '\n\n\n %4.2f\n', MaxBlink);
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
fprintf(fileID, '\r\n %s', 'Acquisition time (ms)');
fprintf(fileID, '\n\n\n %4.2f\n', AcquisitionTime);
fprintf(fileID, '\r\n %s', 'Pixel size (µm)');
fprintf(fileID, '\n\n\n %4.2f\n', PixelSize);

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Initial number of trajectories ');
fprintf(fileID, '\n\n\n %4.2f\n', size(Reconstructed_Traj,1));
fprintf(fileID, '\r\n %s', 'Number of trajectories after applying the filters (blinking ...)');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_ROI);
fprintf(fileID, '\r\n %s', 'Number of trajectories validated for MSD calculation');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_MSD);
fprintf(fileID, '\r\n %s', 'Number of trajectories validated for the Dapp calculation (method based on the fit)');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_Diff_Method1);
fprintf(fileID, '\r\n %s', 'Number of trajectories validated for the Dapp calculation (method based on the weighted average)');
fprintf(fileID, '\n\n\n %4.2f\n', NTraj_Diff_Method2);
