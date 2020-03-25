%*****************************
%
% Simulated_Trajectory_analysis_v3.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function was written in order to analyze the simulated data
% and check that the theoretical results match the parameters used for the
% simulation (ratio of population, diffusion coefficient).
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function fileID = Simulated_Trajectory_analysis_v3(hPlot, FolderName, MaxBlink, MinNPoint, MinTrajLength_MSDCalculation, AcquisitionTime, MaxDisplayTime, p, MinNPointMSD, Reconstructed_Traj, PixelSize)

%% Parameter for the analysis
%% ==========================

FontSize = 15;

%% Filter the trajectories according to the parameters selected.
%% =============================================================

[Reconstructed_Traj_ROI, NTraj_ROI] = Filter_Trajectories_simulation_v1(Reconstructed_Traj, MaxBlink, MinTrajLength_MSDCalculation, MinNPoint);

fprintf('\r\n');
fprintf('------------------------------------------------------\r\n');
fprintf('%i trajectories have been selected for the calculation', NTraj_ROI)
fprintf('\r\n');

%% Analyze the trajectories and calculate the MSD
%% for each trajectory (added as the fifth row)
%% ================================================

[MSD_all,MSD_weight,Reconstructed_Traj_MSD] = MSD_calculation(Reconstructed_Traj_ROI,NTraj_ROI,MinNPointMSD,p);

NTraj_MSD = size(MSD_all,1);
fprintf('%i trajectories have been selected for the calculation (after MSD calculation)', NTraj_MSD)
fprintf('\r\n');

%% The apparent diffusion coefficient is calculated for each single trajectory
%% using the first "p" points of the MSD curves
%% ============================================

% For each point, the apparent coefficient diffusion is calculated
% using either a linear fit on the first p points of the average of the MSD
% for a lagtime of 1 frame
% -------------------------

[MSD_all_Method1,Reconstructed_Traj_MSD_accepted_Method1,Dapp_Method1] = Diff_calculation(1,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD);
[MSD_all_Method2,Reconstructed_Traj_MSD_accepted_Method2,Dapp_Method2] = Diff_calculation(2,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD);

NTraj_Diff_Method1 = size(Dapp_Method1,1);
NTraj_Diff_Method2 = size(Dapp_Method2,1);

fprintf('\n');
fprintf('%i trajectories have been selected for the calculation (method based on the fit) \r', NTraj_Diff_Method1)
fprintf('%i trajectories have been selected for the calculation (method based on the weighted average) \r', NTraj_Diff_Method2)
fprintf('\n');

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
        
    LogDapp = log10(Dapp);  
    [NbrGaussianFit, D_mean, varargout] = Simulation_FitGaussianDistribution_v1(LogDapp, MSD_all, FontSize, hPlot, round(MaxDisplayTime*1000/AcquisitionTime), Reconstructed_Traj_MSD_accepted);
    
    if Method == 1
        MSD_FIT_Method1 = varargout{1};
        D_mean_Method1 = D_mean;
        saveas(hPlot, 'Diffusion_distribution_Fit_Method.png');
    else
        MSD_FIT_Method2 = varargout{1};
        D_mean_Method2 = D_mean;
        saveas(hPlot, 'Diffusion_distribution_Weighted_Average_Method.png');
    end
end

%% Plot the MSD curve for the method based on the fits
%% ==================================================s

for Method = 1 : 2
    
    if Method == 1
        MSD_FIT = MSD_FIT_Method1;
        D_mean = D_mean_Method1;
    else
        MSD_FIT = MSD_FIT_Method2;
        D_mean = D_mean_Method2;
    end
    
    Lag = 1 : 1 : size(MSD_FIT,1);
    
    figure(hPlot)
    ax = gca;
    hold off
    cla
    
    if NbrGaussianFit == 2
        
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
        ylabel('MSD (um²)')
        
        Title = sprintf('log(D_1) = %.2f um²/s -- log(D_2) = %.2f um²/s', D_mean(1,1), D_mean(1,2));
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
        ylabel('MSD (um²/s)')
        Title = sprintf('log(D) = %.2f um²/s', D_mean);
        title(Title);
        legend('Mean values', 'Median values', 'Location', 'northwest')
        
    end
    
    if Method == 1
        saveas(hPlot, 'MSD_Curves_Fit_Method.png');
    else
        saveas(hPlot, 'MSD_Curves_Weighted_Average_Method.png');
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

fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n %s', 'Acquisition time (ms)');
fprintf(fileID, '\n\n\n %4.2f\n', AcquisitionTime);
fprintf(fileID, '\r\n %s', 'Pixel size (um)');
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
