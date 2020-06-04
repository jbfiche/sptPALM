%*****************************
%
% ReconstructTraj_Display_v5.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the MTT data associated to a movie previously loaded 
% in order to visualize the trajectories and the detection of single 
% molecule events
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function h = ReconstructTraj_Display_v5(FileToAnalyse, h)

PixelSize = str2double(get(h.PixelSize, 'String')); % in µm

%% For each MTT file, the trajectories are analyzed and saved in "Reconstructed_Traj"
%% ==================================================================================

Reconstructed_Traj = {};
Frame_Traj = [];

m = matfile(FileToAnalyse); % Load the results of the MTT analysis
Xmatrix = m.Xmatrix;
Ymatrix = m.Ymatrix;
NTrajectory = size(Xmatrix, 1);

clc 
disp('Loading MTT tracking data ...')

% The arrays Xmatrix and Ymatrix are analyzed in order to detect the
% presence of NaN. They usually means that for a specific frame, there
% was non detection and the columns is filled with NaN and 0... not
% clear why there is a mixture of both though!
% In order to avoid the creation of artefactual trajectories (very long
% ones), all the 0 of the columns are converted to NaN
% ----------------------------------------------------

for ncol = 1 : size(Xmatrix,2)
    
    NaN_X = sum(isnan(Xmatrix(:,ncol)));
    NaN_Y = sum(isnan(Ymatrix(:,ncol)));
    
    if NaN_X>0 || NaN_Y>0
        
        Xmatrix(:,ncol) = NaN;
        Ymatrix(:,ncol) = NaN;
    end
end

for ntraj = 1 : NTrajectory
    
    X = PixelSize*Xmatrix(ntraj,:);
    Y = PixelSize*Ymatrix(ntraj,:);
    
    % The way MTT works,can be described as follows:
    %
    % 1- When a particle is detected at the frame #n, the value assigned
    % to Xmatrix and Ymatrix for ALL the frames before n (1 : n-1) is
    % set to 0.
    % 2- When a particle is lost, the soft keep adding the same
    % positions to the output .txt file.
    
    % As a result, when calculating the distance between successive
    % detections, most of the trajectories have a D equal to 0. Therefore,
    % below, each trajectory is defined by the points comprised between
    % the very first detection (the first non-zero values of X and Y)
    % and the last non-zero detection (that is, the last non-zero value
    % of D)
    
    % In parallel, the step length between two consecutive detections
    % is also calculated for each trajectory and saved in
    % "SingleStep_Length".
    % --------------------
    
    Idx_FirstDetection = find(X>0 & Y>0,1);
    D = sqrt((X(Idx_FirstDetection+1:end)-X(Idx_FirstDetection:end-1)).^2 + (Y(Idx_FirstDetection+1:end)-Y(Idx_FirstDetection:end-1)).^2);
    Idx = find(D>0);
    
    if size(Idx,2)>=2 && ~isempty(Idx_FirstDetection)
        
        D = D(1:Idx(end));
        Reconstructed_Traj = cat(1, Reconstructed_Traj,...
            cat(1, Idx_FirstDetection:Idx_FirstDetection+Idx(end), X(Idx_FirstDetection:Idx_FirstDetection+Idx(end)), Y(Idx_FirstDetection:Idx_FirstDetection+Idx(end)), cat(2, 0, D)));
        Frame_Traj = cat(1, Frame_Traj, [Idx_FirstDetection, Idx_FirstDetection+Idx(end)]);

    end
end

h.Reconstructed_Traj_MovieDisplay = Reconstructed_Traj;
h.Frame_Traj_MovieDisplay = Frame_Traj;

disp('Data loaded.')
disp(' ')