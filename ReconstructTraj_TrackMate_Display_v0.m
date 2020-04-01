%*****************************
%
% ReconstructTraj_TrackMate_Display_v0.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the TrackMate data associated to a movie previously loaded 
% in order to visualize the trajectories and the detection of single 
% molecule events
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function h = ReconstructTraj_TrackMate_Display_v0(FileToAnalyse, h)

PixelSize = str2double(get(h.PixelSize, 'String')); % in µm

%% Load the TrackMate file selected 
%% ================================

[m, NTrajectory] = importTrackMateTracks(FileToAnalyse, 1); % Load the results of the TrackMate analysis file

%% For each TrackMate file, the trajectories are analyzed and saved in "Reconstructed_Traj"
%% ==================================================================================

Reconstructed_Traj = cell(NTrajectory,1);
Frame_Traj = zeros(NTrajectory,2);

clc
fprintf('\n Formating the trajectories ...     ')

for ntraj = 1 : NTrajectory
    
    fprintf('\b\b\b\b%03i%%', round(100*ntraj/ntraj))
    
    T = m{ntraj}(:,1);
    X = PixelSize*m{ntraj}(:,2);
    Y = PixelSize*m{ntraj}(:,3);
    D = sqrt( (X(2:end) - X(1:end-1)).^2 + (Y(2:end) - Y(1:end-1)).^2 );
    
    Reconstructed_Traj{ntraj} = cat(1,T',X',Y',cat(2, 0, D'));
    Frame_Traj(ntraj,:) = [T(1), T(end)];
end

fprintf('\n')

h.Reconstructed_Traj_MovieDisplay = Reconstructed_Traj;
h.Frame_Traj_MovieDisplay = Frame_Traj;