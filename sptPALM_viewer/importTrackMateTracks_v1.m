%*****************************
%
% importTrackMateTracks_v1.m
%
% ****************************
%
% JB Fiche
% July, 2022
% Last update : 2022/07/19
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the .csv files output by trackmate batch  
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function [tracks, nTracks] = importTrackMateTracks_v1(file, filenumber)

%% load the spots file. All the headers are removed and the following columns are
%% kept : Spot ID / Track ID / X / Y / Frame
%% =========================================

fprintf('\n Importing trajectories from file # %i ...     ', filenumber);

spots = readcell(file, 'Delimiter', ',',  'HeaderLines', 4);
spots = spots(:,[2,3,5,6,9]);
tracks = cell2mat(spots);
nTracks = max(tracks(:,2));


