%*****************************
%
% Load_TrackMate_Tracking_Files_batch_v1.m
%
% ****************************
%
% JB Fiche
% June, 2023
% Last update : 2023/06/27
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the data files obtained as output from the newest version
% of the TrackMate software (using the batch tool). This version is
% dedicated to the batch version of the analysis.
% -------------------------------------------------------------------------
% Specific:
% -------------------------------------------------------------------------
% To fix:
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function h = Load_TrackMate_Tracking_Files_batch_v1(h, DirectoryName)

%% Check whether an analysis was already run for this folder and
%% if the results would be used again
%% ==================================
cd(DirectoryName)
set(h.FolderPath_Text, 'String', DirectoryName); % Display on the GUI the folder path

set(h.PlotPreviousAnalysis, 'Enable', 'off')

%% Analyse the containt of the folder and look for all the .csv files.
%% For the newest version of TrackMate *_spots.csv is returning the
%% positions of all the detected spots and to which track it is associated.
%% =======================================================================

FileToAnalyse = LookForDirectories_spt(DirectoryName, h.Data_FileName.String);
NFiles = size(FileToAnalyse, 1);

%% Check the .csv files and import the trajectories
%% ================================================

Tracks = cell(NFiles,1);
Validated_files = zeros(NFiles,1);
Total_tracks = 0;

for nfile = NFiles : -1 : 1
    try
        n_file = NFiles - nfile + 1;
        [m, ntracks] = importTrackMateTracks_v1(FileToAnalyse{nfile}, n_file); % Load the results of the TrackMate analysis file
        Total_tracks = Total_tracks + ntracks;
        
        if size(m,1)>1
            Tracks{nfile} = m;
            Validated_files(nfile) = 1;
        end
    catch error
        fprintf('\n%s\n', FileToAnalyse{nfile})
        fprintf('%s\n', error.identifier)
        fprintf('%s\n', error.message)
    end
end

Tracks = Tracks(Validated_files==1);
FileToAnalyse = FileToAnalyse(Validated_files==1);

%% Check whether there are files to analyse or not
%% ===============================================

if Total_tracks>0
    h.FileToAnalyse = FileToAnalyse;
    h.DirectoryName = DirectoryName;
    h.TrackMate = Tracks;
    h.Total_tracks = Total_tracks;
    
    set(h.NMovies, 'String', num2str(size(FileToAnalyse, 1))); % Display on the GUI front pannel the number of files analyzed
    set(h.AnalyseTrajectories, 'Enable', 'on')
    set(h.DiffusionCalculationMethod, 'Enable', 'on')
else
    warndlg('No TrackMate files were found. No analysis could be performed.')
end
