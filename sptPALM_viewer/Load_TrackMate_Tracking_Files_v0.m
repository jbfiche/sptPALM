%*****************************
%
% Load_TrackMate_Tracking_Files_v0.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the excel files obtained as output from the TrackMate software.  
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function [h, Launch_Analysis] = Load_TrackMate_Tracking_Files_v0(h)

%% Check whether an analysis was already run for this folder and
%% if the results would be used again
%% ==================================

DirectoryName = uigetdir;
cd(DirectoryName)
set(h.FolderPath_Text, 'String', DirectoryName); % Display on the GUI the folder path

PreviousAnalysis = ~isempty(dir(h.ResultsFileName));

if PreviousAnalysis
    Launch_Analysis = questdlg('A previous analysis with the same name was found. Do you want to proceed anyway (previous data will be lost)?',...
        'Previous data found?', 'Proceed', 'Load previous analysis', 'Cancel', 'Proceed');
else
    Launch_Analysis = 'Proceed';
end

switch Launch_Analysis
    
    case 'Cancel'

        return
        
    case 'Load previous analysis'
        
        set(h.PlotPreviousAnalysis, 'Enable', 'on')
        h = sptPALM_initialize(h, 'LoadPreviousData');
        
        set(h.AnalyseTrajectories, 'Enable', 'on')
        set(h.DiffusionCalculationMethod, 'Enable', 'on')
        
        if isfield(h, 'FittedDiffDistribution') && isfield(h, 'DiffDistribution') && isfield(h, 'D_mean') && isfield(h, 'MSD_FIT')
            set(h.SaveForTesseler, 'Enable', 'on')
            set(h.PlotPreviousAnalysis, 'Enable', 'on')
            set(h.DataTypePlot, 'Enable', 'on')
        end
    
    case 'Proceed'
        
        set(h.PlotPreviousAnalysis, 'Enable', 'off')
        
        %% Analyse the containt of the folder and look for all the .mat files
        %% ==================================================================
               
        FileToAnalyse = LookForDirectories_spt(DirectoryName, h.MTT_FileName.String);
        NFiles = size(FileToAnalyse,1);

        %% Check the .xls files and import the trajectories
        %% ================================================
        
        Tracks = cell(NFiles,1);
        Validated_files = zeros(NFiles,1);
        Total_tracks = 0;
        
        for nfile = NFiles : -1 : 1
            try
                n_file = NFiles - nfile + 1;
                [m, ntracks] = importTrackMateTracks(FileToAnalyse{nfile}, n_file); % Load the results of the TrackMate analysis file
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
end

%% Save the list of folders in a mat file
%% ======================================

cd(DirectoryName)
