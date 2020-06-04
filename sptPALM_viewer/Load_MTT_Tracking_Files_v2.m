%*****************************
%
% Load_MTT_Tracking_Files_v2.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Load the mat files obtained as output from the MTT software.  
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function [h, Launch_Analysis, FileToAnalyse] = Load_MTT_Tracking_Files_v2(h)

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

        %% Check the .mat files are related to the MTT analysis or something else.
        %% Only the MTT .mat files are kept. The others are discarded.
        %% ===========================================================
        
        hwb = waitbar(0, 'Looking for the MTT files ...');
        
        for Nfiles = numel(FileToAnalyse) : -1 : 1
            waitbar((numel(FileToAnalyse)-Nfiles)/numel(FileToAnalyse));
            m = matfile(FileToAnalyse{Nfiles}); % Load the results of the MTT analysis
            variables = whos(m);
            first_variable = variables(1).name;

            if ~isequal(first_variable, 'Xmatrix')
                FileToAnalyse(Nfiles) = [];
            end
        end
        close(hwb)
        
        %% Check whether there are files to analyse or not
        %% ===============================================
        
        if ~isempty(FileToAnalyse)
            h.FileToAnalyse = FileToAnalyse;
            h.DirectoryName = DirectoryName;
            set(h.NMovies, 'String', num2str(size(FileToAnalyse, 1))); % Display on the GUI front pannel the number of files analyzed
            set(h.AnalyseTrajectories, 'Enable', 'on')
            set(h.DiffusionCalculationMethod, 'Enable', 'on')
        else
            warndlg('No MTT files were found. No analysis could be performed.')
        end
end

%% Save the list of folders in a mat file
%% ======================================

cd(DirectoryName)
