function [h, Repeat_Analysis] = Load_MTT_Tracking_Files_v2(h)

%% Check whether an analysis was already run for this folder and
%% if the results would be used again
%% ==================================

DirectoryName = uigetdir;
cd(DirectoryName)
set(h.FolderPath_Text, 'String', DirectoryName); % Display on the GUI the folder path

PreviousAnalysis = ~isempty(dir('*MTT_sptPALM_analysis*.mat'));

if PreviousAnalysis
    Repeat_Analysis = questdlg('A previous analysis was found. Do you want to used it?', 'Repeat?', 'Yes', 'No', 'Yes');
else
    Repeat_Analysis = 'No';
end

switch Repeat_Analysis
    
    case 'Yes'
        set(h.PlotPreviousAnalysis, 'Enable', 'on')
        h = sptPALM_initialize(h, 'LoadPreviousData');
        
        if isempty(h.FileToAnalyse)
            warndlg('No MTT files were found. No analysis could be performed.')
        else
            set(h.AnalyseTrajectories, 'Enable', 'on')
            set(h.DiffusionCalculationMethod, 'Enable', 'on')
            
            if isfield(h, 'NbrLorentzianFit') && isfield(h, 'FittedDiffDistribution') && isfield(h, 'DiffDistribution') && isfield(h, 'D_mean') && isfield(h, 'MSD_FIT')
                set(h.SaveForTesseler, 'Enable', 'on')
                set(h.PlotPreviousAnalysis, 'Enable', 'on')
                set(h.DataTypePlot, 'Enable', 'on')
            end
        end
    
    case 'No'
        set(h.PlotPreviousAnalysis, 'Enable', 'off')
        
        %% Analyse the containt of the folder and look for all the .mat files
        %% ==================================================================
               
        FileToAnalyse = LookForDirectories_spt(DirectoryName);

        %% Check the .mat files are related to the MTT analysis or something else.
        %% Only the MTT .mat files are kept. The others are discarded.
        %% ===========================================================
        
        hwb = waitbar(0, 'Looking for the MTT files ...');
        
        for Nfiles = numel(FileToAnalyse) : -1 : 1
            waitbar((numel(FileToAnalyse)-Nfiles)/numel(FileToAnalyse));
            m = matfile(FileToAnalyse{Nfiles}); % Load the results of the MTT analysis
            try
                size(m.Xmatrix, 1);
            catch
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
