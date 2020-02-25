function h = sptPALM_initialize(h, Action)

switch Action
    
    case 'Reset_all'
        
        % Initialize all the callback to their default values and states
        % ---------------------------------------------------------------
        
        set(h.FolderPath_Text, 'String', 'Directory name')
        
        set(h.AnalyseTrajectories, 'Enable', 'off');
        set(h.DiffusionCalculationMethod, 'Enable', 'off');
        set(h.SaveForTesseler, 'Enable', 'off');
        set(h.PlotPreviousAnalysis, 'Enable', 'off');
        set(h.LoadTraj, 'Enable', 'off');
        set(h.DataTypePlot, 'Enable', 'off');
        set(h.CreateROI, 'Enable', 'off');
        set(h.DeleteROI, 'Enable', 'off');
        
        set(h.MaxStepLength, 'String','');
        
        set(h.NMovies, 'String','');
        set(h.NTrajectories, 'String','');
        set(h.NTrajectoriesFiltered, 'String','');
        set(h.NTrajectoriesROI, 'String','');
        set(h.NTrajectoriesMSD, 'String','');
        set(h.NTrajectoriesDiff, 'String','');
        set(h.TrackDensity, 'String','');
        
        set(h.Slider_SelectFrame, 'Enable', 'off');
        set(h.Slider_UpperContrast, 'Enable', 'off');
        set(h.Slider_LowerContrast, 'Enable', 'off');
        
        set(h.Edit_SelectFrame, 'Enable', 'off');
        set(h.Edit_UpperContrast, 'Enable', 'off');
        set(h.Edit_LowerContrast, 'Enable', 'off');
        set(h.MinTrackLength, 'Enable', 'off');
        
        set(h.Slider_UpperContrast,'Min', 0, 'Max', 1, 'Value', 0.01,...
            'SliderStep', [1/2000, 1/2000]);
        set(h.Slider_LowerContrast,'Min', 0, 'Max', 1, 'Value', 0,...
            'SliderStep', [1/2000, 1/2000]);
        set(h.Edit_UpperContrast, 'String','0.01');
        set(h.Edit_LowerContrast, 'String','0');
        
        if isfield(h, 'AvIm')
            h = rmfield(h, 'AvIm');
        end
        
        if isfield(h, 'ROI')
            h = rmfield(h, 'ROI');
        end
        
    case 'Reset_h'
        
        if isfield(h, 'Reconstructed_Traj_ROI')
            h = rmfield(h, 'Reconstructed_Traj_ROI');
        end
        
        if isfield(h, 'Reconstructed_Traj_Filtered')
            h = rmfield(h, 'Reconstructed_Traj_Filtered');
        end
        
        if isfield(h, 'Reconstructed_Traj_MSD')
            h = rmfield(h, 'Reconstructed_Traj_MSD');
        end
        
        if isfield(h, 'Reconstructed_Traj_MSD_accepted')
            h = rmfield(h, 'Reconstructed_Traj_MSD_accepted');
        end
        
        if isfield(h, 'Dapp')
            h = rmfield(h, 'Dapp');
        end
        
        if isfield(h, 'Reconstructed_Traj_Diff')
            h = rmfield(h, 'Reconstructed_Traj_Diff');
        end
        
        if isfield(h, 'FittedDiffDistribution')
            h = rmfield(h, 'FittedDiffDistribution');
        end
        
        if isfield(h, 'DiffDistribution')
            h = rmfield(h, 'DiffDistribution');
        end
        
        if isfield(h, 'D_mean')
            h = rmfield(h, 'D_mean');
        end
        
        if isfield(h, 'MSD_FIT')
            h = rmfield(h, 'MSD_FIT');
        end
        
        if isfield(h, 'NbrLorentzianFit')
            h = rmfield(h, 'NbrLorentzianFit');
        end
        
        if isfield(h, 'Reconstructed_Traj_DiffPop1') && isfield(h, 'Reconstructed_Traj_DiffPop2')
            h = rmfield(h, 'Reconstructed_Traj_DiffPop1');
            h = rmfield(h, 'Reconstructed_Traj_DiffPop2');
        end
        
        if isfield(h, 'Fraction')
            h = rmfield(h, 'Fraction');
        end
        
        if isfield(h, 'Density')
            h = rmfield(h, 'Density');
        end
         
    case 'LoadPreviousData'
        
        NFiles = dir(h.ResultsFileName);
        Results = load(NFiles(1).name);
        h.ResultsFileName = NFiles(1).name;

%         if size(NFiles,1) == 1
%             Results = load(NFiles(1).name);
%             h.ResultsFileName = NFiles(1).name;
%         elseif size(NFiles,1) > 1
%             hwarn = warndlg('More than one file were found. Manually select the file you want to use:');
%             uiwait(hwarn);
%             [FileName,~,~] = uigetfile('*MTT_sptPALM_analysis*.mat');
%             Results = load(FileName);
%             h.ResultsFileName = FileName;
%         end

        h.DirectoryName = cd;
        
        if isfield(Results, 'AcquisitionTime') && isfield(Results, 'PixelSize')
            set(h.AcquisitionTime, 'String', num2str(Results.AcquisitionTime));
            set(h.PixelSize, 'String', num2str(Results.PixelSize));
        end
        
        if isfield(Results, 'MaxStepLength')
            set(h.MaxStepLength, 'String', num2str(Results.MaxStepLength));
            set(h.MaxBlinks, 'String', num2str(Results.MaxBlinks));
            set(h.MinTrajLength, 'String', num2str(Results.MinTrajLength_MSDCalculation));
            set(h.MinNumberPoints, 'String', num2str(Results.MinNPoint));
            set(h.MinimumNumberPointsMSD, 'String', num2str(Results.MinNPointMSD));
            set(h.NumberPointsMSDFit, 'String', num2str(Results.p));
            set(h.MaxDisplayTime, 'String', num2str(Results.MaxDisplayTime));
        end
        
        if isfield(Results, 'FileToAnalyse') && isfield(Results, 'DirectoryName')
            set(h.NMovies, 'String', num2str(size(Results.FileToAnalyse,1)));
            h.FileToAnalyse = Results.FileToAnalyse;
            h.DirectoryName = Results.DirectoryName;
        end
        
        if isfield(Results, 'Reconstructed_Traj')
            set(h.NTrajectories, 'String', num2str(size(Results.Reconstructed_Traj,1)));
            h.Reconstructed_Traj = Results.Reconstructed_Traj;
            h.SingleStep_Length = Results.SingleStep_Length;
            h.Length_Traj = Results.Length_Traj;
        end
        
        if isfield(Results, 'Reconstructed_Traj_Filtered')
            set(h.NTrajectoriesFiltered, 'String', num2str(size(Results.Reconstructed_Traj_Filtered,1)));
            h.Reconstructed_Traj_Filtered = Results.Reconstructed_Traj_Filtered;
        end
        
        if isfield(Results, 'Reconstructed_Traj_ROI')
            set(h.NTrajectoriesROI, 'String', num2str(size(Results.Reconstructed_Traj_ROI,1)));
            h.Reconstructed_Traj_ROI = Results.Reconstructed_Traj_ROI;
        end
        
        if isfield(Results, 'Reconstructed_Traj_MSD')
            set(h.NTrajectoriesMSD, 'String', num2str(size(Results.Reconstructed_Traj_MSD,1)));
            h.Reconstructed_Traj_MSD = Results.Reconstructed_Traj_MSD;
        end
        
        if isfield(Results, 'Reconstructed_Traj_MSD_accepted')
            h.Reconstructed_Traj_MSD_accepted = Results.Reconstructed_Traj_MSD_accepted;
        end
        
        if isfield(Results, 'Dapp')
            set(h.NTrajectoriesDiff, 'String', num2str(size(Results.Dapp,1)));
            h.Dapp = Results.Dapp;
        end
        
        if isfield(Results, 'Density')
            if ~isnan(Results.Density)
                set(h.TrackDensity, 'String', num2str(Results.Density));
                h.Density = Results.Density;
            else
                set(h.TrackDensity, 'String', '');
            end
        end
        
        if isfield(Results, 'Reconstructed_Traj_DiffPop1') && isfield(Results, 'Reconstructed_Traj_DiffPop2')
            h.Reconstructed_Traj_DiffPop1 = Results.Reconstructed_Traj_DiffPop1;
            h.Reconstructed_Traj_DiffPop2 = Results.Reconstructed_Traj_DiffPop2;
            h.Fraction = Results.Fraction;
        end
        
        if isfield(Results, 'Reconstructed_Traj_Diff')
            h.Reconstructed_Traj_Diff = Results.Reconstructed_Traj_Diff;
        end
        
        if isfield(Results, 'FittedDiffDistribution')
            h.FittedDiffDistribution = Results.FittedDiffDistribution;
        end
        
        if isfield(Results, 'DiffDistribution')
            h.DiffDistribution = Results.DiffDistribution;
        end
        
        if isfield(Results, 'D_mean')
            h.D_mean = Results.D_mean;
        end
        
        if isfield(Results, 'MSD_FIT')
            h.MSD_FIT = Results.MSD_FIT;
        end
        
        if isfield(Results, 'NbrLorentzianFit')
            h.NbrLorentzianFit = Results.NbrLorentzianFit;
        end
        
        if isfield(Results, 'Fraction')
            h.Fraction = Results.Fraction;
        end
        
        if isfield(Results, 'ROI')
            h.ROI = Results.ROI;
        end
        
        if isfield(Results, 'AvIm')
            h.AvIm = Results.AvIm;
        end
end