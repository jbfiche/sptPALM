%*****************************
%
% sptPALM_viewer.m
%
% ****************************
%
% JB Fiche
% Feb, 2015
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Analyze tracking data obtained from MTT software (Serge et al.
% 2008, Nat. Methods) or TrackMate (Tinevez et al. 2018, Methods). 
% The program can be used to sort the tracks, calculate the MSD and 
% diffusion coefficients Dinst, plot the trajectories and the
% distribution of Dinst. Data can also be simulated. 
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function sptPALM_viewer(varargin)


%% Initialize the GUI
%% ==================

clc
defaultBackground = get(0,'defaultUicontrolBackgroundColor');
scrsz = get(0,'ScreenSize');

if scrsz(4)<1000
    ScreenOption = 'small';
else
    ScreenOption = 'large';
end

% Calculate the dimension of the GUI window according to the disposition of
% the displayed windows.
% ----------------------

Side_Border = 15;
Height_Border = 5;

switch ScreenOption
    case 'small'
        ControlPanel_Width = 580;
        ControlPanel_Height = 800;
        DisplayMovie_Width = 680;
        DisplayMovie_Height = 720;
        Axis_Size = 550;
    case 'large'
        ControlPanel_Width = 600;
        ControlPanel_Height = 950;
        DisplayMovie_Width = 750;
        DisplayMovie_Height = 750;
        Axis_Size = 600;
end

h.sptPALM_DisplayMovie = figure('Visible','on',...
    'Units','pixels',...
    'Position',[1,1,DisplayMovie_Width,DisplayMovie_Height],...
    'MenuBar','none',...
    'Name','sptPALM_Display_Panel',...
    'NumberTitle','off',...
    'Toolbar','none',...
    'Resize','off',...
    'Color',defaultBackground,...
    'DefaultUicontrolUnits','Pixels',...
    'DefaultUipanelUnits','Pixels',...
    'CloseRequestFcn',@sptPALM_DisplayMovie_closefunction);

h.sptPALM_ControlPanel = figure('Visible','on',...
    'Units','pixels',...
    'Position',[1,1,ControlPanel_Width,ControlPanel_Height],...
    'MenuBar','none',...
    'Name','sptPALM_Control_Panel',...
    'NumberTitle','off',...
    'Toolbar','none',...
    'Resize','off',...
    'Color',defaultBackground,...
    'CloseRequestFcn',@sptPALM_viewer_closefunction);

% Construct the components
% ------------------------

h = sptPALM_components(h, ScreenOption, Side_Border, Height_Border, ControlPanel_Width, ControlPanel_Height);

% Define the axis on the display figure
% -------------------------------------

Corner_x = (DisplayMovie_Width - Axis_Size)/(2*DisplayMovie_Width);
Corner_y = (DisplayMovie_Height - Axis_Size)/(2*DisplayMovie_Height);
Lx = Axis_Size/DisplayMovie_Width;
Ly = Axis_Size/DisplayMovie_Height;

h.MainAxes = axes('Parent', h.sptPALM_DisplayMovie, 'box', 'on', 'XTick', [], 'YTick', [], ...
    'Position', [Corner_x,Corner_y,Lx,Ly]);

% Define default parameters for the MSD analysis (this parameters were
% previously availabe on the GUI but for the sake of clarity, they were
% removed and defined during software initialization).
% ----------------------------------------------------

h.MinNumberPoints = 0.75; % Minimum percentage of points to consider a trajectory complete
h.MinimumNumberPointsMSD = 3; % Minimum number of points needed to calculate a MSD value

% Define the default parameters for the simulation
% ------------------------------------------------

SimulationParameters.NProteinsTot_Init = 10000; % Define the total pool of proteins that can be activated and imaged
SimulationParameters.MeanProbActivation = 0.5; % The probability is defined as the mean number of protein activated per frame
SimulationParameters.MeanPhotoBleachingTime = 1; % in s
SimulationParameters.MeanPhotons = 200; % Average number of photons emitted
SimulationParameters.Ton = 0.3; % in s
SimulationParameters.Toff1 = 0.3; % in s
SimulationParameters.Toff2 = 4; % in s

SimulationParameters.MaxBlink = 3; % in frames
SimulationParameters.MinTrajLength = 7; %in frames

SimulationParameters.ImageSize = 252; % px
SimulationParameters.PixelSize = 0.102; % in um
SimulationParameters.AcquisitionTime = 0.02; % in s
SimulationParameters.ExposureTime = 0.02; % in s
SimulationParameters.Gain = 200; % electronic gain of the camera
SimulationParameters.Offset = 98; % intensity offset for the image
SimulationParameters.QY = 0.9; % parameter lambda for the Poisson distribution describing the noise of the EMCCD camera
SimulationParameters.CCDsensitivity = 63; % parameter describing the normal readout noise for each pixel and linked to the photon detection
SimulationParameters.ReadoutNoise = 72; % parameters defining the half-size of the window where the 2D gaussian will be calculated
SimulationParameters.LimitResolution = 0.15; % in um
SimulationParameters.GaussStampSize = 7; % in px

h.SimulationParameters = SimulationParameters;

% Set the Callbacks
% -----------------

set(h.Saving_file_name, 'callback', @UpdateSavingFileName);
set(h.TrackingSoftware, 'callback', @SelectTrackingSoftware);
set(h.LoadMTT, 'callback', @LoadMTT);
set(h.Load_Previous_analysis, 'callback', @Load_Previous_analysis);
set(h.MinTrajLength, 'callback', @CheckMinTrajSize);
set(h.AnalyseTrajectories, 'callback', @AnalyseTrajectories);
set(h.SaveForTesseler, 'callback', @SaveForTesseler);
set(h.PlotPreviousAnalysis, 'callback', @PlotPreviousAnalysis);
set(h.LoadMovie, 'callback', @LoadMovie);
set(h.LoadTraj,'callback', @LoadTraj);
set(h.CreateROI,'callback', @CreateROI);
set(h.DeleteROI,'callback', @DeleteROI);
set(h.Slider_SelectFrame,'callback', @Slider_SelectFrame);
set(h.Edit_SelectFrame,'callback', @Edit_SelectFrame);
set(h.Slider_UpperContrast,'callback', @Slider_UpperContrast);
set(h.Edit_UpperContrast,'callback', @Edit_UpperContrast);
set(h.Slider_LowerContrast,'callback', @Slider_LowerContrast);
set(h.Edit_LowerContrast,'callback', @Edit_LowerContrast);
set(h.SaveMovie,'callback', @Save_Movie);
set(h.Simulation_EmissionParameters,'callback', @Change_EmissionParameters);
set(h.Simulation_AquisitionParameters,'callback', @Change_AquisitionParameters);
set(h.LaunchSimulation,'callback', @LaunchSimulation);

% Define a copy of the handle h (called h_backup) in case we need to reset
% the handle completely. Finally initialize the control panel.
% ------------------------------------------------------------

h.FontSize = 10;
h.ResultsFileName = 'MTT_sptPALM_analysis.mat';

global h_backup
global h_backup_analysis
h_backup = h;

h = sptPALM_initialize(h, 'Reset_all');


%% Select the tracking software used for analyzing the raw data
%% ============================================================

    function SelectTrackingSoftware(~,~)

        Soft = get(h.TrackingSoftware, 'Value');
        switch Soft
            case 1
                set(h.MTT_FileName, 'String', '*.mat');
                set(h.Saving_file_name, 'String', 'MTT_sptPALM_analysis.mat');
                set(h.LoadMTT, 'String', 'Load MTT files');
                set(h.LoadMTT, 'callback', @LoadMTT);
                
            case 2
                set(h.MTT_FileName, 'String', '*.xml');
                set(h.Saving_file_name, 'String', 'TrackMate_sptPALM_analysis.mat');
                set(h.LoadMTT, 'String', 'Load TrackMate files');
                set(h.LoadMTT, 'callback', @LoadTrackMate);
        end
    end

%% Change the saving file name
%% ===========================

    function UpdateSavingFileName(~,~)
        h.ResultsFileName = h.Saving_file_name.String;
        h_backup.ResultsFileName = h.Saving_file_name.String;
        h_backup_analysis.ResultsFileName = h.Saving_file_name.String;
    end

%% Load MTT tracking files
%% ===================

    function LoadMTT(~,~)
        
        clc
        h = sptPALM_initialize(h, 'Reset_all');
        h.ResultsFileName = h.Saving_file_name.String;
        
        [h, Repeat_Analysis] = Load_MTT_Tracking_Files_v2(h);
        SelectTrackingSoftware
        
        if isequal(Repeat_Analysis, 'Proceed')
            clear_display_axis
            h = ReconstructTraj_v6(h);
            h_backup_analysis = h;
        end
    end

%% Load TrackMate files
%% ====================

    function LoadTrackMate(~,~)
        
        clc
        h = sptPALM_initialize(h, 'Reset_all');
        h.ResultsFileName = h.Saving_file_name.String;
        
        [h, Repeat_Analysis] = Load_TrackMate_Tracking_Files_v0(h);
        
        if isequal(Repeat_Analysis, 'Proceed')
            clear_display_axis
            h = ReconstructTraj_TrackMate_v6(h);
            h_backup_analysis = h;
        end
    end
        
%% Load previous analysis
%% ======================

    function Load_Previous_analysis(~,~)
        
        clc
        clear_display_axis
        h = sptPALM_initialize(h, 'Reset_all');
        
        [FileName, PathName] = uigetfile('*.mat');
        
        if FileName ~= 0
            
            cd(PathName);
            h.ResultsFileName = FileName;
            set(h.Saving_file_name, 'String', FileName);
            
            set(h.PlotPreviousAnalysis, 'Enable', 'on')
            h = sptPALM_initialize(h, 'LoadPreviousData');
            
            SelectTrackingSoftware
            
            set(h.AnalyseTrajectories, 'Enable', 'on')
            set(h.DiffusionCalculationMethod, 'Enable', 'on')
            
            if isfield(h, 'FittedDiffDistribution') && isfield(h, 'DiffDistribution') && isfield(h, 'D_mean') && isfield(h, 'MSD_FIT')
                set(h.SaveForTesseler, 'Enable', 'on')
                set(h.PlotPreviousAnalysis, 'Enable', 'on')
                set(h.DataTypePlot, 'Enable', 'on')
            end
        end
    end

%% Check that the value chosen for the minimum trajectory length is
%% in agreement with the value of the parameter "p"
%% ================================================

    function CheckMinTrajSize(~,~)
        clc
        p = str2double(get(h.NumberPointsMSDFit, 'String'));
        MinTrajLength = str2double(get(h.MinTrajLength, 'String'));
        
        if MinTrajLength<p
            hwarn = warndlg('The minimum trajectory length cannot be lower than the number of points used for the calculation of Dapp!');
            uiwait(hwarn)
            set(h.MinTrajLength, 'String', num2str(p))
        end
    end

%% Analyze the trajectories
%% ========================

    function AnalyseTrajectories(~,~)
        clc
        clear_display_axis
        h = Trajectory_analysis_v4(h);
    end

%% Save the analysis for the Tesseler clustering analysis
%% ======================================================

    function SaveForTesseler(~,~)
        clc
        h = Save_For_Tesseler_sptPALM(h);
    end

%% Plot the previous analysis
%% ==========================

    function PlotPreviousAnalysis(~,~)
        clc
        clear_display_axis
        PlotPreviousAnalysis_v1(h)
    end

%% Load movie
%% ==========

    function LoadMovie(~,~)
        
        clc
        clear_display_axis
        set(h.sptPALM_DisplayMovie, 'Visible', 'on');
        
        [ImageName, ImageDirectory] = uigetfile('*.tif');
        h.MovieDisplayFullName = strcat(ImageDirectory, ImageName);
        ImInfo = imfinfo(h.MovieDisplayFullName);
        h.MovieDisplayNImages = length(ImInfo);
        h.ImSize = [ImInfo(1).Width, ImInfo(1).Height];
        
        set(h.Slider_SelectFrame, 'Min', 1)
        set(h.Slider_SelectFrame, 'Max', h.MovieDisplayNImages)
        set(h.Slider_SelectFrame, 'Value', 1)
        set(h.Slider_SelectFrame, 'SliderStep', [1/h.MovieDisplayNImages,1/h.MovieDisplayNImages])
        set(h.Edit_SelectFrame, 'String', '1')
        
        ChangeContrast_sptPALM_v2(h)
        
        if isfield(h, 'MovieDisplayROI')
            h = rmfield(h, 'MovieDisplayROI');
        end
        
        set(h.LoadTraj, 'Enable', 'on')
        set(h.MinTrackLength, 'Enable', 'on')
        set(h.MinTrackLength, 'String', '1')
        set(h.CreateROI, 'Enable', 'on')
        set(h.DeleteROI, 'Enable', 'on')
        set(h.Slider_SelectFrame, 'Enable', 'on')
        set(h.Edit_SelectFrame, 'Enable', 'on')
        set(h.Slider_UpperContrast, 'Enable', 'on')
        set(h.Edit_UpperContrast, 'Enable', 'on')
        set(h.Slider_LowerContrast, 'Enable', 'on')
        set(h.Edit_LowerContrast, 'Enable', 'on')
        set(h.SaveMovie, 'Enable', 'on');
        set(h.FirstFrame, 'Enable', 'on');
        set(h.LastFrame, 'Enable', 'on');
        set(h.MovieName, 'Enable', 'on');
    end

%% Load trajectories associated to the movie loaded
%% ================================================

    function LoadTraj(~,~)
        
        Soft = get(h.TrackingSoftware, 'Value');
        FileName_template = get(h.MTT_FileName, 'String');

        switch Soft
            case 1
                
                [MTTFileName, MTTFileDirectory] = uigetfile(FileName_template);
                FileToAnalyse = strcat(MTTFileDirectory, MTTFileName);
                h = ReconstructTraj_Display_v5(FileToAnalyse, h);
                
            case 2
                
                [TrackMateFileName, TrackMateFileDirectory] = uigetfile(FileName_template);
                FileToAnalyse = strcat(TrackMateFileDirectory, TrackMateFileName);
                h = ReconstructTraj_TrackMate_Display_v0(FileToAnalyse, h);
        end

        PlotTrajectories_sptPALM_v1(h)
    end

%% Create an ROI for the displayed movie
%% =====================================

    function CreateROI(~,~)
        
        axes(h.MainAxes)
        h.MovieDisplayROI = getrect(h.MainAxes);
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Delete the ROI create for the displayed movie
%% =============================================

    function DeleteROI(~,~)
        
        axes(h.MainAxes)
        if isfield(h, 'MovieDisplayROI')
            h = rmfield(h, 'MovieDisplayROI');
        end
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change frame using "Slider_SelectFrame"
%% =======================================

    function Slider_SelectFrame(~,~)
        
        clc
        N = round(get(h.Slider_SelectFrame, 'Value'));
        set(h.Edit_SelectFrame, 'String', num2str(N));
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change frame using "Edit_SelectFrame"
%% =======================================

    function Edit_SelectFrame(~,~)
        
        clc       
        N = round(str2double(get(h.Edit_SelectFrame, 'String')));
        if N<=0
            N = 1;
            set(h.Edit_SelectFrame, 'String', num2str(N));
        elseif N>h.MovieDisplayNImages
            N = h.MovieDisplayNImages;
            set(h.Edit_SelectFrame, 'String', num2str(N));
        end

        set(h.Slider_SelectFrame, 'Value', N);
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change the contrast using "Slider_UpperContrast"
%% ================================================

    function Slider_UpperContrast(~,~)
        
        clc
        Up = get(h.Slider_UpperContrast, 'Value');
        Low = get(h.Slider_LowerContrast, 'Value');
        
        if Up<=Low
            Up = Low+0.001;
            set(h.Slider_UpperContrast, 'Value', Up)
        end
        
        set(h.Edit_UpperContrast, 'String', num2str(Up))
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change the contrast using "Edit_UpperContrast"
%% ==============================================

    function Edit_UpperContrast(~,~)
        
        clc
        Up = str2double(get(h.Edit_UpperContrast, 'String'));
        Low = get(h.Slider_LowerContrast, 'Value');
        
        if Up<Low
            Up = Low+0.001;
            set(h.Edit_UpperContrast, 'String', num2str(Up))
        elseif Up>1
            Up = 1;
            set(h.Edit_UpperContrast, 'String', num2str(Up))
        end
        
        set(h.Slider_UpperContrast, 'Value', Up)
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change the contrast using "Slider_LowerContrast"
%% ================================================

    function Slider_LowerContrast(~,~)
       
        clc
        Up = get(h.Slider_UpperContrast, 'Value');
        Low = get(h.Slider_LowerContrast, 'Value');
        
        if Low>Up
            Low = Up-0.001;
            set(h.Slider_LowerContrast, 'Value', Low)
        end
        
        set(h.Edit_LowerContrast, 'String', num2str(Low))
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change the contrast using "Edit_LowerContrast"
%% ==============================================

    function Edit_LowerContrast(~,~)
        
        clc
        Up = get(h.Slider_UpperContrast, 'Value');
        Low = str2double(get(h.Edit_LowerContrast, 'String'));
        
        if Low>Up
            Low = Up-0.001;
            set(h.Edit_LowerContrast, 'String', num2str(Low))
        elseif Low<0
            Low = 0;
            set(h.Edit_LowerContrast, 'String', num2str(Low))
        end
        
        set(h.Slider_LowerContrast, 'Value', Low)
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Save the movie (as a avi file)
%% ===============================

    function Save_Movie(~,~)
        
        clc
        SaveDisplayedMovie(h)
    end

%% Change the parameters describing the emission properties of the emitters
%% ========================================================================

    function Change_EmissionParameters(~,~)
        
        clc
        % h.SimulationParameters = SimulationParameters;
        
        prompt = {'Total number of emitters available', ...
            'Average number of proteins activated per frame', ...
            'Mean Photo bleaching time (s)', ...
            'Mean number of emitted photons', ...
            'Average ON time (s)', ...
            'Average SHORT OFF time (s)', ...
            'Average LONG OFF time (s)', ...
            'Maximum number of frames allowed for the blinking', ...
            'Minimum trajectory length in frames'};
        dlg_title = 'Define emission parameters';
        num_lines= 1;
        default_answer = {num2str(h.SimulationParameters.NProteinsTot_Init), ...
            num2str(h.SimulationParameters.MeanProbActivation), ...
            num2str(h.SimulationParameters.MeanPhotoBleachingTime), ...
            num2str(h.SimulationParameters.MeanPhotons), ...
            num2str(h.SimulationParameters.Ton), ...
            num2str(h.SimulationParameters.Toff1), ...
            num2str(h.SimulationParameters.Toff2), ...
            num2str(h.SimulationParameters.MaxBlink), ...
            num2str(h.SimulationParameters.MinTrajLength)};
        NewEmissionParameters = inputdlg(prompt,dlg_title,num_lines,default_answer);
        
        if ~isempty(NewEmissionParameters)
            h.SimulationParameters.NProteinsTot_Init = str2double(NewEmissionParameters{1}); % Define the total pool of proteins that can be activated and imaged
            h.SimulationParameters.MeanProbActivation = str2double(NewEmissionParameters{2}); % The probability is defined as the mean number of protein activated per frame
            h.SimulationParameters.MeanPhotoBleachingTime = str2double(NewEmissionParameters{3}); % in s
            h.SimulationParameters.MeanPhotons = str2double(NewEmissionParameters{4}); % in s
            h.SimulationParameters.Ton = str2double(NewEmissionParameters{5}); % in s
            h.SimulationParameters.Toff1 = str2double(NewEmissionParameters{6}); % in s
            h.SimulationParameters.Toff2 = str2double(NewEmissionParameters{7}); % in s
            h.SimulationParameters.MaxBlink = str2double(NewEmissionParameters{8}); % in frames
            h.SimulationParameters.MinTrajLength = str2double(NewEmissionParameters{9}); %in frames
        end
    end

%% Change the parameters describing the aquisition parameters
%% ==========================================================

    function Change_AquisitionParameters(~,~)
        
        clc
        
        prompt = {'Size of the image in pixels (image is square by default)', ...
            'Pixel size(um)', ...
            'Aquisition time (s)', ...
            'Exposure time (s)', ...
            'Electronic gain', ...
            'Intensity offset', ...
            'Quantum yield', ...
            'CCD sensitivity (e/counts)', ...
            'Readout noise (e)', ...
            'Enter the half size of the window where the intensity of a single emitter is plot (in px)', ...
            'Enter the diffraction limit (in um and defined as the standard deviation)'};
        dlg_title = 'Define aquisition parameters';
        num_lines= 1;
        default_answer = {num2str(h.SimulationParameters.ImageSize), ...
            num2str(h.SimulationParameters.PixelSize), ...
            num2str(h.SimulationParameters.AcquisitionTime),...
            num2str(h.SimulationParameters.ExposureTime),...
            num2str(h.SimulationParameters.Gain), ...
            num2str(h.SimulationParameters.Offset), ...
            num2str(h.SimulationParameters.QY), ...
            num2str(h.SimulationParameters.CCDsensitivity), ...
            num2str(h.SimulationParameters.ReadoutNoise), ...
            num2str(h.SimulationParameters.GaussStampSize), ...
            num2str(h.SimulationParameters.LimitResolution)};
        NewAquisitionParameters = inputdlg(prompt,dlg_title,num_lines,default_answer);
        
        if ~isempty(NewAquisitionParameters)
            h.SimulationParameters.ImageSize = str2double(NewAquisitionParameters{1}); % px
            h.SimulationParameters.PixelSize = str2double(NewAquisitionParameters{2}); % um
            h.SimulationParameters.AcquisitionTime = str2double(NewAquisitionParameters{3}); % s
            h.SimulationParameters.ExposureTime = str2double(NewAquisitionParameters{4}); % s
            h.SimulationParameters.Gain = str2double(NewAquisitionParameters{5}); % um
            h.SimulationParameters.Offset = str2double(NewAquisitionParameters{6}); % intensity offset for the image
            h.SimulationParameters.QY = str2double(NewAquisitionParameters{7}); % quantum yield of the detector
            h.SimulationParameters.CCDsensitivity = str2double(NewAquisitionParameters{8}); % CCD sensitivity (conversion electrons to AD counts)
            h.SimulationParameters.ReadoutNoise = str2double(NewAquisitionParameters{9}); % readout noise of the detector
            h.SimulationParameters.GaussStampSize = str2double(NewAquisitionParameters{10}); % in um
            h.SimulationParameters.LimitResolution = str2double(NewAquisitionParameters{11}); % Mean intensity of single emitters
        end
    end

%% Launch the simulation of a sptPALM experiment
%% ===============================================

    function LaunchSimulation(~,~)
       clc
       clear_display_axis
       SimuDiff_v5(h)
    end

%% Function controlling the closing of the display window
%% =======================================================

    function sptPALM_DisplayMovie_closefunction(~,~)
        set(h.sptPALM_DisplayMovie, 'Visible', 'off');
    end

%% Function controlling the closing of the main window
%% ===================================================

    function sptPALM_viewer_closefunction(~,~)
        delete(h.sptPALM_DisplayMovie)
        delete(h.sptPALM_ControlPanel)
    end

%% Function for clearing the display axis
%% ======================================

    function clear_display_axis(~,~)
        figure(h.sptPALM_DisplayMovie)
        clf
        h.MainAxes = axes('Parent', h.sptPALM_DisplayMovie, 'box', 'on', 'XTick', [], 'YTick', [], ...
    'Position', [Corner_x,Corner_y,Lx,Ly]);
    end
end