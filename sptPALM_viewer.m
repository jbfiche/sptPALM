%*****************************
%
% sptPALM_viewer.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Analyze tracking data obtained from MTT software (Serge et al.
% 2008, Nat. Methods). The program can be used to sort the track, calculate
% the MSD and diffusion coefficients Dinst, plot the trajectories and the
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

MinEditBox = 0;
MaxEditBox = 1;

% Calculate the dimension of the GUI window according to the disposition of
% the displayed windows.
% ----------------------

Side_Border = 15;
Height_Border = 5;
ControlPanel_Width = 600;
ControlPanel_Height = 865;
DisplayMovie_Width = 720;
DisplayMovie_Height = 750;
Axis_Size = 600;

h.sptPALM_DisplayMovie = figure('Visible','on',...
    'Units','pixels',...
    'Position',[round(scrsz(3)/2)-350,round(scrsz(4)/2)-350,DisplayMovie_Width,DisplayMovie_Height],...
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
    'Position',[round(scrsz(3)/2)-250,round(scrsz(4)/2)-425,ControlPanel_Width,ControlPanel_Height],...
    'MenuBar','none',...
    'Name','sptPALM_Control_Panel',...
    'NumberTitle','off',...
    'Toolbar','none',...
    'Resize','off',...
    'Color',defaultBackground,...
    'CloseRequestFcn',@sptPALM_viewer_closefunction);

% Construct the components
% ------------------------

h = sptPALM_components(h, Side_Border, Height_Border, ControlPanel_Width, ControlPanel_Height);

% Reinitialize the values of all the callbacks and define several
% parameters that will be used for the display or saving the results
% ------------------------------------------------------------------

h = sptPALM_initialize(h, 'Reset_all');
h.FontSize = 10;
h.ResultsFileName = 'MTT_sptPALM_analysis.mat';

% Define the axis on the display figure
% -------------------------------------

Corner_x = (DisplayMovie_Width - Axis_Size)/(2*DisplayMovie_Width);
Corner_y = (DisplayMovie_Height - Axis_Size)/(2*DisplayMovie_Height);
Lx = Axis_Size/DisplayMovie_Width;
Ly = Axis_Size/DisplayMovie_Height;

h.MainAxes = axes('Parent', h.sptPALM_DisplayMovie, 'box', 'on', 'XTick', [], 'YTick', [], ...
    'Position', [Corner_x,Corner_y,Lx,Ly]);

% Define the default parameters for the simulation
% ------------------------------------------------

SimulationParameters.NProteinsTot_Init = 10000; % Define the total pool of proteins that can be activated and imaged
SimulationParameters.MeanProbActivation = 0.5; % The probability is defined as the mean number of protein activated per frame
SimulationParameters.MeanPhotoBleachingTime = 1; % in s
SimulationParameters.Ton = 0.3; % in s
SimulationParameters.Toff1 = 0.3; % in s
SimulationParameters.Toff2 = 4; % in s

SimulationParameters.MaxBlink = 3; % in frames
SimulationParameters.MinTrajLength = 5; %in frames

SimulationParameters.ImageSize = 252; % px
SimulationParameters.PixelSize = 0.102; % in �m
SimulationParameters.AcquisitionTime = 0.02; % in s
SimulationParameters.Offset = 90; % intensity offset for the image
SimulationParameters.Noise1 = 60; % parameter lambda for the Poisson distribution describing the noise of the EMCCD camera
SimulationParameters.Noise2 = 0.25; % parameter describing the normal noise for each pixel and linked to the photon detection
SimulationParameters.GaussStampSize = 5; % parameters defining the half-size of the window where the 2D gaussian will be calculated
SimulationParameters.LimitResolution = 0.15; % in �m
SimulationParameters.MeanSingleEventIntensity = 80; % Mean intensity of single emitters

h.SimulationParameters = SimulationParameters;

% Set the Callbacks
% -----------------

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
set(h.Simulation_EmissionParameters,'callback', @Change_EmissionParameters);
set(h.Simulation_AquisitionParameters,'callback', @Change_AquisitionParameters);
set(h.LaunchSimulation,'callback', @LaunchSimulation);

%% Load MTT files
%% ==============

    function LoadMTT(~,~)
        clc
        h = sptPALM_initialize(h, 'Reset_all');
        h = sptPALM_initialize(h, 'Reset_h');
        h.ResultsFileName = h.Saving_file_name.String;
        
        [h, Repeat_Analysis] = Load_MTT_Tracking_Files_v2(h);
        
        if isequal(Repeat_Analysis, 'Start new analysis')
            h = ReconstructTraj_v6(h);
        end
    end

%% Load previous analysis
%% ======================

    function Load_Previous_analysis(~,~)
        clc
        clear_display_axis
        h = sptPALM_initialize(h, 'Reset_all');
        h = sptPALM_initialize(h, 'Reset_h');
        [FileName, PathName] = uigetfile('*.mat');
        
        if FileName ~= 0
            
            cd(PathName);
            h.ResultsFileName = FileName;
            set(h.Saving_file_name, 'String', FileName);
            
            set(h.PlotPreviousAnalysis, 'Enable', 'on')
            h = sptPALM_initialize(h, 'LoadPreviousData');
            
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
        h = MTT_Trajectory_analysis_v4(h);
    end

%% Save the analysis for the Tesseler clustering analysis
%% ======================================================

    function SaveForTesseler(~,~)
        clc
        if h.Parallel_computing.Value == 0
            h = Save_For_Tesseler_sptPALM(h);
        else
            h = Save_For_Tesseler_sptPALM_parallel_computing(h);
        end
    end

%% Plot the previous analysis
%% ==========================

    function PlotPreviousAnalysis(~,~)
        clc
        clear_display_axis
        sptPALM_CBS_PlotPreviousAnalysis_v1(h)
    end

%% Load movie
%% ==========

    function LoadMovie(~,~)
        
        clc
        set(h.sptPALM_DisplayMovie, 'Visible', 'on');
        
        [ImageName, ImageDirectory] = uigetfile('*.tif');
        h.MovieDisplayFullName = strcat(ImageDirectory, ImageName);
        ImInfo = imfinfo(h.MovieDisplayFullName);
        h.MovieDisplayNImages = length(ImInfo);
        h.ImSize = [ImInfo(1).Width, ImInfo(1).Height];
        
        set(h.Slider_SelectFrame, 'Min', 1)
        set(h.Slider_SelectFrame, 'Max', h.MovieDisplayNImages)
        set(h.Slider_SelectFrame, 'Value', 1)
        set(h.Slider_SelectFrame, 'SliderStep', [1/h.MovieDisplayNImages,10/h.MovieDisplayNImages])
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
        
    end

%% Load trajectories associated to the movie loaded
%% ================================================

    function LoadTraj(~,~)
        
        [MTTFileName, MTTFileDirectory] = uigetfile('*.mat');
        FileToAnalyse = strcat(MTTFileDirectory, MTTFileName);
        h = ReconstructTraj_Display_v5(FileToAnalyse, h);
        PlotTrajectories_sptPALM_v1(h)
    end

%% Create an ROI for the displayed movie
%% =====================================

    function CreateROI(~,~)
        
        axes(h.MainAxes)
        h.MovieDisplayROI = getrect(h.MainAxes);
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'MovieDisplayReconstructed_Traj') && isfield(h, 'MovieDisplayFrame_Traj')
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
        
        if isfield(h, 'MovieDisplayReconstructed_Traj') && isfield(h, 'MovieDisplayFrame_Traj')
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
        
        if Up<Low
            Up = Low+0.001;
            set(h.Slider_UpperContrast, 'Value', num2str(Up))
        elseif Up>MaxEditBox
            Up = MaxEditBox;
            set(h.Slider_UpperContrast, 'Value', num2str(Up))
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
        elseif Up>MaxEditBox
            Up = MaxEditBox;
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
        
        if Low<MinEditBox
            Low = MinEditBox;
            set(h.Slider_LowerContrast, 'Value', num2str(Low))
        elseif Low>Up
            Low = Up-0.001;
            set(h.Slider_LowerContrast, 'Value', num2str(Low))
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
        
        if Low<MinEditBox
            Low = MinEditBox;
            set(h.Edit_LowerContrast, 'String', num2str(Low))
        elseif Low>Up
            Low = Up-0.001;
            set(h.Edit_LowerContrast, 'String', num2str(Low))
        end
        
        set(h.Slider_LowerContrast, 'Value', Low)
        
        ChangeContrast_sptPALM_v2(h);
        
        if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
            PlotTrajectories_sptPALM_v1(h)
        end
    end

%% Change the parameters describing the emission properties of the emitters
%% ========================================================================

    function Change_EmissionParameters(~,~)
        
        clc
        % h.SimulationParameters = SimulationParameters;
        
        prompt = {'Enter the total number of emitters available', 'Enter the average number of proteins activated per frame', ...
            'Enter the mean Photo bleaching time (s)', 'Enter the average ON time (s)', 'Enter the average SHORT OFF time (s)', ...
            'Enter the average LONG OFF time (s)', 'Enter the maximum number of frames allowed for the blinking', ...
            'Enter the minimum trajectory length in frames'};
        dlg_title = 'Define emission parameters';
        num_lines= 1;
        default_answer = {num2str(h.SimulationParameters.NProteinsTot_Init), num2str(h.SimulationParameters.MeanProbActivation), ...
            num2str(h.SimulationParameters.MeanPhotoBleachingTime), num2str(h.SimulationParameters.Ton), ...
            num2str(h.SimulationParameters.Toff1), num2str(h.SimulationParameters.Toff2), ...
            num2str(h.SimulationParameters.MaxBlink), num2str(h.SimulationParameters.MinTrajLength)};
        NewEmissionParameters = inputdlg(prompt,dlg_title,num_lines,default_answer);
        
        if ~isempty(NewEmissionParameters)
            h.SimulationParameters.NProteinsTot_Init = str2double(NewEmissionParameters{1}); % Define the total pool of proteins that can be activated and imaged
            h.SimulationParameters.MeanProbActivation = str2double(NewEmissionParameters{2}); % The probability is defined as the mean number of protein activated per frame
            h.SimulationParameters.MeanPhotoBleachingTime = str2double(NewEmissionParameters{3}); % in s
            h.SimulationParameters.Ton = str2double(NewEmissionParameters{4}); % in s
            h.SimulationParameters.Toff1 = str2double(NewEmissionParameters{5}); % in s
            h.SimulationParameters.Toff2 = str2double(NewEmissionParameters{6}); % in s
            h.SimulationParameters.MaxBlink = str2double(NewEmissionParameters{7}); % in frames
            h.SimulationParameters.MinTrajLength = str2double(NewEmissionParameters{8}); %in frames
        end
    end

%% Change the parameters describing the aquisition parameters
%% ==========================================================

    function Change_AquisitionParameters(~,~)
        
        clc
        
        prompt = {'Enter the size of the image in pixels (image is square by default)', 'Enter the aquisition time (s)',...
            'Enter the pixel size(�m)', 'Enter the intensity offset', ...
            'Enter the parameter LAMBDA describing the dark noise of the EMCCD camera (Poisson distribution)',...
            'Enter the parameter SIGMA describing the amplification noise of the EMCCD camera (Normal distribution)',...
            'Enter the half size of the window where the intensity of a single emitter is plot (in px)', ...
            'Enter the diffraction limit (in �m and defined as the standard deviation)',...
            'Enter the mean intensity associated to a single emitters (intensity distribution described by a Poisson distribution)'};
        dlg_title = 'Define aquisition parameters';
        num_lines= 1;
        default_answer = {num2str(h.SimulationParameters.ImageSize), num2str(h.SimulationParameters.AcquisitionTime),...
            num2str(h.SimulationParameters.PixelSize),num2str(h.SimulationParameters.Offset), ...
            num2str(h.SimulationParameters.Noise1), num2str(h.SimulationParameters.Noise2), ...
            num2str(h.SimulationParameters.GaussStampSize), num2str(h.SimulationParameters.LimitResolution), ...
            num2str(h.SimulationParameters.MeanSingleEventIntensity)};
        NewAquisitionParameters = inputdlg(prompt,dlg_title,num_lines,default_answer);
        
        if ~isempty(NewAquisitionParameters)
            h.SimulationParameters.ImageSize = str2double(NewAquisitionParameters{1}); % px
            h.SimulationParameters.AcquisitionTime = str2double(NewAquisitionParameters{2}); % s
            h.SimulationParameters.PixelSize = str2double(NewAquisitionParameters{3}); % �m
            h.SimulationParameters.Offset = str2double(NewAquisitionParameters{4}); % intensity offset for the image
            h.SimulationParameters.Noise1 = str2double(NewAquisitionParameters{5}); % parameter lambda for the Poisson distribution describing the noise of the EMCCD camera
            h.SimulationParameters.Noise2 = str2double(NewAquisitionParameters{6}); % parameter describing the normal noise for each pixel and linked to the photon detection
            h.SimulationParameters.GaussStampSize = str2double(NewAquisitionParameters{7}); % parameters defining the half-size of the window where the 2D gaussian will be calculated
            h.SimulationParameters.LimitResolution = str2double(NewAquisitionParameters{8}); % in �m
            h.SimulationParameters.MeanSingleEventIntensity = str2double(NewAquisitionParameters{9}); % Mean intensity of single emitters
        end
    end

%% Launch the simulation of a sptPALM experiment
%% ===============================================

    function LaunchSimulation(~,~)
       clc
       SimuDiff_v3(h)
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