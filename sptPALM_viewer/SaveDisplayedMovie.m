%*****************************
%
% SaveDisplayedMovie.m
%
% ****************************
%
% JB Fiche
% Last update : 19-05-2020
%
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This option was suggested by C.Royer and first implemented in 
% her lab. The code was modified and adapted to the latest version of 
% sptPALM_viewer. It allows the user to save the currently displayed movie 
% (with or without ROI) as an .avi movie. 
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.

function SaveDisplayedMovie(h)

% Read the values of the first and last frames
% --------------------------------------------

FirstFrame = get(h.FirstFrame, 'String');
LastFrame = get(h.LastFrame, 'String');

if ~isempty(FirstFrame) && ~isempty(LastFrame)
    FirstFrame = str2double(FirstFrame);
    LastFrame = str2double(LastFrame);
else
    hwarn = warndlg('You need to define first the first & last frames');
    uiwait(hwarn)
    delete(hwarn)
    return
end

% Read the name of the movie we want to save
% ------------------------------------------

MovieName = get(h.MovieName, 'String');

% Calculate the image according to the contrast and whether or not a ROI
% was defined.
% ------------

v = VideoWriter(MovieName, 'Uncompressed AVI');
open(v);

for nframe = FirstFrame : 1 : LastFrame
    
    set(h.Slider_SelectFrame, 'Value', nframe);
    set(h.Edit_SelectFrame, 'String', num2str(nframe));
    
    ChangeContrast_sptPALM_v2(h);
    
    if isfield(h, 'Reconstructed_Traj_MovieDisplay') && isfield(h, 'Frame_Traj_MovieDisplay')
        PlotTrajectories_sptPALM_v1(h)
    end
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v)