function [h, Area, Reconstructed_Traj_ROI, NTraj] = Select_Trajectories_ROI(h, Results, ax, Reconstructed_Traj, PixelSize)

%% In the previous version of the software, a "ROI.mat" was created in order
%% to save the ROI. It was abandonned later but, in order to avoid redifining
%% the ROIs for the older analysis, the program is checking whether there is 
%% a "ROI.mat" file already in the folder. If there is, the ROI are loaded and
%% saved in the structure called "Results"
%% =======================================

if ~isempty(dir('ROI.mat'))
    load('ROI.mat')
    Results.ROI = ROI;
    h.ROI = ROI;
    delete('ROI.mat')
end

if ~isempty(dir('Average_Image.mat'))
    load('Average_Image.mat')
    Results.AvIm = AvIm;
    h.AvIm = AvIm;
    delete('Average_Image.mat')
end

%% When needed ROIs can be defined manually in order to
%% select only the trajectories that are within the ROIs
%% =====================================================

if isfield(h, 'ROI') && isfield(h, 'AvIm')
    ROIFound = 1;
    AvImFound = 1;
else
    ROIFound = 0;
    AvImFound = 0;
end

if ~ROIFound || ~AvImFound
    DefineROI = questdlg('Do you want to define an ROI for the analysis?', 'Define ROI', 'Define ROI', 'No', 'Define ROI');
else
    DefineROI = questdlg('Do you want to load previously saved ROI for the analysis?', 'Load ROI', 'Load ROI', 'Define ROI', 'No', 'Load ROI');
end

switch DefineROI
    
    case 'Define ROI'
        
        axes(ax)
        hold off
        cla
        
        % Select the image and, if it is a movie,
        % calculate the average movie
        % ----------------------------
        
        if ~AvImFound
            AvIm = CalculateAverageImage;
        else
            AvIm = h.AvIm;
        end
        
        imagesc(AvIm)
        axis image
        axis off
        colormap('Gray')
        hold on
        
        % Plot the average position of the trajectories
        % ---------------------------------------------
        
        NTraj_init = size(Reconstructed_Traj, 1);
        Traj_average = zeros(3,NTraj_init);
        
        for ntraj = 1 : NTraj_init
            Traj_average(:,ntraj) = mean(Reconstructed_Traj{ntraj},2);
        end
        plot(Traj_average(3,:)/PixelSize, Traj_average(2,:)/PixelSize, 'or')
        
        % Define the ROI by clicking directly on the image
        % ------------------------------------------------
        
        ROI = {};
        nROI = 1;
        Clicking = 1;
        NPoint = 0;
        X = [];
        Y = [];
        Color = {[1 0 0], [0 0 1], [0 1 0], [1 0.5 0]};
        
        while Clicking
            
            if NPoint == 0
                hdlg = helpdlg({'Clik left on the image to define the ROI', 'Click right to stop the process.'});
                uiwait(hdlg);
                [x,y, Button] = ginput(2);
            else
                [x,y, Button] = ginput(1);
            end
            
            NPoint = NPoint + 1;
            X = cat(1, X, x);
            Y = cat(1, Y, y);
            
            if Button ==1
                line(X, Y, 'Color', Color{nROI})
            elseif Button == 3
                X(end+1) = X(1);
                Y(end+1) = Y(1);
                ROI{nROI} = [X, Y];
                line(X, Y, 'Color', Color{nROI}, 'LineWidth', 2)
                DefineROI = questdlg('Do you want to define another ROI for the analysis?', 'Define ROI', 'Yes', 'No', 'Yes');
                switch DefineROI
                    case 'Yes'
                        X = [];
                        Y = [];
                        nROI = nROI + 1;
                        
                    case 'No'
                        Clicking = 0;
                end
            end
        end
        
        Results.ROI = ROI; % Save the coordinates of the ROI
        Results.AvIm = AvIm; % Save the average image
        h.ROI = ROI; % Save the coordinates of the ROI
        h.AvIm = AvIm; % Save the average image
        
    case 'Load ROI'
        
        NTraj_init = size(Reconstructed_Traj, 1);
        Traj_average = zeros(3,NTraj_init);
        
        for ntraj = 1 : NTraj_init
            Traj_average(:,ntraj) = mean(Reconstructed_Traj{ntraj},2);
        end
        
        ROI = h.ROI;
        AvIm = h.AvIm;
        nROI = size(ROI,2);
        
        axes(ax)
        hold off
        cla
        
        imagesc(AvIm)
        hold on
        axis image
        axis off
        colormap('Gray');
        
    case 'No'
        
        nROI = 0;
        if ROIFound && AvImFound
            h = rmfield(h, 'ROI');
            h = rmfield(h, 'AvIm');
        end
end

% Calculate the area of the ROI
% ----------------------------

if nROI>0
    Area = zeros(1, nROI);
    
    for nroi = 1 : nROI
        
        roi = ROI{nroi};
        Area(nroi) = round(0.5*abs(sum(roi(1:end-1,1).*roi(2:end,2)) - sum(roi(1:end-1,2).*roi(2:end,1))));
    end
    
else
    Area = NaN;
end

% Select only the trajectories whose average positions are within the ROI
% ------------------------------------------------------------------------

if nROI > 0
    
    IN_all_ROI = zeros(1,NTraj_init);
    for roi = 1 : nROI
        
        X_ROI = ROI{roi}(:,1);
        Y_ROI = ROI{roi}(:,2);
        
        IN = inpolygon(Traj_average(3,:)/PixelSize, Traj_average(2,:)/PixelSize, X_ROI, Y_ROI);
        IN_all_ROI = IN_all_ROI + IN;
    end
    
    Reconstructed_Traj_ROI = Reconstructed_Traj(IN_all_ROI==1);
    Traj_average = Traj_average(:,IN_all_ROI==1);
    NTraj = size(Reconstructed_Traj_ROI, 1);
    
    plot(Traj_average(3,:)/PixelSize, Traj_average(2,:)/PixelSize, 'ob')
    
else
    Reconstructed_Traj_ROI = Reconstructed_Traj;
    NTraj = size(Reconstructed_Traj_ROI, 1);
end