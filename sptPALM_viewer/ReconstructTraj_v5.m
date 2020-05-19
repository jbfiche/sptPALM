%*****************************
%
% ReconstructTraj_v5.m
%
% ****************************
%
% JB Fiche
% Feb, 2020
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function is reading all the selected mat files created by
% MTT and creating a variable called "Reconstructed_Traj" where all the
% trajectories (even the single events) are saved. 
% Also, if the "Save trajectories for Tesseler" is checked, two txt files
% are created. One with the X,Y positions of ALL the events detected bt
% MTT. The other with the mean X,Y positions of each trajectories. 
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.



function h = ReconstructTraj_v5(h)

FileToAnalyse = h.FileToAnalyse;
AcquisitionTime = str2double(get(h.AcquisitionTime, 'String')); % in ms
PixelSize = str2double(get(h.PixelSize, 'String')); % in �m
FontSize = h.FontSize;

%% Check whether all the detections should be saved
%% in a separate .txt file for Tesseler analysis
%% =============================================

CreateTxtFile = h.Save_traj_txt.Value;

%% Open the figure where the results will be saved
%% ===============================================

hPlot = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen

%% For each MTT file, the trajectories are analyzed and saved in "Reconstructed_Traj"
%% ==================================================================================

Reconstructed_Traj = {};
Localizations_all = [];
Localizations_all_average = [];
Length_Traj = [];
SingleStep_Length = [];

tic
hwaitbar = waitbar(0,strcat('Calculating the trajectories for file #1'));

for Nfiles = 1 : numel(FileToAnalyse)
    
    waitbar(Nfiles/numel(FileToAnalyse), hwaitbar, strcat('Calculating the trajectories for file #', num2str(Nfiles)));
    
    clear('m', 'Traj', 'X', 'Y', 'Frame', 'Blink');
    m = matfile(FileToAnalyse{Nfiles}); % Load the results of the MTT analysis
    Xmatrix = m.Xmatrix;
    Ymatrix = m.Ymatrix;
    Imatrix = m.alphamatrix;
    NTrajectory = size(Xmatrix, 1);
    NFrames = size(Xmatrix, 2);
    
    % The arrays Xmatrix and Ymatrix are analyzed in order to detect the
    % presence of NaN. They usually means that for a specific frame, there
    % was non detection and the columns is filled with NaN and 0... not
    % clear why there is a mixture of both though!
    % In order to avoid the creation of artefactual trajectories (very long
    % ones), all the 0 of the columns are converted to NaN
    % ----------------------------------------------------
    
    for ncol = 1 : NFrames
        
        NaN_X = sum(isnan(Xmatrix(:,ncol)));
        NaN_Y = sum(isnan(Ymatrix(:,ncol)));
        
        if NaN_X>0 || NaN_Y>0
            
            Xmatrix(:,ncol) = NaN;
            Ymatrix(:,ncol) = NaN;
        end
    end
    
    % The NaN can introduced errors when calculating the distance D. For
    % particular experiments where the number of detected particles is
    % really low, there will be a lot of events were a NaN will be added to
    % the X/Ymatrix files (this is because a NaN is added when no particle
    % was detected). In that case, when a NaN is detected at frame
    % n>1, it will be automatically replaced by the value detected for
    % frame n-1
    % ---------
    
    NaN_frames = find(isnan(Xmatrix(1,:)));
    if ~isempty(NaN_frames)
        for nframes = 1 : size(NaN_frames,2)
            
            if NaN_frames(nframes) == 1
                Xmatrix(:,NaN_frames(nframes)) = zeros(size(Xmatrix,1),1);
            else
                Xmatrix(:,NaN_frames(nframes)) = Xmatrix(:,NaN_frames(nframes)-1);
            end
        end
    end
    
    NaN_frames = find(isnan(Ymatrix(1,:))==1);
    if ~isempty(NaN_frames)
        for nframes = 1 : size(NaN_frames,2)
            
            if NaN_frames(nframes) == 1
                Ymatrix(:,NaN_frames(nframes)) = zeros(size(Ymatrix,1),1);
            else
                Ymatrix(:,NaN_frames(nframes)) = Ymatrix(:,NaN_frames(nframes)-1);
            end
        end
    end
    
    % For the files generated on the MARS platforme, they are usually
    % called AC0, AC1, AC2, ... and can therefore be orderd as a function
    % of acquisition time.
    % -------------------
    
    [~,MTTFileName,~] = fileparts(FileToAnalyse{Nfiles});
    if isequal(MTTFileName(1:2), 'AC')
        FileNumber = str2double(MTTFileName(3:4));
    else
        FileNumber = NFiles-1;
    end
    
    for ntraj = 1 : NTrajectory
        
        X = PixelSize*Xmatrix(ntraj,:);
        Y = PixelSize*Ymatrix(ntraj,:);
        I = Imatrix(ntraj,:);
        
        % The way MTT works,can be described as follows:
        %
        % 1- When a particle is detected at the frame #n, the value assigned
        % to Xmatrix and Ymatrix for ALL the frames before n (1 : n-1) is
        % set to 0.
        % 2- When a particle is lost, the soft keep adding the same
        % positions to the output .txt file except if there was no detection
        % at all. In that case, a NaN is added... for this specific
        % situation, all the NaN are then replaced by the value of the
        % preciding column.
        
        % As a result, when calculating the distance between successive
        % detections, most of the trajectories have a D equal to 0. Therefore,
        % below, each trajectory is defined by the points comprised between
        % the very first detection (the first non-zero values of X and Y)
        % and the last non-zero detection (that is, the last non-zero value
        % of D)
        
        % In parallel, the step length between two consecutive detections
        % is also calculated for each trajectory and saved in
        % "SingleStep_Length".
        % --------------------
        
        Idx_FirstDetection = find(X>0 & Y>0,1);
        D = sqrt((X(Idx_FirstDetection+1:end)-X(Idx_FirstDetection:end-1)).^2 + (Y(Idx_FirstDetection+1:end)-Y(Idx_FirstDetection:end-1)).^2);
        Idx = find(D>0);
        
        if size(Idx,2)>=2 && ~isempty(Idx_FirstDetection)
            
            D = D(1:Idx(end));
            Reconstructed_Traj = cat(1, Reconstructed_Traj,...
                cat(1, NFrames*FileNumber+(Idx_FirstDetection:Idx_FirstDetection+Idx(end)), X(Idx_FirstDetection:Idx_FirstDetection+Idx(end)), Y(Idx_FirstDetection:Idx_FirstDetection+Idx(end)), cat(2, 0, D)));
            Length_Traj = cat(1, Length_Traj, (Idx(end)+1)*AcquisitionTime/1000);
            
            for k = 1 : size(D,2)
                if k==1 && D(k)>0
                    SingleStep_Length = cat(1, SingleStep_Length, D(k));
                elseif k>1 && D(k)>0 && D(k-1)>0
                    SingleStep_Length = cat(1, SingleStep_Length, D(k));
                end
            end
        end
        
        % In case it is necessary, all the
        % detections are saved below in a .txt file
        % -----------------------------------------
        
        if CreateTxtFile
            
            if size(Idx,2)>0 && ~isempty(Idx_FirstDetection)
                
                x = transpose(X(Idx_FirstDetection:Idx_FirstDetection+Idx(end)));
                y = transpose(Y(Idx_FirstDetection:Idx_FirstDetection+Idx(end)));
                i = transpose(I(Idx_FirstDetection:Idx_FirstDetection+Idx(end)));
                t = transpose(NFrames*FileNumber+(Idx_FirstDetection:Idx_FirstDetection+Idx(end)));
                idx = i(:)>0;
                
                L = sum(idx);
                x_av = sum(x(idx))/L;
                y_av = sum(y(idx))/L;
                i_av = sum(i(idx))/L;
                
                Localizations_all = cat(1, Localizations_all, [x(idx),y(idx), i(idx),t(idx)]);
                Localizations_all_average = cat(1, Localizations_all_average, [x_av, y_av, i_av, L]);
                
            elseif size(Idx,2)==0 && ~isempty(Idx_FirstDetection)
                
                x = X(Idx_FirstDetection);
                y = Y(Idx_FirstDetection);
                i = I(Idx_FirstDetection);
                t = NFrames*FileNumber+Idx_FirstDetection;
                idx = i(:)>0;
                
                Localizations_all = cat(1, Localizations_all, [x(idx),y(idx), i(idx),t(idx)]);
                Localizations_all_average = cat(1, Localizations_all_average, [x, y, i, 1]);
                
            end
        end
    end
end

close(hwaitbar);
toc

h.SingleStep_Length = SingleStep_Length;
h.Length_Traj = Length_Traj;
h.Reconstructed_Traj = Reconstructed_Traj;
set(h.NTrajectories, 'String', num2str(size(Reconstructed_Traj,1))); % Display the # of trajectories detected

%% If needed, all the positions are saved in a .txt file
%% =====================================================

if CreateTxtFile
    
    FileID = fopen('Localizations.txt', 'w+');
    fprintf(FileID,'x,y,intensity,frame\r\n');
    
    [~, Idx] = sort(Localizations_all(:,4));
    Localizations_all = Localizations_all(Idx,:);
    dlmwrite('Localizations.txt', Localizations_all,'-append','precision', '%f', 'newline','pc')
    
    fileID_av = fopen('Localizations_average.txt', 'w+');
    fprintf(fileID_av,'x_mean,y_mean,intensity_mean,traj_length\r\n');
    
    [~, Idx] = sort(Localizations_all_average(:,4));
    Localizations_all_average = Localizations_all_average(Idx,:);
    dlmwrite('Localizations_average.txt', Localizations_all_average,'-append','precision', '%f', 'newline','pc')
end

%% Calculate and plot the empirical cumulative distribution of the length
%% step in order to check whether the tracking parameters are properly
%% defined or they tend to artificially shorten the trajectories
%% =============================================================

figure(hPlot)
hold off
cla
ax = gca;

[f1,x1] = ecdf(SingleStep_Length);

Idx_0p995 = find(f1>0.995,1);
MaxStepLength_0p995 = x1(Idx_0p995);
Idx_0p99 = find(f1>0.99,1);
MaxStepLength_0p99 = x1(Idx_0p99);
Idx_0p985 = find(f1>0.985,1);
MaxStepLength_0p985 = x1(Idx_0p985);
Idx_0p98 = find(f1>0.98,1);
MaxStepLength_0p98 = x1(Idx_0p98);
[f2,x2] = ecdf(SingleStep_Length(SingleStep_Length<MaxStepLength_0p995));
[f3,x3] = ecdf(SingleStep_Length(SingleStep_Length<MaxStepLength_0p99));
[f4,x4] = ecdf(SingleStep_Length(SingleStep_Length<MaxStepLength_0p985));
[f5,x5] = ecdf(SingleStep_Length(SingleStep_Length<MaxStepLength_0p98));

plot(x1,f1, '-r', 'LineWidth', 1)
hold on
plot(x2,f2, '-', 'Color', [1 0.5 0], 'LineWidth', 1)
plot(x3,f3, '-g', 'LineWidth', 1)
plot(x4,f4, '-', 'Color', [0 0.5 1], 'LineWidth', 1)
plot(x5,f5, '-b', 'LineWidth', 1)

axis square
axis([0 max(x2) 0 1])
box on
ax.FontSize = FontSize;
xlabel('Step length (um)')
ylabel('Cumulative distribution')
title('Cumulative distribution of the step length')
legend(sprintf('All values, Lmax = %.2f um', max(x1)), ...
    sprintf('All values without the 0.5%% longest, Lmax = %.2f um', max(x2)), ...
    sprintf('All values without the 1%% longest, Lmax = %.2f um', max(x3)), ...
    sprintf('All values without the 1.5%% longest, Lmax = %.2f um', max(x4)), ...
    sprintf('All values without the 2%% longest, Lmax = %.2f um', max(x5)),'Location', 'southeast');

saveas(hPlot, 'Cumulative_Distribution_LengthStep.png');

% % Replot the zoom on the part representing the 10 last percents of the
% % cumulative distribution
% % -----------------------
% 
% axis([0 max(x2) 0.9 1])
% axis square
% 
% saveas(hPlot, 'Cumulative_Distribution_LengthStep_ZOOM.png');

%% Plot the distribution of the trajectories lengths. Since some
%% trajectories can be very long as compared to the vast majority of the
%% trajectories, the graph can sometimes be streched out and make the
%% visualization of the distribution quite complicated. To avoid this issue,
%% the binning is limited to the range of lengths representing 99% of the
%% trajectories. The 1% longest remaining are artificially combined in the
%% last bin.
%% ========

figure(hPlot)
hold off
cla
ax = gca;

Length_Traj = sort(Length_Traj);
Max99p = Length_Traj(round(length(Length_Traj)*0.99));
Bin = 0 : AcquisitionTime/1000 : Max99p;
[counts,~] = hist(Length_Traj, Bin);

Cumul = 0;
for n1 = 1 : size(Bin,2)
    Cumul = Cumul + counts(n1)/sum(counts);
    if Cumul>0.8
        T80p = Bin(n1);
        break
    end
end

for n2 = n1 : size(Bin,2)
    Cumul = Cumul + counts(n2)/sum(counts);
    if Cumul>0.9
        T90p = Bin(n2);
        break
    end
end

% Idx_NonZero = find(counts>0);
% X = Bin(Idx_NonZero);
% Y = counts(Idx_NonZero);
% fitobject = fit(X', Y', 'exp1');
% BinFit = Bin(Idx_NonZero(1)) : 0.01 : Max99p;

plot(Bin, counts, '--ob', 'MarkerSize', 5, 'LineWidth', 0.5)
hold on
plot([T80p, T80p], [0, max(counts)], '--g', 'LineWidth', 0.5)
plot([T90p, T90p], [0, max(counts)], '--r', 'LineWidth', 0.5)
% plot(BinFit, fitobject(BinFit), 'Color', [0.3 0.3 0.3])

ax.FontSize = FontSize;
axis([0 Max99p 0 max(counts)])
axis square
title('Trajectories duration distribution')
xlabel('Trajectories duration (s)')
ylabel('Counts')
legend('Length distribution', '80% limit', '90% limit', 'Location', 'northeast')
saveas(hPlot, 'Trajectories_duration.png');

close(hPlot)
