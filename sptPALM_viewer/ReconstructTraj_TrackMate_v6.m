%*****************************
%
% ReconstructTraj_v6.m
%
% ****************************
%
% JB Fiche
% Creation 2020
% Last update : 2020/06/04
%
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: This function is reading all the selected TrackMate xml output 
% files and creating a variable called "Reconstructed_Traj" where all the
% trajectories (even the single events) are saved.
% Also, if the "Save trajectories in txt files" is checked, two txt files
% are created. One with the X,Y positions of ALL the events detected bt
% MTT. The other with the mean X,Y positions of each trajectories.
% -------------------------------------------------------------------------
% Specific:
% -------------------------------------------------------------------------
% To fix:
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function h = ReconstructTraj_TrackMate_v6(h)

AcquisitionTime = str2double(get(h.AcquisitionTime, 'String')); % in ms
PixelSize = str2double(get(h.PixelSize, 'String')); % in �m
FontSize = h.FontSize;
ax = h.MainAxes;

%% Check whether all the detections should be saved
%% in a separate .txt file for Tesseler analysis
%% =============================================

CreateTxtFile = h.Save_traj_txt.Value;

%% Ask for the number of frames 
%% ============================

prompt = {'Enter number of frames analyzed for each TrackMate files :'};
dlgtitle = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput);

Nframe = str2double(answer{1});

%% Load all the trajectories and save them in "Reconstructed_Traj" with the right format
%% ======================================================================================

TrackMate = h.TrackMate;
NTrajectory = h.Total_tracks;

h = rmfield(h,{'TrackMate', 'Total_tracks'});

Reconstructed_Traj = cell(NTrajectory,1);
Length_Traj = zeros(NTrajectory,1);
Accepted_tracks = zeros(NTrajectory,1);
SingleStep_Length = [];
Localizations_all = [];
Localizations_all_average = zeros(NTrajectory,4);

SavedTracks = 0;

for nfile = 1 : size(TrackMate,1)
    
    fprintf('\n Formating the trajectories of file #%i ...     ',nfile)
    m = TrackMate{nfile};
    NTraj = size(m,1);
    
    for ntraj = 1 : NTraj
        
        fprintf('\b\b\b\b%03i%%', round(100*ntraj/NTraj))
        
        T = m{ntraj}(:,1) + Nframe*(nfile-1);
        X = PixelSize*m{ntraj}(:,2);
        Y = PixelSize*m{ntraj}(:,3);
        D = sqrt( (X(2:end) - X(1:end-1)).^2 + (Y(2:end) - Y(1:end-1)).^2 );
        
        % Sometimes TrackMate is assigning wrong positions to the events
        % detected too close to the image boundaries, leading to errors in
        % the values of D (either equal to 0 or 1). These errors appear as
        % steps in the cumulative distribution of the single step-lengths. 
        % These events are automatically removed from the analysis below.
        % 
        % Several trajectories are also displaying several times the same
        % localizations (either in X or Y). Those trajectories are also
        % removed since the localizations are wrong.
        % ------------------------------------------
        
        if sum(D==0)==0
            
            L = T(end)-T(1);
            x_av = sum(X)/length(X);
            y_av = sum(Y)/length(Y);
            
            Reconstructed_Traj{ntraj+SavedTracks} = cat(1,T',X',Y',cat(2, 0, D'));
            Length_Traj(ntraj+SavedTracks) = L*AcquisitionTime/1000;
            Accepted_tracks(ntraj+SavedTracks) = 1;
            
            SingleStep_Length = cat(1, SingleStep_Length, D);
            
            if CreateTxtFile
                Localizations_all = cat(1, Localizations_all, cat(2, m{ntraj}, zeros(size(m{ntraj},1),1)));
                Localizations_all_average(ntraj+SavedTracks,:) = [x_av, y_av, 0, L];
            end
        end
    end
    
    SavedTracks = SavedTracks + NTraj;
    fprintf('\n')
end

h.SingleStep_Length = SingleStep_Length;
h.Length_Traj = Length_Traj(Accepted_tracks==1);
h.Reconstructed_Traj = Reconstructed_Traj(Accepted_tracks==1);
set(h.NTrajectories, 'String', num2str(size(Reconstructed_Traj,1))); % Display the # of trajectories detected

%% If needed, all the positions are saved in a .txt file
%% =====================================================

if CreateTxtFile
    
    hwarn = warndlg('Keep in mind that TrackMate does not return the localizations of events localized only once. The minimum number of events for the trajectories is equal to 2.');
    uiwait(hwarn)
    delete(hwarn)
    
    [~, Idx] = sort(Localizations_all(:,4));
    Localizations_all = Localizations_all(Idx,:);
    Localizations_all = array2table(Localizations_all, 'VariableNames', {'x','y','intensity','frame'});
    
    writetable(Localizations_all, 'Localizations.txt', 'Delimiter', 'space')
    
    [~, Idx] = sort(Localizations_all_average(:,4));
    Localizations_all_average = Localizations_all_average(Idx,:);
    Localizations_all_average = array2table(Localizations_all_average, 'VariableNames', {'x_mean','y_mean','intensity_mean','traj_length'});
    
    writetable(Localizations_all_average, 'Localizations_average.txt', 'Delimiter', 'space')
end

%% Calculate and plot the empirical cumulative distribution of the length
%% step in order to check whether the tracking parameters are properly
%% defined or they tend to artificially shorten the trajectories
%% =============================================================

axes(ax)
hold off
cla

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

saveas(ax, 'Cumulative_Distribution_LengthStep.png');

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

axes(ax)
hold off
cla

Length_Traj = sort(Length_Traj);
Max99p = Length_Traj(round(length(Length_Traj)*0.99));
bin = 0 : AcquisitionTime/1000 : Max99p;
hist = histogram(Length_Traj, bin, 'Normalization', 'probability');

counts = hist.Values;
bin = (hist.BinEdges(2:end) + hist.BinEdges(1:end-1))/2;

Cumul = 0;
for n1 = 1 : size(bin,2)
    Cumul = Cumul + counts(n1)/sum(counts);
    if Cumul>0.8
        T80p = bin(n1);
        break
    end
end

for n2 = n1 : size(bin,2)
    Cumul = Cumul + counts(n2)/sum(counts);
    if Cumul>0.9
        T90p = bin(n2);
        break
    end
end

% Idx_NonZero = find(counts>0);
% X = Bin(Idx_NonZero);
% Y = counts(Idx_NonZero);
% fitobject = fit(X', Y', 'exp1');
% BinFit = Bin(Idx_NonZero(1)) : 0.01 : Max99p;

% plot(Bin, counts, '--ob', 'MarkerSize', 5, 'LineWidth', 0.5)
hold on
plot([T80p, T80p], [0, max(counts)], '--g', 'LineWidth', 0.5)
plot([T90p, T90p], [0, max(counts)], '--r', 'LineWidth', 0.5)
% plot(BinFit, fitobject(BinFit), 'Color', [0.3 0.3 0.3])

ax.FontSize = FontSize;
axis([0 Max99p 0 max(counts)])
axis square
title('Trajectories duration distribution')
xlabel('Trajectories duration (s)')
ylabel('Fraction of trajectories')
legend('Length distribution', '80% limit', '90% limit', 'Location', 'northeast')
saveas(ax, 'Trajectories_duration.png');
