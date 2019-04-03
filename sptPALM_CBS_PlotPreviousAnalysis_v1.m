function sptPALM_CBS_PlotPreviousAnalysis_v1(h)

PlotType = get(h.DataTypePlot, 'Value');
NbrLorentzianFit = h.NbrLorentzianFit;

hPlot = figure;
set(0,'Units','pixels'); %Define the type of units used later for the position (here in pixels)
scnsize = get(0,'ScreenSize');%Get the size of the screen in pixels
set(hPlot,'OuterPosition',scnsize);%Display fig1 in order to completely fill the screen
ax = gca;

switch PlotType
    
    % Replot the distribution of diffusion coefficient
    % ------------------------------------------------
    
    case 1
        
        DiffDistribution = h.DiffDistribution;
        FittedDiffDistribution = h.FittedDiffDistribution;
        D_mean = h.D_mean;
        
        switch NbrLorentzianFit
            case 1
                              
                Bin = DiffDistribution(1,:);
                N = DiffDistribution(2,:);
                BinFit = FittedDiffDistribution(1,:);
                LorentzFit = FittedDiffDistribution(2,:);
                
                bar(Bin', N', 'FaceColor', [0 0.4 1]);
                hold on
                plot(BinFit, LorentzFit, '--k', 'LineWidth',1);
                
            case 2
                
                F = h.Fraction;
                
                Bin = DiffDistribution(1,:);
                N1_norm = DiffDistribution(2,:);
                N2_norm = DiffDistribution(3,:);
                BinFit = FittedDiffDistribution(1,:);
                LorentzFit1 = FittedDiffDistribution(2,:);
                LorentzFit2 = FittedDiffDistribution(3,:);
                LorentzFitAll = FittedDiffDistribution(4,:);
                
                b = bar(Bin', cat(2, N1_norm', N2_norm'), 'stacked');
                b(1).FaceColor = [1 0.5 0];
                b(2).FaceColor = [0 0.4 1];
                
                hold on
                plot(BinFit, LorentzFit1, '-r', 'LineWidth',1)
                plot(BinFit, LorentzFit2, '-b', 'LineWidth',1)
                plot(BinFit, LorentzFitAll, '--k', 'LineWidth',1)
        end
        
        ax.FontSize = h.FontSize;
        axis square
        box on
        xlabel('Log of apparent diffusion coefficient (µm²/s)')
        ylabel('Fraction of molecule (%)')
        
        for n = 1 : NbrLorentzianFit
            
            NewLine = sprintf('log(D_%d) = %.2f µm²/s', n, D_mean(1,n));
            if n>1
                Title = strvcat(Title, NewLine);
            else
                Title = NewLine;
            end
            
            if NbrLorentzianFit>1 && n>1
                NewLine = sprintf('Mobile fraction : %d%%', F);
                Title = strvcat(Title, NewLine);
            end
        end
        title(Title);
        
        % Replot the MSD curve
        % --------------------
        
    case 2
        
        MSD_FIT = h.MSD_FIT;
        D_mean = h.D_mean;
        AcquisitionTime = str2num(get(h.AcquisitionTime, 'String'));
        FontSize = h.FontSize;
        
        Lag = 1 : 1 : size(MSD_FIT,1);
        
        
        if NbrLorentzianFit == 2
            
            MSD_1 = MSD_FIT(:,1:3);
            MSD_2 = MSD_FIT(:,4:6);
            
            errorbar(Lag*AcquisitionTime/1000, MSD_1(:,1), MSD_1(:,2), '-o', 'Color', [1 0.5 0])
            hold on
            errorbar(Lag*AcquisitionTime/1000, MSD_2(:,1), MSD_2(:,2), '-o', 'Color', [0 0.4 1])
            %     [fitobject1,gof1] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_1(Idx(1:4),1), 'poly1');
            %     [fitobject2,gof2] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_2(Idx(1:4),1), 'poly1');
            
            errorbar(Lag*AcquisitionTime/1000, MSD_1(:,3), MSD_1(:,2), '-s', 'Color', [1 0.5 0], 'LineWidth', 2)
            errorbar(Lag*AcquisitionTime/1000, MSD_2(:,3), MSD_2(:,2), '-s', 'Color', [0 0.4 1], 'LineWidth', 2)
            %     [fitobject1,gof1] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_1(Idx(1:4),1), 'poly1');
            %     [fitobject2,gof2] = fit(Lag(Idx(1:4))'*AcquisitionTime/1000, MSD_2(Idx(1:4),1), 'poly1');
            
            axis square
            ax.FontSize = FontSize;
            box on
            xlabel('Time(s)')
            ylabel('MSD (µm²)')
            
            Title = sprintf('log(D_1) = %.2f µm²/s -- log(D_2) = %.2f µm²/s', D_mean(1,1), D_mean(1,2));
            title(Title);
            legend('D1 average', 'D2 average', 'D1 median', 'D2 median', 'Location', 'northwest');
            
        else
            
            MSD = MSD_FIT;
            
            errorbar(Lag*AcquisitionTime/1000, MSD(:,1),  MSD(:,2), '-o', 'Color', [0 0.4 1])
            hold on
            errorbar(Lag*AcquisitionTime/1000, MSD(:,3),  MSD(:,2), '-s', 'Color', [0 0.4 1], 'LineWidth', 2)
            [fitobject,~] = fit(Lag(1:4)'*AcquisitionTime/1000, MSD(1:4,1), 'poly1');
            
            Fit = fitobject.p1 * Lag'*AcquisitionTime/1000 + fitobject.p2;
            plot(Lag*AcquisitionTime/1000, Fit, '-r')
            axis square
            axis tight
            box on
            ax.FontSize = FontSize;
            xlabel('Time(s)')
            ylabel('MSD (µm²/s)')
            legend('Mean values', 'Median values', 'Location', 'northwest')
            
        end
        
        
        % Replot the trajectories (all the trajectories selected within the
        % ROI and after applying the filters)
        % -----------------------------------
        
    case 3
        
        Reconstructed_Traj_ROI = h.Reconstructed_Traj_ROI;
        NTraj_ROI = size(Reconstructed_Traj_ROI,1);
        
        if isfield(h, 'AvIm')
            imagesc(h.AvIm)
            colormap('gray')
        else
            fig = gcf;
            fig.Color = [0.6 0.6 0.6];
        end
        
        axis image
        axis off
        legend off
        title('')
        
        Color = jet;
        ntraj_color = ceil(NTraj_ROI/size(Color,1));
        
        for ntraj = 1 : NTraj_ROI
            if isfield(h, 'AvIm')
                PixelSize = str2num(get(h.PixelSize, 'String'));
                X = Reconstructed_Traj_ROI{ntraj}(2,:)/PixelSize;
                Y = Reconstructed_Traj_ROI{ntraj}(3,:)/PixelSize;
            else
                X = Reconstructed_Traj_ROI{ntraj}(2,:);
                Y = Reconstructed_Traj_ROI{ntraj}(3,:);
                if ~isempty(find(X==0))
                    ntraj
                end
            end
            line(Y, X, 'Color', Color(ceil(ntraj/ntraj_color),:),'LineWidth',1)
        end
        
        if isfield(h, 'AvIm')
            
            ScaleBar = [5, 5, 1/PixelSize, 1];
            rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
            
        else
            AxisLimits = axis;
            Box = [AxisLimits(1), AxisLimits(3), AxisLimits(2)-AxisLimits(1), AxisLimits(4)-AxisLimits(3)];
            rectangle('Position', Box, 'EdgeColor', [0 0 0], 'LineWidth', 2)
            
            ScaleBar = [AxisLimits(1)+1, AxisLimits(3)+1, 1, 0.1];
            rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
        end
        
        % Replot the trajectories by selecting only the two populations
        % selected during the fit of the diffusion coefficients
        % distributions
        % -------------
        
    case 4
        
        if NbrLorentzianFit == 1
            hwarn = warndlg('This kind of plot cannot be plotted for this analysis, only one population was fitted');
            uiwait(hwarn)
            close(hPlot)
            return
        end
        
        Reconstructed_Traj_DiffPop1 = h.Reconstructed_Traj_DiffPop1;
        Reconstructed_Traj_DiffPop2 = h.Reconstructed_Traj_DiffPop2;
        NTraj_Pop1 = size(Reconstructed_Traj_DiffPop1,1);
        NTraj_Pop2 = size(Reconstructed_Traj_DiffPop2,1);
        
        if isfield(h, 'AvIm')
            imagesc(h.AvIm)
            colormap('gray')
        else
            fig = gcf;
            fig.Color = [0.6 0.6 0.6];
        end
        axis image
        axis off
        box on
        legend off
        title('')
        
        for ntraj = 1 : NTraj_Pop1
            if isfield(h, 'AvIm')
                PixelSize = str2num(get(h.PixelSize, 'String'));
                X = Reconstructed_Traj_DiffPop1{ntraj}(2,:)/PixelSize;
                Y = Reconstructed_Traj_DiffPop1{ntraj}(3,:)/PixelSize;
            else
                X = Reconstructed_Traj_DiffPop1{ntraj}(2,:);
                Y = Reconstructed_Traj_DiffPop1{ntraj}(3,:);
            end
            line(Y, X, 'Color', [1 0.7 0],'LineWidth',1)
        end
        
        for ntraj = 1 : NTraj_Pop2
            if isfield(h, 'AvIm')
                X = Reconstructed_Traj_DiffPop2{ntraj}(2,:)/PixelSize;
                Y = Reconstructed_Traj_DiffPop2{ntraj}(3,:)/PixelSize;
            else
                X = Reconstructed_Traj_DiffPop2{ntraj}(2,:);
                Y = Reconstructed_Traj_DiffPop2{ntraj}(3,:);
            end
            line(Y, X, 'Color', [0 0.4 1],'LineWidth',1)
        end
        
        if isfield(h, 'AvIm')
            
            ScaleBar = [5, 5, 1/PixelSize, 1];
            rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
            
        else
            AxisLimits = axis;
            Box = [AxisLimits(1), AxisLimits(3), AxisLimits(2)-AxisLimits(1), AxisLimits(4)-AxisLimits(3)];
            rectangle('Position', Box, 'EdgeColor', [0 0 0], 'LineWidth', 2)
            
            ScaleBar = [AxisLimits(1)+1, AxisLimits(3)+1, 1, 0.1];
            rectangle('Position', ScaleBar, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2)
        end
end
