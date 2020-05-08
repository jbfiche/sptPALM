function [NbrGaussianFit, D_mean, varargout] = FitGaussianDistribution_v3(LogDapp, MSD_all, FontSize, ax, Lmax, Traj, DiffCalculationMethod)

% Plot the histogram
% ------------------

MinBin = floor(20*min(LogDapp))/20;
MaxBin = floor(20*max(LogDapp))/20;
Edges = MinBin : 0.05 : MaxBin;

hist = histogram(LogDapp, Edges, 'Visible', 'off');

hist_Values = hist.Values;
Ntot = sum(hist_Values);
hist_Values = hist_Values'*100/Ntot;
hist_Edges = hist.BinEdges;
hist_Bin = (hist_Edges(1:end-1) + hist_Edges(2:end)) / 2;
hist_Bin = hist_Bin';

axes(ax)
hold off
cla
bar(hist_Bin,hist_Values)
ax.FontSize = FontSize;
axis square
box on
xlabel('Log of apparent diffusion coefficient (um^2/s)')
ylabel('Fraction of molecule (%)')

% Ask whether the distributions should be fitted with one or two gaussians
% -------------------------------------------------------------------------

Gaussian = fittype( @(s,x0,A,x) A*exp(-((x - x0)/(2*s)).^2) );
NbrGaussianFit = questdlg('Do you want to fit the distribution with one or two Gaussian?', 'Fit', '1', '2', '1');
NbrGaussianFit = str2double(NbrGaussianFit);

if NbrGaussianFit == 2
    
    % Fit the two distributions with a Gaussian function and plot the two
    % fits as well as the sum of the two on the graph.
    % ------------------------------------------------
    
    warnh = warndlg('Click on the approximate threshold separating the two populations.');
    uiwait(warnh)
    [X,~] = ginput(1);
  
    D0 = LogDapp(LogDapp>X);
    D1 = LogDapp(LogDapp<=X);
    
    s0 = std(D0);
    X0 = mean(D0);
    A0 = hist_Values(find(hist_Bin>X0, 1, 'first'));
    s1 = std(D1);
    X1 = mean(D1);
    A1 = hist_Values(find(hist_Bin>X1, 1, 'first'));
    
    Gaussian2 = fittype( @(s1,x01,A1,s2,x02,A2,x) A1*exp(-((x - x01)/(2*s1)).^2) + A2*exp(-((x - x02)/(2*s2)).^2));
    [fitobject2,gof2] = fit(hist_Bin, hist_Values, Gaussian2, 'start', [s0,X0,A0,s1,X1,A1]);
    disp(strcat('For the Gaussian fit, R^2=', num2str(100*gof2.rsquare), '%'))
    
    % Find the intersection point between the two fits (estimation). By
    % convention, population 1 is the one with the lowest diffusion
    % coefficient.
    % -----------
    
    BinFit = min(hist_Bin) : 0.01 : max(hist_Bin);
    
    if fitobject2.x01 < fitobject2.x02
        D_mean = [fitobject2.x01, fitobject2.x02];
        Idx_Bin = find(BinFit >= fitobject2.x01 & BinFit <= fitobject2.x02);
        GaussianFit1 = Gaussian(fitobject2.s1, fitobject2.x01, fitobject2.A1, BinFit);
        GaussianFit2 = Gaussian(fitobject2.s2, fitobject2.x02, fitobject2.A2, BinFit);
    else
        D_mean = [fitobject2.x02, fitobject2.x01];
        Idx_Bin = find(BinFit >= fitobject2.x02 & BinFit <= fitobject2.x01);
        GaussianFit1 = Gaussian(fitobject2.s2, fitobject2.x02, fitobject2.A2, BinFit);
        GaussianFit2 = Gaussian(fitobject2.s1, fitobject2.x01, fitobject2.A1, BinFit);
    end
    
    Diff = abs(GaussianFit1(Idx_Bin) - GaussianFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter = BinFit(Idx_Bin(Idx));
    
    % Analyse separatly the two distributions and return the fraction of
    % events belonging to each of the two distributions. It also splits the
    % trajectories and MSD between the two populations.
    % -------------------------------------------------
    
    GaussianFitAll = Gaussian2(fitobject2.s1, fitobject2.x01, fitobject2.A1, fitobject2.s2, fitobject2.x02, fitobject2.A2, BinFit);
    
    Idx_1 = LogDapp(:,1) < x_Inter;
    MSD_all_1 = MSD_all(Idx_1);
    Traj_1 = Traj(Idx_1);
    Idx_2 = LogDapp(:,1) > x_Inter;
    MSD_all_2 = MSD_all(Idx_2);
    Traj_2 = Traj(Idx_2);
    MSD_1 = zeros(Lmax, 3); % MSD_1 corresponds to the distribution with the highest diffusion coeff
    MSD_2 = zeros(Lmax, 3); % MSD_2 corresponds to the distribution with the lowest diffusion coeff
    
    LogDapp_1 = LogDapp(Idx_1,1);
    LogDapp_2 = LogDapp(Idx_2,1);
    F = round(100*size(LogDapp_2,1)/(size(LogDapp_1,1)+size(LogDapp_2,1)));
%     disp(strcat('The mobile fraction is equal to f=',num2str(F), '%'));
    
    axes(ax)
    hist1 = histogram(LogDapp_1(:,1), hist_Edges, 'Visible', 'off');
    hist1_Values = hist1.Values;
    hist1_Values = hist1_Values'*100/Ntot;

    hist2 = histogram(LogDapp_2(:,1), hist_Edges, 'Visible', 'off');
    hist2_Values = hist2.Values;
    hist2_Values = hist2_Values'*100/Ntot;
    
    % replot the two distributions, showing now the two separate
    % populations
    % -----------
    
    axes(ax)
    hold off
    cla
    
    b = bar(hist_Bin, cat(2, hist1_Values, hist2_Values), 'stacked');   
    b(1).FaceColor = [0 0.4470 0.7410];
    b(2).FaceColor = [0.8500 0.3250 0.0980];
    hold on
    plot(BinFit, GaussianFit1, '-b', 'LineWidth',1)
    plot(BinFit, GaussianFit2, '-r', 'LineWidth',1)
    plot(BinFit, GaussianFitAll, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel('Log_1_0(D_i_n_s_t)')
    ylabel('Fraction of molecule (%)')
    
else
    
    [A0,Idx] = max(hist_Values); % Estimate the starting parameters for the first fit
    A0 = double(A0);
    X0 = double(hist_Bin(Idx));
    s0 = double(1);
    
    [fitobject,gof1] = fit(hist_Bin, hist_Values, Gaussian, 'start', [s0,X0,A0]);
    GaussianFit = Gaussian(fitobject.s, fitobject.x0, fitobject.A, hist_Bin);
    
    axes(ax)
    hold on
    plot(hist_Bin, GaussianFit, '--k', 'LineWidth',1)
    if isunix
        disp(strcat('For the gaussian fit, R^2=', num2str(100*gof1.rsquare), '%'))
    else
        disp(strcat('For the gaussian fit, R^2=', num2str(100*gof1.rsquare), '%'))
    end
    
    MSD_all_1 = MSD_all( LogDapp(:,1)>fitobject.x0-3*abs(fitobject.s) & LogDapp(:,1)<fitobject.x0+3*abs(fitobject.s) );
    MSD = zeros(Lmax, 3);
    D_mean = fitobject.x0;
    
%     disp(strcat('X0=', num2str(fitobject.x0), '_ s0=', num2str(fitobject.s)));
end

% Display the title on the graph, indicating when there are two populations
% the two diffusion coefficients as well as the fraction of mobile
% particles.
% ----------

for n = 1 : NbrGaussianFit
    
    NewLine = sprintf('log_1_0(D_%d) = %.2f -- D_%d = %.3f um^2/s', n, D_mean(1,n),n,10^(D_mean(1,n)));
    if n>1
        Title = char(Title, NewLine);
    else
        Title = NewLine;
    end
    
    if NbrGaussianFit>1 && n>1
        NewLine = sprintf('Mobile fraction : %d%%', F);
        Title = char(Title, NewLine);
    end
end

title(Title);

% Save the plot
% -------------

if DiffCalculationMethod == 2
    saveas(ax, 'Diffusion_distribution_Fit_Method.png');
else
    saveas(ax, 'Diffusion_distribution_Weighted_Average_Method.png');
end

% When detecting two populations, it splits the MSD into two arrays,
% according to the 2 populations selected. Also, it makes sure that the MSD
% values that are selected are not equal to zero.
% ----------------------------------------------

for n = 1 : Lmax
    
    if NbrGaussianFit == 2
        
        msd = [];
        for k = 1 : size(MSD_all_1,1)
            if size(MSD_all_1{k},2) >= n
                if MSD_all_1{k}(n) > 0
                    
                    msd = cat(1, msd, MSD_all_1{k}(n));
                end
            end
        end
        
        MSD_1(n,1) = mean(msd); % Average value of the MSD
        MSD_1(n,2) = std(msd)/sqrt(length(msd)); % Standard error of the mean
        MSD_1(n,3) = median(msd); % Median value for the MSD
        
        msd = [];
        for k = 1 : size(MSD_all_2,1)
            if size(MSD_all_2{k},2) >= n
                if MSD_all_2{k}(n) > 0
                    
                    msd = cat(1, msd, MSD_all_2{k}(n));
                end
            end
        end
        
        MSD_2(n,1) = mean(msd); % Average value of the MSD
        MSD_2(n,2) = std(msd)/sqrt(length(msd)); % Standard error of the mean
        MSD_2(n,3) = median(msd); % Median value for the MSD
        
    else
        
        msd = [];
        for k = 1 : size(MSD_all_1,1)
            if size(MSD_all_1{k},2) >= n
                if MSD_all_1{k}(n) > 0
                    
                    msd = cat(1, msd, MSD_all_1{k}(n));
                end
            end
        end
        
        MSD(n,1) = mean(msd);
        MSD(n,2) = std(msd)/sqrt(length(msd));
        MSD(n,3) = median(msd);
    end
end

if NbrGaussianFit == 2
    varargout = {{cat(2, MSD_1, MSD_2), Traj_1, Traj_2, F, cat(2, BinFit', GaussianFit1', GaussianFit2', GaussianFitAll'), cat(2, hist_Bin, hist1_Values, hist2_Values)}};
else
    varargout = {{MSD, Traj, cat(2, hist_Bin, GaussianFit), cat(2, hist_Bin, hist_Values)}};
end