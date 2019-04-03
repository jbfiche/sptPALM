function [NbrLorentzianFit, D_mean, varargout] = Simulation_FitLorentzianDistribution_v1(N, Bin, LogDapp, MSD_all, FontSize, hPlot, Lmax, Traj)

NbrLorentzianFit = questdlg('Do you want to fit the distribution with one or two Lorentzian?', 'Fit', '1', '2', '1');
NbrLorentzianFit = str2num(NbrLorentzianFit);

[A0,Idx] = max(N); % Estimate the starting parameters for the first fit
A0 = double(A0);
X0 = double(Bin(Idx));
g0 = double(1);

Lorentz = fittype( @(g,x0,A,x) A./(1 + ((x - x0)/g).^2) );
[fitobject,gof1] = fit(double(Bin)', double(N)', Lorentz, 'start', [g0,X0,A0]);
LorentzFit = Lorentz(fitobject.g, fitobject.x0, fitobject.A, Bin);

if NbrLorentzianFit == 2
    
    % Fit the two distributions with a Lorentzian function and plot the two
    % fits as well as the sum of the two on the graph. 
    % ------------------------------------------------
    
    warnh = warndlg('Click on the approximate positions of the two maxima.');
    uiwait(warnh)
    [X,Y] = ginput(2);
    
    g0 = 0.2;
    X0 = X(1);
    A0 = Y(1);
    g1 = 0.2;
    X1 = X(2);
    A1 = Y(2);
    
    Lorentz2 = fittype( @(g1,x01,A1,g2,x02,A2,x) A1./(1 + ((x - x01)/g1).^2) + A2./(1 + ((x - x02)/g2).^2) );
    [fitobject2,gof2] = fit(double(Bin)', double(N)', Lorentz2, 'start', [g0,X0,A0,g1,X1,A1]);
    disp(strcat('for the fit, R²=', num2str(100*gof2.rsquare), '%'))
    
    % Find the intersection point between the two fits (estimation)
    % -------------------------------------------------------------
    
    BinFit = min(Bin) : 0.01 : max(Bin);
    LorentzFit1 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, BinFit);
    LorentzFit2 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, BinFit);
    
    if fitobject2.x01 < fitobject2.x02
        Idx_Bin = find(BinFit >= fitobject2.x01 & BinFit <= fitobject2.x02);
    else
        Idx_Bin = find(BinFit >= fitobject2.x02 & BinFit <= fitobject2.x01);
    end
    
    Diff = abs(LorentzFit1(Idx_Bin) - LorentzFit2(Idx_Bin));
    [~,Idx] = min(Diff);
    x_Inter = BinFit(Idx_Bin(Idx));
    
    % Analyse separatly the two distributions and return the fraction of
    % events belonging to each of the two distributions. It also splits the
    % trajectories and MSD between the two populations.
    % -------------------------------------------------
       
    if fitobject2.x01 > fitobject2.x02
        D_mean = [fitobject2.x01, fitobject2.x02];
        LorentzFit1 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, BinFit);
        LorentzFit2 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, BinFit);
    else
        D_mean = [fitobject2.x02, fitobject2.x01];
        LorentzFit1 = Lorentz(fitobject2.g2, fitobject2.x02, fitobject2.A2, BinFit);
        LorentzFit2 = Lorentz(fitobject2.g1, fitobject2.x01, fitobject2.A1, BinFit);
    end
    LorentzFitAll = Lorentz2(fitobject2.g1, fitobject2.x01, fitobject2.A1, fitobject2.g2, fitobject2.x02, fitobject2.A2, BinFit);
    
    Idx_1 = LogDapp(:,1) > x_Inter;
    MSD_all_1 = MSD_all(Idx_1);
    Traj_1 = Traj(Idx_1);
    Idx_2 = LogDapp(:,1) < x_Inter;
    MSD_all_2 = MSD_all(Idx_2);
    Traj_2 = Traj(Idx_2);
    MSD_1 = zeros(Lmax, 3); % MSD_1 corresponds to the distribution with the highest diffusion coeff
    MSD_2 = zeros(Lmax, 3); % MSD_2 corresponds to the distribution with the lowest diffusion coeff
    
    LogDapp_1 = LogDapp(Idx_1,1);
    LogDapp_2 = LogDapp(Idx_2,1);
    F = round(100*size(LogDapp_1,1)/(size(LogDapp_1,1)+size(LogDapp_2,1)));
    disp(strcat('The mobile fraction is equal to f=',num2str(F), '%'));
    
    [N1,Bin] = hist(LogDapp_1(:,1), Bin);
    [N2,Bin] = hist(LogDapp_2(:,1), Bin);
    N1_norm = 100*N1/sum(N1+N2);
    N2_norm = 100*N2/sum(N1+N2);
    
    % replot the two distributions, showing now the two separate
    % populations
    % -----------
    
    figure(hPlot)
    ax = gca;
    hold off
    cla
    
    b = bar(Bin', cat(2, N1_norm', N2_norm'), 'stacked');
    b(1).FaceColor = [1 0.5 0];
    b(2).FaceColor = [0 0.4 1];
    
    hold on
    plot(BinFit, LorentzFit1, '-r', 'LineWidth',1)
    plot(BinFit, LorentzFit2, '-b', 'LineWidth',1)
    plot(BinFit, LorentzFitAll, '--k', 'LineWidth',1)
    
    ax.FontSize = FontSize;
    axis square
    box on
    xlabel('Log of apparent diffusion coefficient (µm²/s)')
    ylabel('Fraction of molecule (%)')
    
else
    hold on
    plot(Bin, LorentzFit, '--k', 'LineWidth',1)
    disp(strcat('for the fit, R²=', num2str(100*gof1.rsquare), '%'))
    
    Idx = find(LogDapp(:,1)>fitobject.x0-3*abs(fitobject.g) & LogDapp(:,1)<fitobject.x0+3*abs(fitobject.g));
    MSD_all_1 = MSD_all(Idx);
    MSD = zeros(Lmax, 3);
    D_mean = fitobject.x0;
    
    disp(strcat('X0=', num2str(fitobject.x0), '_ s0=', num2str(fitobject.g)));
end

% Display the title on the graph, indicating when there are two populations
% the two diffusion coefficients as well as the fraction of mobile
% particles.
% ----------

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
    
% When detecting two populations, it splits the MSD into two arrays, 
% according to the 2 populations selected. Also, it makes sure that the MSD
% values that are selected are not equal to zero.
% ----------------------------------------------

for n = 1 : Lmax
    
    if NbrLorentzianFit == 2
        
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

if NbrLorentzianFit == 2
    varargout = {{cat(2, MSD_1, MSD_2), Traj_1, Traj_2, F, cat(1, BinFit, LorentzFit1, LorentzFit2, LorentzFitAll), cat(1, Bin, N1_norm, N2_norm)}};
else
    varargout = {{MSD, Traj, cat(1, Bin, LorentzFit), cat(1, Bin, N)}};
end