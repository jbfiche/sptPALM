function [MSD_all,Reconstructed_Traj_MSD_accepted,Dapp] = Diff_calculation_parallel_computing(DiffCalculationMethod,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD)

disp(' ')
disp('Calculating the apparent diffusion coefficient ...')

Lag = 1 : 1 : p;
Lag = Lag*AcquisitionTime/1000;

Dapp = zeros(size(MSD_all,1),1);
Idx_MSD_accepted = zeros(size(MSD_all,1),1);

parfor nMSD = 1 : size(MSD_all,1)
    
    MSD = MSD_all{nMSD}(1:p);
    Weight = MSD_weight{nMSD}(1:p);
    Nnan = sum(isnan(MSD));
    Nzeros = sum(MSD(:)>0);
    
    if DiffCalculationMethod == 2 && size(MSD,2)>=p
        
        % Calculation of Dapp using the fit
        % ---------------------------------
        
        if Nzeros == p && Nnan == 0
            
            %                 figure(1)
            %                 cla
            %                 hold on
            %                 plot(Lag', MSD_all(nMSD, 1 : Minp+(p-1))', 'o')
            %                 plot(0, 0, 'o')
            
            [fitobject,~] = fit(Lag', MSD', 'Poly1', 'Weights', Weight', 'Upper', [Inf, min(MSD)]);
            if fitobject.p1>0
                Dapp(nMSD) = fitobject.p1/4; % Return the coefficient of diffusion
                %                     Dapp(end,2) = fitobject.p2; % Return the dynamic localization uncertainty (4 s^2)
                %                     D(end+1,1) = fitobject.p1; % Return the coefficient of diffusion
                %                     D(end,2) = sqrt(fitobject.p2); % Return the dynamic localization uncertainty (4 s^2)
                %                         fitobject = ezfit(Lag', MSD_all(nMSD, 1 : Minp+(p-1) )', 'poly1'); % Replaced the function fit by ezfit to avoid licence problem
                %                         if fitobject.m(2)>0
                %                             D(end+1,1) = fitobject.m(2)/4; % Return the coefficient of diffusion
                %                             D(end,2) = sqrt(fitobject.m(1)/4); % Return the dynamic localization uncertainty (4 s^2)
                Idx_MSD_accepted(nMSD) = 1; % Save the index of the MSD values that are used for the calculation of the apparent diffusion coefficient
            end
        end
        
    else
        
        % Calculation of Dapp using the average
        % -------------------------------------
        
        if Nzeros == p && Nnan == 0
            Dapp(nMSD) = sum(MSD.*Weight./(4*Lag))/sum(Weight);
            Idx_MSD_accepted(nMSD) = 1;
        end
        
    end
end

MSD_all = MSD_all(Idx_MSD_accepted==1); % Keep only the points in MSD_all that were used to calculate the distribution of apparent diffusion coefficient
Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD(Idx_MSD_accepted==1); % Keep only the trajectories that associated to the accepted MSD
disp('Diff calculation done')
disp(' ')