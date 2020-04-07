function [MSD_all,Reconstructed_Traj_MSD_accepted,Dapp] = Diff_calculation(DiffCalculationMethod,MSD_all,MSD_weight,p,AcquisitionTime,Reconstructed_Traj_MSD)

if DiffCalculationMethod == 1
    fprintf('Starting calculation of the diffusion coefficient with the average method ...     ')
else
    fprintf('Starting calculation of the diffusion coefficient with the fit method ...     ')
end

Lag = 1 : 1 : p;
Lag = Lag*AcquisitionTime/1000;

Dapp = zeros(size(MSD_all,1),1);
Idx_MSD_accepted = zeros(size(MSD_all,1),1);

N_MSD = size(MSD_all,1);

for nMSD = 1 : N_MSD
    
    MSD = MSD_all{nMSD}(1:p);
    Weight = MSD_weight{nMSD}(1:p);
    Nnan = sum(isnan(MSD));
    Nzeros = sum(MSD(:)>0);
    
    fprintf('\b\b\b\b%03i%%', round(100*nMSD/N_MSD))
    
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
            %             Dapp(nMSD) = sum(MSD.*Weight./Lag)/(4*sum(Weight));
            Dapp(nMSD) = 1/4 * (MSD(4)-MSD(1))/(Lag(4)-Lag(1));
            if Dapp(nMSD)>0
                Idx_MSD_accepted(nMSD) = 1;
            end
        end
    end
end

MSD_all = MSD_all(Idx_MSD_accepted==1); % Keep only the points in MSD_all that were used to calculate the distribution of apparent diffusion coefficient
Reconstructed_Traj_MSD_accepted = Reconstructed_Traj_MSD(Idx_MSD_accepted==1); % Keep only the trajectories that associated to the accepted MSD
Dapp = Dapp(Idx_MSD_accepted==1);
fprintf('\r')