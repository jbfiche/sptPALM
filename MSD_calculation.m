function [MSD_all,MSD_weight,Reconstructed_Traj_ROI] = MSD_calculation(Reconstructed_Traj_ROI,NTraj_ROI,MinNPointMSD,p);

hwaitbar = waitbar(0,'Calculating the MSD');
MSD_all = cell(NTraj_ROI,1);
MSD_weight = cell(NTraj_ROI,1);
Idx_TrajAcceptedMSD = zeros(NTraj_ROI,1);

for ntraj = 1 : NTraj_ROI
    
    MSD = [];
    Weight = [];
    
    waitbar(ntraj/NTraj_ROI);
    Traj = Reconstructed_Traj_ROI{ntraj};
    LagMax = (Traj(1,end) - Traj(1,1)); % Calculate the maximum lag time for this trajectory
    
    % Calculate the MSD. The lagtime goes from "1" to "LagTime-(MinNPointMSD-1)"
    % since we want, for each value of the MSD, an average over at least 
    % "MinNPointMSD" different values.
    % For each lagtime values, the MSD is kept only if there are at least
    % "MinNPointMSD" distances used for the calculation. Else, the value is
    % discarted.
    % --------- 
    
    for lag = 1 : LagMax-(MinNPointMSD-1)
        
        D_all = [];
        n = 1;
        nMSD = 0;
        MaxLagPoint = Traj(1, end) - lag; % Return the value of time after which there is no point left in the trajectory that could be separated by a lagtime equals to "lag"
        MaxNPoint = find(Traj(1,:)>MaxLagPoint,1); % Return a list of points that could not be used for the MSD calculation for lagtime = lag;
        
        while n < MaxNPoint && MaxNPoint >= MinNPointMSD
            
            Idx = find(Traj(1,:)==Traj(1,n)+lag, 1); % Check that two events have been detected with the rigth lag time
            if ~isempty(Idx)
                
                Xi = Traj(2,n);
                Xj = Traj(2,Idx);
                Yi = Traj(3,n);
                Yj = Traj(3,Idx);
                d = (Xi - Xj)^2 + (Yi - Yj)^2;
                D_all = cat(1, D_all, d);
                nMSD = nMSD+1;
            end
            n = n+1;
        end
        
        % If the number of segments for the calculation of the MSD is lower
        % than MinNPointMSD, the calculation is stopped
        % ---------------------------------------------
        
        if nMSD > MinNPointMSD
            MSD(lag) = mean(D_all);
            Weight(lag) = std(D_all);
        else
            break
        end
    end
    
    % Several trajectories (with blinks) can still gives MSD that have less
    % that p points, which can result in bug/error for the calculation of 
    % the apparent diffusion coefficient. In order to avoid this issue, a 
    % last check is performed here and all the MSD array with less than p 
    % points are discarded.
    % ---------------------
    
    if size(MSD,2)>=p
        MSD_all{ntraj} = MSD;
        MSD_weight{ntraj} = Weight;
        Idx_TrajAcceptedMSD(ntraj) = 1;
    end 
end

MSD_all = MSD_all(Idx_TrajAcceptedMSD==1);
MSD_weight = MSD_weight(Idx_TrajAcceptedMSD==1);
Reconstructed_Traj_ROI = Reconstructed_Traj_ROI(Idx_TrajAcceptedMSD==1);

close(hwaitbar);