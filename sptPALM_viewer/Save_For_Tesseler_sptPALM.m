function h = Save_For_Tesseler_sptPALM(h)

%% Create the txt files for SR-Tesseler clustering algorithm. The first file
%% is indicating the mean position of each trajectory (columns 1&2) as well
%% as the associated diffusion coefficient (column 3). The  fourth column is
%% associated to the first frame of the detection. This txt file is called
%% "Tesseler_diffusion.txt". 
%% The second txt file is keeping all the position and is associating to each
%% position the instant velocity in the file "Tesseler_Instant_Velocity.txt".
%% ========================================================================

clc
disp('Starting calculation for Tesseler files ...')

NTraj_Diff = size(h.Dapp,1);
AcquisitionTime = str2double(get(h.AcquisitionTime, 'String'));

Tesseler_diffusion = zeros(NTraj_Diff,4);
Tesseler_Instant_Velocity = cell(NTraj_Diff,1);

Trajectories = h.Reconstructed_Traj_MSD_accepted;
D = h.Dapp;

for ntraj = 1 : NTraj_Diff
    
    SingleTraj = Trajectories{ntraj};
    Tesseler_diffusion(ntraj,:) = [median(SingleTraj(2,:)), median(SingleTraj(3,:)), D(ntraj), SingleTraj(1,1)];
    
    Instant_velocity = zeros(size(SingleTraj,2)-1,4);
    for ndetection = 2 : size(SingleTraj,2)
        
        Xmean =  mean(SingleTraj(2,ndetection-1:ndetection));
        Ymean =  mean(SingleTraj(3,ndetection-1:ndetection));
        d = sqrt((SingleTraj(2,ndetection-1)-SingleTraj(2,ndetection))^2 + (SingleTraj(3,ndetection-1)-SingleTraj(3,ndetection))^2);
        Velocity = d/AcquisitionTime*1000;
        Instant_velocity(ndetection-1,:) = [Xmean, Ymean, Velocity, SingleTraj(1,ndetection)];
    end
    Tesseler_Instant_Velocity{ntraj} = Instant_velocity;
end

disp('Saving Tesseler_diffusion.txt ...')
T = array2table(Tesseler_diffusion, 'VariableNames', {'X_median', 'Y_median', 'D_apparent', 'Frame'});
writetable(T, 'Tesseler_diffusion.txt');

disp('Saving Tesseler_Instant_Velocity.txt ...')
Tesseler_Instant_Velocity = cell2mat(Tesseler_Instant_Velocity);
T = array2table(Tesseler_Instant_Velocity, 'VariableNames', {'X', 'Y', 'Instant_velocity_um_s', 'Frame'});
writetable(T, 'Tesseler_Instant_Velocity.txt');

disp('Tesseler files saved!')
