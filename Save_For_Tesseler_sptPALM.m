function h = Save_For_Tesseler_sptPALM(h)

%% Create the txt files for SR-Tesseler clustering algorithm. The first file
%% is indicating the mean position of each trajectory (columns 1&2) as well
%% as the associated diffusion coefficient (column 3). The  fourth column is
%% associated to the first frame of the detection. This txt file is called
%% "Tesseler_diffusion.txt". 
%% The second txt file is keeping all the position and is associating to each
%% position the instant velocity in the file "Tesseler_Instant_Velocity.txt".
%% ========================================================================

NTraj_Diff = size(h.Dapp,1);
AcquisitionTime = str2double(get(h.AcquisitionTime, 'String'));

Tesseler_diffusion = zeros(NTraj_Diff,4);
Tesseler_Instant_Velocity = [];

hwb = waitbar(0,'Saving the data for Tesseler');

for ntraj = 1 : NTraj_Diff
    
    waitbar(ntraj/NTraj_Diff*0.5)
    
   SingleTraj = h.Reconstructed_Traj_MSD_accepted{ntraj};
   Tesseler_diffusion(ntraj,:) = [median(SingleTraj(2,:)), median(SingleTraj(3,:)), h.Dapp(ntraj), SingleTraj(1,1)];
   
   for ndetection = 2 : size(SingleTraj,2)
       
       Xmean =  mean(SingleTraj(2,ndetection-1:ndetection));
       Ymean =  mean(SingleTraj(3,ndetection-1:ndetection));
       d = sqrt((SingleTraj(2,ndetection-1)-SingleTraj(2,ndetection))^2 + (SingleTraj(3,ndetection-1)-SingleTraj(3,ndetection))^2);
       Velocity = d/AcquisitionTime*1000;

       Tesseler_Instant_Velocity = cat(1, Tesseler_Instant_Velocity, [Xmean, Ymean, Velocity, SingleTraj(1,ndetection)]);
       
   end
end

T = array2table(Tesseler_diffusion, 'VariableNames', {'X_median', 'Y_median', 'D_apparent', 'Frame'});
writetable(T, 'Tesseler_diffusion.txt');
waitbar(3/4)

T = array2table(Tesseler_Instant_Velocity, 'VariableNames', {'X', 'Y', 'Instant_velocity_um_s', 'Frame'});
writetable(T, 'Tesseler_Instant_Velocity.txt');
waitbar(1)

close(hwb)
