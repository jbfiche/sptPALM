function [Reconstructed_Traj_Filtered, NTraj] = Filter_Trajectories_simulation_v1(Reconstructed_Traj, MaxBlink, MinTrajLength, MinNPoint)

Reconstructed_Traj_Filtered = {};

%% Start the analysis of all the trajectories to calculate the MSD and diffusion
%% distribution
%% ------------

for nTraj = 1 : size(Reconstructed_Traj,1)
    
    Traj = Reconstructed_Traj{nTraj};
    Frame_Traj_Filtered = [];
    
    Frame = Traj(1,:);
    X = Traj(2,:);
    Y = Traj(3,:);
    D = Traj(4,2:end);
    
    % Analyse the distances in order to detect the blinks. Remove the
    % points selected from the trajectory
    % -----------------------------------
    
    Idx_Blink = find(D==0);
    Frame(Idx_Blink+1) = [];
    X(Idx_Blink+1) = [];
    Y(Idx_Blink+1) = [];
    
    % Analyse the array "Frame" and split the trajectory into
    % sub-trajectories according to the value of "MaxBink". When the length
    % of the blinks is below "MaxBlink", the trajectory is kept as it is.
    % However, when the length of the blink is higher than "MaxBlink", the
    % trajectory is splited in two. The resulting sub-trajectories are kept
    % in "Frame_Traj_Filtered".
    % ------------------------
    
    First_Frame = 1;
    Npoint = 1;
    
    for nFrame = 2 : size(Frame,2)
        
        dt = Frame(nFrame)-Frame(nFrame-1);
        if dt > MaxBlink+1 && nFrame < size(Frame,2)
            
            Frame_Traj_Filtered = cat(1, Frame_Traj_Filtered, [First_Frame, nFrame-1, Npoint]);
            First_Frame = nFrame;
            Npoint = 1;
            
        elseif nFrame == size(Frame,2)
            
            Frame_Traj_Filtered = cat(1, Frame_Traj_Filtered, [First_Frame, nFrame, Npoint+1]);
            
        else
            
            Npoint = Npoint+1;
        end
    end
    
    % According to the values of "MinTrajLength" and "MinNPoint", the
    % trajectories are either validated and kept in
    % "Reconstructed_Traj_Filtered" or discarded.
    % -------------------------------------------
    
    for nSubTraj = 1 : size(Frame_Traj_Filtered,1)
        
        dt = Frame_Traj_Filtered(nSubTraj,2)-Frame_Traj_Filtered(nSubTraj,1);
        npoint = Frame_Traj_Filtered(nSubTraj,3);
        
        if dt >= MinTrajLength && npoint/(dt+1) >= MinNPoint
            
            FirstFrame = Frame_Traj_Filtered(nSubTraj,1);
            LastFrame = Frame_Traj_Filtered(nSubTraj,2);
            Reconstructed_Traj_Filtered = cat(1, Reconstructed_Traj_Filtered, [Frame(FirstFrame:LastFrame); X(FirstFrame:LastFrame); Y(FirstFrame:LastFrame)]);
        end
    end
    
    % Return the total number of trajectories selected
    % ------------------------------------------------
    
    NTraj = size(Reconstructed_Traj_Filtered,1);   
end