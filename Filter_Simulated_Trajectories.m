function [Reconstructed_Traj_Filtered, NTraj] = Filter_Simulated_Trajectories(Reconstructed_Traj, MaxBlink, MinTrajLength, MinNPoint, MaxStepLength)

Reconstructed_Traj_Filtered = {};

for nTraj = 1 : size(Reconstructed_Traj,1)
    
    Traj = Reconstructed_Traj{nTraj};
    Frame_Traj_Filtered = [];
    
    Frame = Traj(:,1);
    X = Traj(:,2);
    Y = Traj(:,3);
    D = Traj(2:end,4);
    
    % Filter the trajectories that have single step distances above the
    % threshold
    % ---------
    
    if ~isempty(MaxStepLength)
        
        D_single_step = Traj(2:end,4)./(Frame(2:end)-Frame(1:end-1)); %Return the distance divided by the number of frames separating two successive detections
        
        while ~isempty( D_single_step(D_single_step>MaxStepLength)) && ~isempty(Frame)
            
            Idx = find(D_single_step>MaxStepLength,1);
            if Idx>1
                Frame(Idx+1) = [];
                X(Idx+1) = [];
                Y(Idx+1) = [];
                D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                D_single_step = D ./ (Frame(2:end)-Frame(1:end-1));
            else
                Frame(Idx) = [];
                X(Idx) = [];
                Y(Idx) = [];
                D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                D_single_step = D ./ (Frame(2:end)-Frame(1:end-1));
            end
        end
    end
    
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
    
    for nFrame = 2 : size(Frame,1)
        
       dt = Frame(nFrame)-Frame(nFrame-1);
       if dt > MaxBlink+1 && nFrame < size(Frame,1)
           
           Frame_Traj_Filtered = cat(1, Frame_Traj_Filtered, [First_Frame, nFrame-1, Npoint]);
           First_Frame = nFrame;
           Npoint = 1;
           
       elseif nFrame == size(Frame,1)
           
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
           Reconstructed_Traj_Filtered = cat(1, Reconstructed_Traj_Filtered, [Frame(FirstFrame:LastFrame), X(FirstFrame:LastFrame), Y(FirstFrame:LastFrame)]);
       end
    end
    
    % Return the total number of trajectories selected 
    % ------------------------------------------------
    
    NTraj = size(Reconstructed_Traj_Filtered,1);
    
end