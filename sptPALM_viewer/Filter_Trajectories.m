function [Reconstructed_Traj_Filtered, NTraj] = Filter_Trajectories(Reconstructed_Traj, MaxBlink, MinTrajLength, MinNPoint, MaxStepLength)

Reconstructed_Traj_Filtered = {};

%% Before analyzing the trajectories for the MSD and diffusion, all the
%% trajectories are quickly check for misconnections. Filter the trajectories
%% that have single step distances above a specific threshold defined by the
%% user. There are three different cases:
%%  1- the distance between the two first detections is above the threshold.
%%     In that case, the first point is automatically removed.
%%  2- same case for the last point
%%  3- if the misconnection is detected somewhere in between the first and
%%     last points, there is again two possibilities:
%%       a- the misconnection is only due to a single misdetection and the
%%          following points still belong to the same track. In that case,
%%          if the wrong detection is detected at frame 'n', the distances
%%          between 'n' and the points 'n+1' and 'n-1' should be larger than
%%          the mean distance separating two consecutive points of the
%%          trajectory. Therefore the distance between 'n-1' and 'n+1' should
%%          be smaller than the threshold. In that case the point 'n' is simply
%%          removed from the trajectory.
%%       b- the misconnection is between two different tracks. In that case,
%%          the distance between 'n' and 'n-1' and 'n-1' and 'n+1' should both
%%          be above the threshold. In that case, the initial track is split
%%          into two tracks.
%% -------------------------

NTraj = size(Reconstructed_Traj,1);
nTraj = 1;

if ~isnan(MaxStepLength)
    
    fprintf('Removing misconnexions ...     ')
    
    while nTraj <= NTraj
        
        fprintf('\b\b\b\b%03i%%', round(100*nTraj/NTraj))
        Traj = Reconstructed_Traj{nTraj};
        
        Frame = Traj(1,:);
        X = Traj(2,:);
        Y = Traj(3,:);
        D = Traj(4,2:end);
        
        Misconnections = find(D>MaxStepLength);
        N_misconnection = size(Misconnections,2);
        MeanDistance = mean(D(D<=MaxStepLength));
        StdDistance = std(D(D<=MaxStepLength));
        
        while N_misconnection>0
            
            if Misconnections(1) == 1
                Traj(:,1) = [];
                Frame(1) = [];
                X(1) = [];
                Y(1) = [];
                D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                
            elseif Misconnections(1) == size(D,2)
                Traj(:,end) = [];
                Frame(end) = [];
                X(end) = [];
                Y(end) = [];
                D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                
            else
                n = Misconnections(1);
                D2 = sqrt( (X(n)-X(n+2))^2 + (Y(n)-Y(n+2))^2 );
                
                if D2 < MeanDistance + 3*StdDistance
                    Traj(:,n+1) = [];
                    Frame(n+1) = [];
                    X(n+1) = [];
                    Y(n+1) = [];
                    D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                else
                    Reconstructed_Traj{end+1} = Traj(:,n+1:end);
                    Traj = Traj(:,1:n);
                    Reconstructed_Traj{nTraj} = Traj;
                    
                    Frame = Traj(1,:);
                    X = Traj(2,:);
                    Y = Traj(3,:);
                    D = sqrt( (X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2 );
                end
            end
            
            Misconnections = find(D>MaxStepLength);
            N_misconnection = size(Misconnections,2);
        end
        
        Reconstructed_Traj{nTraj} = Traj;
        NTraj = size(Reconstructed_Traj,1);
        nTraj = nTraj + 1;
    end
    
    fprintf('\r\n')
end

%% Start the analysis of all the trajectories to calculate the MSD and diffusion
%% distribution
%% ------------

for nTraj = 1 : size(Reconstructed_Traj,1)
    
    Traj = Reconstructed_Traj{nTraj};
    Frame_Traj_Filtered = [];
    
    try
    Frame = Traj(1,:);
    X = Traj(2,:);
    Y = Traj(3,:);
    D = Traj(4,2:end);
    catch
        disp('')
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