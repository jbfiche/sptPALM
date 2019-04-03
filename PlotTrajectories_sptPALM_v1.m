function PlotTrajectories_sptPALM_v1(h)

N = round(get(h.Slider_SelectFrame, 'Value'));
PixelSize = str2double(get(h.PixelSize, 'String'));
SelectedTraj = find(h.Frame_Traj_MovieDisplay(:,1)<=N & h.Frame_Traj_MovieDisplay(:,2)>=N); % Select the trajectories according to the frame number
MinTrajLength = str2double(get(h.MinTrackLength, 'String'));

axes(h.MainAxes)
hold on

if isfield(h, 'MovieDisplayROI')
    ROI = h.MovieDisplayROI;
    Xo = ROI(1,1)-1;
    Yo = ROI(1,2)-1;
    Borders = [ROI(1,1), ROI(1,2) ; ROI(1,1)+ROI(1,3), ROI(1,2) ; ROI(1,1)+ROI(1,3), ROI(1,2)+ROI(1,4) ; ROI(1,1), ROI(1,2)+ROI(1,4) ; ROI(1,1), ROI(1,2)];
else
    Xo = 0;
    Yo = 0;
    ImSize = h.ImSize;
    Borders = [1, 1 ; 1 + ImSize(1), 1 ; 1 + ImSize(1), 1 + ImSize(2); 1, 1 + ImSize(2); 1, 1];
end

for ntraj = 1 : size(SelectedTraj,1)
    
    Traj = h.Reconstructed_Traj_MovieDisplay{SelectedTraj(ntraj)};
    TrahLength = size(Traj,2);
    
    if TrahLength >= MinTrajLength
        
%         TrajPoint_Id = find(Traj(1,:)==N);
%         TrajPoint_XY = Traj(2:3,TrajPoint_Id);
        Traj = Traj(2:3, Traj(1,:)<=N)/PixelSize;
        MeanTraj = mean(Traj,2);
        
        IN = inpolygon(MeanTraj(2), MeanTraj(1), Borders(:,1), Borders(:,2));
        
        if IN
            ColorVector = ColorCoding(TrahLength);
%             if ~isempty(TrajPoint_Id)
%                 centers = [TrajPoint_XY(2,1)+2-Xo, TrajPoint_XY(1,1)+2-Yo];
%                 radius = 3;
%                 viscircles(centers,radius,'LineWidth', 2, 'EdgeColor', ColorVector(TrajPoint_Id, :));
%                 %             plot(TrajPoint_XY(1,2)+2-Xo, TrajPoint_XY(1,1)+2-Yo, 'o', 'MarkerEdgeColor', ColorVector(TrajPoint_Id, :), 'MarkerFaceColor', ColorVector(TrajPoint_Id, :))
%             end
            
            if ~isempty(Traj)
                for n = 1 : size(Traj,2)-1
                    plot(Traj(2,n:n+1)+1-Xo, Traj(1,n:n+1)+1-Yo, 'Color', ColorVector(n, :), 'LineWidth',2);
                end
            end
        end
    end
end

hold off