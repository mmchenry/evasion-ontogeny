function skel = remove_branches(skel,B,E,xLast,yLast)
% Removes branches from a skeleton

jump_dist = 6;

% Skeleton where we subtract points as we progress thru image
skel_hist = skel;

% Define center cross im, peripheral im, diagonal im
im_cntr = [0 1 0;1 0 1; 0 1 0];
im_peri = [1 1 1;1 0 1;1 1 1];
im_diag = [1 0 1; 0 0 0;1 0 1];

% Step thru skeletal points, until all branches removed
while sum(B(:))~=0
    
    % Black out current point
    skel_hist(yLast,xLast) = false;
    
    xLast = max([xLast 2]);
    xLast = min([xLast size(skel_hist,2)-1]);
    
    yLast = max([yLast 2]);
    yLast = min([yLast size(skel_hist,1)-1]);
    
    % Define little roi
    imLocal = skel_hist((yLast-1):(yLast+1),(xLast-1):(xLast+1));
    
    % Local endpoints 
    E_local = E((yLast-1):(yLast+1),(xLast-1):(xLast+1));
    
    % Local branchpoints 
    B_local = B((yLast-1):(yLast+1),(xLast-1):(xLast+1));
    
    % Convert erronenous endpoints
    if sum(E_local(:))>1
        %E((yLast-1):(yLast+1),(xLast-1):(xLast+1)) = false;
    end
    
    % Coordinates for branches in local FOR
    [yB,xB] = find(imLocal & im_peri);
    
    % Coordinates for branches in global FOR
    xB = xLast + xB - 2;
    yB = yLast + yB - 2;
    
%     % Changes branch status, if less than 2 paths reside on branch
%     if B(yLast,xLast)==1 && length(xB)<2
%         B(yLast,xLast) = 0;
%     end

    
    % If this new point is a fat branchpoint --------------
    if B(yLast,xLast)==1 && sum(imLocal(:))>3
        
        % Coordinates of branch remaining points
        [yBran,xBran]  = find(B);
        
        % Distance between remaining branchpoints and current point
        xVal = repmat(xLast,length(xBran),1);
        yVal = repmat(yLast,length(yBran),1);
        dists = sqrt((xBran-xVal).^2 + (yBran-yVal).^2);
        
        [dists,idx] = sort(dists);
        xBran = xBran(idx);
        yBran = yBran(idx);
        
        % Index for new tail point as min distance to endpoint
        iVal = find(dists>jump_dist,1,'first');
        
        % If no branches outside of jump, then set all branches to false
        if isempty(iVal)
            B(:) = false;
        end
        
        xLast = xBran(iVal);
        yLast = yBran(iVal);
        
        clear xVal yVal dists xBran yBran
        
    % If this new point is at a simple branch --------------
    elseif B(yLast,xLast)==1
 
        % Initialize masks for each branch
        msk = false([size(skel) length(xB)]);
        
        % Initialize logicals for branch
        isB = ones(size(xB));
        
        % Initialize distances of branch
        Bdist = zeros(size(xB));
        
        % Initialize number of steps along branck
        nSteps = zeros(size(xB));
        
        % Skeleton where all points are substracted as they are considered
        skel_tmp = skel_hist;
        
        % Look thru branches
        for i = 1:length(xB)
            
            % Start current coordinate at branch point
            xD = xB(i);
            yD = yB(i);
            
            % Black out current branch
            %skel_hist(yB,xB) = false;
            
            k = 1;
            
            % While the current point is not on an edge or branch . . .
            while true
                
                % Black out current coord
                %skel_tmp(yD,xD) = false;
                
                % Create current coord in mask
                msk(yD,xD,i) = true;
                
                % Remove current coord from image
                skel_tmp(yD,xD) = false;
                
                % Check if a branch . . .
                if B(yD,xD)==1
                    isB(i) = 0;
                    Bdist(i) = inf;                        
                    nSteps(i) = k;
                    break
                    
                % Or an edge point . . .
                elseif E(yD,xD)==1
                    Bdist(i) = sqrt((xD-xLast)^2 + (yD-yLast)^2);
                    nSteps(i) = k;
                    break
                end
                
                % Define local roi
                imLocal = skel_tmp((yD-1):(yD+1),...
                                   (xD-1):(xD+1));
                
                % Visual check
                if 1
                    warning off
                    subplot(1,2,1)
                    imshow(skel_tmp);hold on
                    plot(xD,yD,'r+'); hold off
                    subplot(1,2,2)
                    imshow(imLocal)
                    pause(0.01);
                end
                
                % If next point on vertical or horizontal
                if max(max(imLocal & im_cntr))
                    
                    % Get next vertical or horizontal
                    [yNext,xNext] = find(imLocal & im_cntr,1,'first');
                    
                    % Otherwise, get next diagonal
                else
                    % Find next pixel in local FOR
                    [yNext,xNext] = find(imLocal & im_diag,1,'first');
                end
                
                % Check that we have a point
                if isempty(xNext)
                    error('Oops');
                end
                
                % Displace for next iteration
                xD = xD + xNext - 2;
                yD = yD + yNext - 2;
                
                k = k + 1;
                
                % Clear variables
                clear xNext yNext imLocal
            end
            
        end
        
        % Possible index values
        idx = 1:length(xB);
        
        if max(idx)>1 && max(isB==1)
            % Include only branch with end points
            idx = idx(isB==1);
            
            % Only prune if 2 or more steps
            if max(nSteps(idx))>2
                % Remove all branches from history
                for j = 1:length(idx)
                    skel_hist  = skel_hist - msk(:,:,idx(j));
                end
                
                % Select shortest branch
                idx = idx(find(Bdist(idx)==min(Bdist(idx)),1,'first'));
                
                % Check length
                if length(idx)~=1
                    error('Oops');
                end
                
                % Subtract short branch from skeletons
                skel = skel - msk(:,:,idx);
            end
            
        else
            ttt=3
        end

        % Un-set point as a branch point
        B(yLast,xLast) = false;
        
        % Clear variables
        clear skel_tmp xD yD msk idx isB Bdist xB yB xD yD
        
        %figure;imshow(skel);pause(.01);
        
        
    % If new point is not a branch --------------
    else
        
        % If next point on vertical or horizontal
        if max(max(imLocal & im_cntr))
            
            % Get next vertical or horizontal
            [yNext,xNext] = find(imLocal & im_cntr,1,'first');
            
            % Otherwise, just get next diagonal
        else
            % Find next pixel in local FOR
            [yNext,xNext] = find(imLocal & im_diag,1,'first');
        end
        
        % Check that we have a point
        if isempty(xNext)
            return
        end
        
        % Displace for next iteration
        xLast = xLast + xNext - 2;
        yLast = yLast + yNext - 2;
        
        % Erase point
        skel_hist(yLast,xLast) = false;
        
    end
    
    
    
    % Visual check
    if 0
        warning off
        subplot(1,2,1);
        imshow(skel);hold on;
        plot(xLast,yLast,'r+')
        subplot(1,2,2);
        imshow(skel_hist);hold on;
        plot(xLast,yLast,'r+')
        %title(num2str(i))
        warning on
        pause(0.01)
    end
    
end