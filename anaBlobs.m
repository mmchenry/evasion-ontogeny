function B = anaBlobs(vPath,tPath,dPath,p,overwrite)
% Runs analysis of fish blobs


%% Parameters

% Initial visualizing predator acquisition
visSteps = 1;

% Minimum number of midline points to work with
minMid = 10;


%% Get path of data file, load data

% Suffix for filename
nameSuffix = 'mat';

% Load filenames for frames
a = dir([tPath  filesep '*.' nameSuffix]);

if isempty(a)
    warning('No thumbnail data files found');
    return
end

% Get indicies for video frames
for i = 1:length(a)
    
    % Read filename
    frNum = str2num(a(i).name(end-p.num_digit_frame-length(nameSuffix):...
                    end-length(nameSuffix)-1));
                
    % Check ordering    
    if (i>1) && (frNum ~= p.frNums(i-1)+1)
        error('Frame numbers not consecutive')
    end
    
    % Store
    p.frNums(i,1) = frNum;
    p.filename{i} = a(i).name;
end


%% Analyze blobs

% Look for data file
a3 = dir([dPath filesep 'blob data.mat']);


% If data file not present . . .
if isempty(a3) || overwrite
    
    % Update status
    disp('      Analyzing blob thumbnails . . .')   
    
    % Status update
    if visSteps
        % Figure window
        f = figure;
        set(f,'DoubleBuffer','on')
    else
        hW = waitbar(0,'1','Name','Analyzing thumbnails',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(hW,'canceling',0)
    end
    
    % Start timer
    tstart = tic;
    
    % Step though frames
    for i = 1:length(p.filename) %1107:length(p.filename)
        
        % Load current blob 'imBlob'
        load([tPath filesep p.filename{i}])
        
        % Get current filename & frame number
        B(i).filename = p.filename{i};
        B(i).framenum = p.frNums(i);
        
        % If blob is a not a nan . . .
        if ~isnan(max(imBlob(:)))
            
            % Return the body as a binary
            imBW = giveBigBlob(imBlob);
            
            % If first logged frame . . .
            if isempty(whos('B')) || ~isfield(B(i),'xHead') || (i-1)>length(B) || max(isnan(B(i-1).xMid))
                % Find first head and tail points
                [B(i).xHead,B(i).yHead,B(i).xTail,B(i).yTail] = findHeadTail(imBW);

            % Otherwise . . .
            else
                % Use prevous to find current head and tail
                [B(i).xHead,B(i).yHead,B(i).xTail,B(i).yTail] = findHeadTail(imBW,B);
                
            end
            
            % Find midline points
            [B(i).sMid,B(i).xMid,B(i).yMid] = findMid(imBW,B(i),minMid);
            
            % Check for enough points
            if isnan(B(i).sMid(1)) || (length(B(i).sMid) < minMid)
                B(i).sMid = nan;
                B(i).xMid = nan;
                B(i).yMid = nan;
                
            % If enough points . . .
            else
                
                % Find coordinate system and eyes
                [B(i).rost,B(i).origin,B(i).S,B(i).eye] = findEyes(imBlob,B(i));
                
                % Quantify curvature
                [B(i).sKappa,B(i).kappa] = calcTailCurve(B(i));
                
                % Visualize frames
                if visSteps
                    
                    % Read full frame
                    [imWholecl,cmap] = imread([vPath  filesep p.filename{i}(1:(end-3)) ...
                        p.nameSuffix]);
                    
                    warning off
                    % subplot(1,2,1)
                    imshow(imBlob,cmap)
                    hold on
                    %plot(yPerim,xPerim,'r-')
                    plot(B(i).xMid,B(i).yMid,'go')
                    h = plot([B(i).origin(1) B(i).rost(1)],...
                             [B(i).origin(2) B(i).rost(2)],'y-', ...
                              B(i).rost(1),B(i).rost(2),'yo');
                    set(h(2),'MarkerFaceColor','y')
                    set(h(1),'LineWidth',2)
                    set(h(2),'MarkerSize',7)
                    h = plot(B(i).eye.R(1),B(i).eye.R(2),'ro',...
                        B(i).eye.L(1),B(i).eye.L(2),'r+');
                    title(['Frame ' num2str(i)])
                    hold off
                    warning on
                    pause(0.01)
                    
                    % Or update status bar
                else
                    % Check for Cancel button press
                    if getappdata(hW,'canceling')
                        close(hW,'force')
                        return
                        
                        % Otherwise, update status
                    else
                        waitbar(i/length(p.frNums),hW,['Frames ' num2str(i) ' of ' ...
                            num2str(length(p.frNums))])
                    end
                end
            end
            
            clear imBlob2 imBW imBW2 LL props branchPts L maxB idx dists xEnd yEnd
            clear endPts
            
        else
            B(i).sMid = nan;
            B(i).xMid = nan;
            B(i).yMid = nan;
        end %isnan  
    end
    
    % Report duration
    telap = toc(tstart)/60;  
    disp(['     . . . completed in ' num2str(telap) ' min'])
    
    % Close waitbar, if there
    if ~visSteps
        close(hW,'force')
    end
    
    % Save data 'B'
    save([dPath filesep 'blob data.mat'],'B')
    beep
    beep
    
% If analysis already done . . .
else
    disp('      Blobs previously analyzed.')
end


function [ti,kappa] = calcTailCurve(B)

num_pts = length(B.xMid)*4;

x = B.xMid;
y = B.yMid;
t = B.sMid;
ti = linspace(t(1),t(end),num_pts)';

% Smoothing spline, calc coords, curvature (kappa)   
tol = 1.e-2;
sp = spaps(t',[x y]',tol);
dsp = fnder(sp);
dspt = fnval(dsp,ti);
ddspt = fnval(fnder(dsp),ti);
kappa = (abs(dspt(1,:).*ddspt(2,:)-dspt(2,:).*ddspt(1,:))./...
    (sum(dspt.^2)).^(3/2))';

if 0
    subplot(2,1,1)
    plot(B.xMid,B.yMid,'ko')
    axis equal
    hold on
    fnplt(sp)
    hold off
    
    subplot(2,1,2)
    plot(ti,kappa,'k-')
    pause(.1)
end


function [rost,origin,R,eye] = findEyes(im,B)

%TODO:save more data in each frame to allow for temporal averaging.  Also,
% it could be useful to share details of head

% Number of radials
num_rad = 200;

% Number of divisions of the head
num_div = 20;

% 
span_frac = 1/5;

% Radial positions
extent = 3*pi/4;
thetaR = [-extent :  extent/round(num_rad/2) : 0]';
thetaL = [extent  : -extent/round(num_rad/2) : 0]';

% Body positions that define margins of cranium and tail
pos_cran = 1/4;

% Indices for head and tail points
iHead = find((B.sMid./max(B.sMid))<pos_cran);

% Linear fit to head and tail points
cHx = polyfit(B.sMid(iHead),B.xMid(iHead),1);
cHy = polyfit(B.sMid(iHead),B.yMid(iHead),1);

iOrigin = round(length(iHead)*.5);
origin(1,1) = polyval(cHx,B.sMid(iOrigin));
origin(1,2) = polyval(cHy,B.sMid(iOrigin));
rost(1,1)   = polyval(cHx,B.sMid(iHead(1)));
rost(1,2)   = polyval(cHy,B.sMid(iHead(1)));

cran_len = hypot(origin(1)-rost(1),origin(2)-rost(2));

% Length of radials
len = cran_len.*1.2;

% Define rotation matrix for local coord system
R = local_system(origin,rost);

% Edge values on right side
[xEdgeR,yEdgeR,xLR,yLR] = findEdge(R,origin,im,thetaR,len);

% Edge values on left side
[xEdgeL,yEdgeL,xLL,yLL] = findEdge(R,origin,im,thetaL,len);

[x,iX] = sort(xLL);
w = mean([abs(yLL(iX)) abs(yLR(iX))],2);

[xPeaks,wPeaks,tmp] = polyScale(x,w,7);

%intvl = round(length(x)/num_div);

% c = polyfit(x,w,7);
% 
% Dc = polyder(c);DDc=polyder(Dc);
% 
% rts= roots(DDc);
% rts = rts((rts>min(x)) & (rts<max(x)));
% 
% %subplot(2,1,1); plot(x,w,'k',x,polyval(c,x),'r--',rts,polyval(c,rts),'ro');axis equal
% %subplot(2,1,2); plot(x,polyval(Dc,x),'r',rts,polyval(Dc,rts),'ro'); grid on


if length(xPeaks)<3 
    xEye = 0.3*cran_len;
else
    xEye = xPeaks(3) + .5*sqrt((xPeaks(2)-xPeaks(3))^2);
end

if (abs(xEye-max(x))/cran_len<0.5) || (abs(xEye-max(x))/cran_len>0.8)
    xEye = 0.3*cran_len;
end

 wEye = 0.8*max(w);

%wEye = polyval(c,xEye).*.8;

%rel_x = abs(xEye-max(x))/cran_len

% clear c Dc DDc
% 
% for i = 1:(num_div-1)
%     
%     idx = [1:round((length(x)*span_frac))] + (intvl*(i-1));
%     %idx = (i*intvl):(i+1)*intvl;
%     
%     if max(idx)>(length(x))
%        break
%     end
%    
%     % Fit for slope
%     c(i,:) = polyfit(x(idx),w(idx),1);
%     
%     xM(i,1) = mean(x(idx));
%     
%     if 0
%         subplot(2,1,1)
%        plot(x,w)
%        hold on
%         plot(x(idx),polyval(c(i,:),x(idx)),'k-') 
%         axis equal
%         hold off
%         
%         subplot(2,1,2)
%         plot(xM,c(:,1),'ko-')
%         grid on
%     end
%     
% end

% Find position of eyes in local FOR
%iEye = max([1 find(c(:,1)==max(c(:,1)),1,'first')+2]);
%tmp = abs(x-xM(iEye));
%xEye = x(tmp==min(tmp));
%wEye = 0.85*w(tmp==min(tmp));

% Determine global coordinates for the two eyes
[xEyeR,yEyeR] = local_to_global(origin,R,xEye,-wEye);
[xEyeL,yEyeL] = local_to_global(origin,R,xEye,wEye);

% Store values for the eye
eye.R = [xEyeR yEyeR];
eye.L = [xEyeL yEyeL];
eye.R_local = [xEye -wEye];
eye.L_local = [xEye wEye];

if 0
    figure
    subplot(2,1,1)
    plot(x,w,'k',xEye,wEye,'r+')
    
    subplot(2,1,2)
    plot(xM,c(:,1),'ko-')

end


%Show results
if 0
    subplot(2,1,1)
    plot(xLR,yLR,'b',xLL,yLL,'r')
    axis equal
    
    subplot(2,1,2)
    warning off
    imshow(im)
    warning on
    hold on
    plot(B.xMid(iHead),B.yMid(iHead),'yo')
    plot(B.xMid(iHead),polyval(cH,B.xMid(iHead)),'g-')
    plot(origin(1),origin(2),'g+')
    plot(xEdgeR,yEdgeR,'b-',xEdgeL,yEdgeL,'r-')
    plot(xEyeR,yEyeR,'y+')
    plot(xEyeL,yEyeL,'y+')
    hold off
    pause(0.1)
end


function [xEdgeG,yEdgeG,xLocal,yLocal] = findEdge(R,origin,im,thetas,len)
% Find edge in image at radials defined by theta

% Loop thru radials
for i = 1:length(thetas)
    
    % Radius values
    r = [0:ceil(len)]';
    
    % Coordinates for a radial (local FOR)
    x = r.*cos(thetas(i));
    y = r.*sin(thetas(i));
    
    % Transform radials into global FOR
    [xG,yG] = local_to_global(origin,R,x,y);
    
    % Get image values at radial points
    grayVal = interp2(double(im),xG,yG);
     
    % Eliminate nans
    iNoNan = ~isnan(grayVal);
    grayVal = grayVal(iNoNan);
    r = r(iNoNan);
    x = x(iNoNan);
    y = y(iNoNan);
    
    % Find radius of edge
    [tmp1,tmp2,rEdge] = polyScale(r,grayVal,9);
    
    % Coordinates of edge in local FOR
    xLocal(i,1) = rEdge*cos(thetas(i));
    yLocal(i,1) = rEdge*sin(thetas(i));
    
    % Coordinates of edge coords in global FOR
    [xEdgeG(i,1),yEdgeG(i,1)] = local_to_global(origin,R,...
                                                xLocal(i,1),yLocal(i,1));
    
    % Plot the selection of an edge
    if 0
        subplot(4,1,1:2)
        warning off
        imshow(im)
        warning on
        hold on
        %plot(B.xMid(iHead),B.yMid(iHead),'yo')
        %plot(B.xMid(iHead),polyval(c,B.xMid(iHead)),'g-')
        %plot(ptsT(:,1),ptsT(:,2),'g--')
        plot(x,y,'y-')
        plot(xEdgeG,yEdgeG,'ro');
        hold off
        
        subplot(3,1,3)
        plot(r,grayVal,'o')
        hold on
        plot(rEdge.*[1 1],ylim,'r-')
        hold off

    end
    
end


function [xPeaks,yPeaks,xHighPeak] = polyScale(x,y,ordr)
% Runs high-order polynomial fit with scaled data and returns peak values

% Set to 9th-order as default
if nargin < 3
    ordr = 9;
end

% Scale x-values
xS = (x-mean(x))./std(x);
 
% Polynomial fit
c = polyfit(xS,y,ordr); 
    
% Rates of change in gray values
Dc = polyder(c);
DDc = polyder(Dc);
    
% Get changes in the change in gray values
peak_vals = real(roots(DDc));
    
% Sort peak values
peak_vals = sort(peak_vals);
peak_vals = peak_vals(end:-1:1);

% Filter out-of-range values
peak_vals = peak_vals(peak_vals>min(xS) & peak_vals<max(xS));
    
% Get values of rate of change in gray at those points
peak_rate = polyval(Dc,peak_vals);
    
% Peak with biggest positive change
xHighPeakS = peak_vals(find(peak_rate==max(peak_rate),1,'first'));
     
% Transform peaks
xPeaks = peak_vals .* std(x) + mean(x);
yPeaks = polyval(c,peak_vals);
xHighPeak = xHighPeakS .* std(x) + mean(x);

if 0
    subplot(3,1,1)
    plot(x,y,'o')
    hold on
    plot(xHighPeak.*[1 1],ylim,'r-')
    plot(xPeaks,yPeaks,'r+')
    xlabel('x');ylabel('y')
    hold off
    
    subplot(3,1,2)
    plot(xS,y,'o')
    hold on
    plot(xS,polyval(c,xS),'r-')
    xlabel('xS');ylabel('y')
    hold off
    
    subplot(3,1,3)
    plot(xS,polyval(Dc,xS),'r-')
    hold on
    plot(peak_vals,polyval(Dc,peak_vals),'ro')
    plot(xHighPeakS,polyval(Dc,xHighPeakS),'r+')
    hold off
end


function R = local_system(origin,rost)

% Check dimensions
if size(origin,1)~=1 || size(origin,2)~=2 || size(rost,1)~=1 || size(rost,2)~=2 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = rost(1) - origin(1);
xaxis(1,2) = rost(2) - origin(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis(1:2); yaxis(1:2)];


function [xT,yT] = global_to_local(origin,R,x,y)
% Assumes columns vectors for coordinates

pts = [x y];

% Translate
pts(:,1) = pts(:,1) - origin(1);
pts(:,2) = pts(:,2) - origin(2);

% Rotate points
ptsT = [R * pts']';

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);


function [xT,yT] = local_to_global(origin,R,x,y)
% Assumes columns vectors for coordinates

pts = [x y];

% Rotate points
ptsT = [inv(R) * pts']';

% Translate global coordinates wrt origin
ptsT(:,1) = ptsT(:,1) + origin(1);
ptsT(:,2) = ptsT(:,2) + origin(2);

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);


function [xRT,yRT,zRT] = convert_eye(xR,yR,zR,verg_ang,fov)

% psi - forward tilt angle of an eye relative to the body
psi = pi/2 + verg_ang - fov/2;

% unit vector axes
yaxis = [-sin(psi) cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxis = [cos(psi) sin(psi) 0]./norm([cos(psi) sin(psi) 0]);
zaxis = [0 0 1];

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis' yaxis' zaxis'];

% Package points in a matrix
pts = [xR yR zR]';

% Rotate points
ptsT = [R * pts]';

% Extract points
xRT = ptsT(:,1);
yRT = ptsT(:,2);
zRT = ptsT(:,3);

if 1
   figure
   plot3(xRT,yRT,zRT,'.'); %axis equal
   xlabel('X'); ylabel('Y');zlabel('Z')
   view(2)  
end


function  [xHead,yHead,xTail,yTail] = findHeadTail(imBW,B)

% Skeletonize blob
skel = bwmorph(imBW,'skel',Inf);

% Find branch and endpoints
branchPts = bwmorph(skel, 'branchpoints');
endPts = bwmorph(skel, 'endpoints');

% Get coordinates for blob, branch & end points
[yEnd,xEnd]    = find(endPts);
[yBran,xBran]  = find(branchPts);
[ySkel,xSkel]  = find(skel);

% Assume first run, if 'B' is missing
if nargin<2 || isempty(B(end).xTail)
    
    % Compute distances between points
    for j = 1:length(xEnd)
        dists(:,j) = sqrt((xEnd-xEnd(j)).^2 + (yEnd-yEnd(j)).^2);
    end
    
    % Get indices of points furthest apart
    [iR,iC] = find(dists==max(dists(:)));
    
    % Theta & radius for ROI
    theta = linspace(0,2*pi,100);
    r = max(size(branchPts))/8;
    
    % Step thru each of those points
    for j = 1:length(iR)
        xROI   = r.*cos(theta) + xEnd(iR(j));
        yROI   = r.*sin(theta) + yEnd(iR(j));
        cBW    = roipoly(imBW,xROI,yROI) & imBW;
        
        % Number of white pixels in ROI
        bwVal(j,1) = sum(cBW(:));
        
        % Test visually
        if 0
            %figure;
            subplot(1,2,1)
            imshow(imBW);hold on; plot(xROI,yROI,'r-')
            subplot(1,2,2)
            imshow(cBW)
            title(['Sum = ' num2str(bwVal)])
        end
    end
    
    % Index for tail point
    iVal = find(bwVal==min(bwVal),1,'first');
    
    % Coordinate for tail point
    xTail = xEnd(iR(iVal));
    yTail = yEnd(iR(iVal));
    
    % Index for head point
    iVal = find(bwVal==max(bwVal),1,'first');
    
    % Coordinates for head point
    xHead = xEnd(iR(iVal));
    yHead = yEnd(iR(iVal));
    
    % Clear variables not needed below
    clear j iR iC theta r xROI yROI iVal bwVal cBW
    
% If not the first run    
else
    if length(B)==1 || isempty(B(end-2).xTail)
        xPredTail = B.xTail;
        yPredTail = B.yTail;
        xPredHead = B.xHead;
        yPredHead = B.yHead;
    else
        % Tail prediction
        xPredTail = B(end).xTail;
        yPredTail = B(end).yTail;
        
        % Prior displacement of head
        xDispHead = B(end).xHead - B(end-1).xHead;
        yDispHead = B(end).yHead - B(end-1).yHead;
        
        % Head prediction
        xPredHead = B(end).xHead + xDispHead;
        yPredHead = B(end).yHead + yDispHead;
        
    end
    
    % Distance between prior and tail point and endpoints
    xTail = repmat(xPredTail,length(xEnd),1);
    yTail = repmat(yPredTail,length(yEnd),1);
    dists = sqrt((xEnd-xPredTail).^2 + (yEnd-yPredTail).^2);
    
    % Index for new tail point as min distance to endpoint
    iVal = find(dists==min(dists),1,'first');
    xTail = xEnd(iVal);
    yTail = yEnd(iVal);
    
    % Distance bewteen predicted and current head point and endpoints
    xHead = repmat(xPredHead,length(xEnd),1);
    yHead = repmat(yPredHead,length(yEnd),1);
    dists = sqrt((xEnd-xPredHead).^2 + (yEnd-yPredHead).^2);
    
    % Index for new head point as min distance to endpoint
    iVal = find(dists==min(dists),1,'first');
    xHead = xEnd(iVal);
    yHead = yEnd(iVal);
end


function  [s,xMid,yMid] = findMid(imBW,B,minMid)
% Finds a smooth midline from skeletonized image

% Skeletonize blob
skel = bwmorph(imBW,'skel',Inf);

% Prune branches from skeleton
skel = prune_branches(skel,B.xHead,B.yHead,B.xTail,B.yTail);

% Number of pixels between midline points
mid_space = 5;

tol = 10;

% Store iniital image
skel_start = skel;

% Pad image
skel = pad_im(skel);

% Pad coordinates
xHead = B.xHead + 1;
yHead = B.yHead + 1;
xTail = B.xTail + 1;
yTail = B.yTail + 1;

% Define center cross im, peripheral im, diagonal im
im_cntr = [0 1 0;1 0 1; 0 1 0];
im_diag = [1 0 1; 0 0 0;1 0 1];

% Arc length
s = 0;

% Index
i = 2;

% Start at head
xMid = xHead;
yMid = yHead;

% Loop until done
while sum(skel(:))>mid_space-1
    
    % Remaining skeleton points
    [ySkel,xSkel]  = find(skel);
    
    % Distance between last midpoint and skel points
    xTmp = repmat(xMid(end),length(xSkel),1);
    yTmp = repmat(yMid(end),length(ySkel),1);
    dists = sqrt((xTmp-xSkel).^2 + (yTmp-ySkel).^2);
    
    dist_diff = abs(dists-mid_space);
     
    % Index for new tail point as min distance to endpoint
    iVal = find(dist_diff==min(dist_diff),1,'first');
    xMid(i,1) = xSkel(iVal);
    yMid(i,1) = ySkel(iVal);
    
    % Black out interval points
    idx = dists<=mid_space;
    skel(ySkel(idx),xSkel(idx)) = false;

    if 0
        imshow(skel)
        hold on
        plot(xMid(end),yMid(end),'r+')
        hold off
        pause(0.1);
    end
    
    % Check if at end
    if xMid(i)==xTail && yMid(i)==yTail
        break
    end
    
    % Iterate index
    i = i + 1;
end

% Remove effects of padding
xMid = xMid - 1;
yMid = yMid - 1;

% Position along curve
s = [0; cumsum(sqrt(diff(xMid).^2+diff(yMid).^2))];

% Check number of points
if length(xMid)<minMid
    % Store nans
    B.sMid = nan;
    B.xMid = nan;
    B.yMid = nan;
else 
    % Spline smoothed coordinates
    [xSp,xSm] = spaps(s,xMid,tol);
    [ySp,ySm] = spaps(s,yMid,tol);
    
    % Store coordinates
    B.sMid = s;
    B.xMid = xSm';
    B.yMid = ySm';
end

% Visual check
if 0
    warning off
    imshow(skel_start);hold on
    plot(xSm,ySm,'r-');hold off
    warning on
end


function imBW = giveBigBlob(imBlob)

% Guess for initial threshold value
tVal = min([0.95 graythresh(imBlob)+0.1]);

% Apply threshold
imBW    = ~im2bw(imBlob,tVal);

% Close gaps with dilation and erosion
se   = strel('disk',4,4);
imBW = imdilate(imBW,se);
imBW = imerode(imBW,se);

% Identify blobs
LL    = bwlabel(imBW);
props = regionprops(LL,'Centroid','Area');

% Get peripheral shapes
[bb,L] = bwboundaries(imBW,'noholes');

% Select blob with greatest periphery
maxB = 0;
idx = [];
for j = 1:length(bb)
    if length(bb{j}) > maxB
        maxB = length(bb{j});
        perim = bb{j};
        xPerim = perim(:,1);
        yPerim = perim(:,2);
        idx = j;
    end
end

% Define image as having only largest blob
imBW = LL==idx;


function skel = prune_branches(skel,xHead,yHead,xTail,yTail)
% Gets rid of spurs from skeleton

% Pad image
skel = pad_im(skel);

% Make clubs at the head and tail
ptsX = [-1 0 1 1 1 0 -1 -1];
ptsY = [-1 -1 -1 0 1 1 1 0];

% Get image of pixels around head and tail
imTail = skel(ptsY+yTail+1,ptsX+xTail+1);
imHead = skel(ptsY+yHead+1,ptsX+xHead+1);

% Create clubs to prevent erosion
skel(ptsY+yTail+1,ptsX+xTail+1) = true;
skel(ptsY+yHead+1,ptsX+xHead+1) = true;

% Initialize constants
chng = inf;
nLast = inf;

% Loop until no change in number of pixels
while chng>0
    % Remove spurs
    skel = bwmorph(skel,'spur');
    
    % Change in num of pixels
    chng = abs(sum(skel(:))-nLast);
    
    % Num of pixels in current frame
    nLast = sum(skel(:));
end

% Remove clubs at the head and tail
skel(ptsY+yTail+1,ptsX+xTail+1) = imTail;
skel(ptsY+yHead+1,ptsX+xHead+1) = imHead;

% Remove padding
skel = unpad_im(skel);


function im = pad_im(im)
% Pad images with single border of black pixels

% Top row
im = [false(1,size(im,2)); im];

% Bottom row
im = [im; false(1,size(im,2))];

% Left column
im = [false(size(im,1),1) im];

% Right column
im = [im false(size(im,1),1)];


function im = unpad_im(im)
% UnPad images with single border of black pixels

% Top row
im = im(2:end,:);

% Bottom row
im = im(1:end-1,:);

% Left column
im = im(:,2:end);

% Right column
im = im(:,1:end-1);   