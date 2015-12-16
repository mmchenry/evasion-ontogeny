function  [blob, eye] = findLandmarks(blob,pBlob,im)
% Finds the midline and eye of a fish

% Return the body as a binary
%[blob.BW,imBWsmall] = giveBlobs(im);
                
% Number of radials for head points
%num_rad = 200;

% The extent of radial positions toward the posterior of the head
%extent = 3*pi/4;

% Realtive body position of the end of the cranium
pos_cran = 1/5;

% Maximum change in orientation between midline segments at the head
Dangl_max = pi/3;

% Get centroid of small blob
statsLarge = regionprops(blob.BW,'Area');

% Get centroid of nano blob
stats = regionprops(blob.BWsmall,'Centroid','Area');

% Circle radius (in pix), scaled to blob area
radius0 = round(sqrt(statsLarge(1).Area)/5);

%TODO: Use pBlob to use previous frame to find head

%% FIND TAIL POINT 

% Skeletonize blob
skel = bwmorph(blob.BW,'skel',Inf);

% Find branch and endpoints
%branchPts = bwmorph(skel, 'branchpoints');
endPts = bwmorph(skel, 'endpoints');

% Get coordinates for blob, branch & end points
[yEnd,xEnd]    = find(endPts);
%[yBran,xBran]  = find(branchPts);
%[ySkel,xSkel]  = find(skel);

% Distance between centroid and endpoints
xCent = repmat(stats.Centroid(1),length(xEnd),1);
yCent = repmat(stats.Centroid(2),length(yEnd),1);
dists = sqrt((xEnd-xCent).^2 + (yEnd-yCent).^2);

% Index for new tail point as max distance to endpoint
iVal = find(dists==max(dists),1,'first');
xTail = xEnd(iVal);
yTail = yEnd(iVal);

% Coordinates of image
[xAll,yAll] = meshgrid(1:size(blob.BW,2),1:size(blob.BW,1));


%% FIND MIDLINE POINTS 

% Current point (starts at the tail)
cX = xTail;
cY = yTail;

% Index that advances with the below loop
i = 1;

% Image of the blob to substract as the code moves anteriorly
imSketch = blob.BW;

% Boolean that switches the focus from the large to small blob
gosmall = 0;

radius = radius0;

% Loop that moves anteriorly from the tail, finding midline points
while true
   
    % Circle image
    circ = sqrt((xAll-cX).^2 + (yAll-cY).^2) <= radius & ...
           sqrt((xAll-cX).^2 + (yAll-cY).^2) > (radius-2);
       
    % Filled circle image   
    circFill = sqrt((xAll-cX).^2 + (yAll-cY).^2) <= radius;
    
    % If still using large image
    if gosmall==0
        
        % Create intersection image with small blob
        cIMsmall = circ & blob.BWsmall;

        % If there is any white . . .
        if max(cIMsmall(:))==1
            % Set sketch image to small blob
            imSketch = blob.BWsmall;
            
            % Change indicator variable
            gosmall = 1;
        end
    end
    
    % Image of intersection of blob and circle    
    cIM = circ & imSketch;
    
    % Check for break
    if max(cIM)==0
        break;
    end
    
    % Get centroid of small blob
    stats = regionprops(cIM,'MajorAxisLength','Centroid');   
     
    % If more than one blob
    if length(stats)>1
        maxLen = [0 0];
        
        % Choose the one with greater major AxisLength
        for j = 1:length(stats)
            if stats(j).MajorAxisLength>maxLen(2)
                maxLen = [j stats(j).MajorAxisLength];
                
                radius = max([stats(j).MajorAxisLength/2 radius0]);
            end
        end
    else
        maxLen = 1;
        radius = max([stats.MajorAxisLength/2 radius0]);
    end
    
    % Store coodinates of center point
    blob.xMid(i,1) = stats(maxLen(1)).Centroid(1);
    blob.yMid(i,1) = stats(maxLen(1)).Centroid(2);
    
    % Store arclength positon
    if i==1
        blob.sMid(i,1) = 0;
    else
        blob.sMid(i,1) = blob.sMid(i-1) + ...
                        hypot(blob.xMid(i)-blob.xMid(i-1),blob.yMid(i)-blob.yMid(i-1));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Visualize for debugging
    if 0
        % plot fish blob binary image with anterior sections subtracted
        subplot(1,2,1)
        imshow(imSketch,[],'InitialMagnification','fit')
        hold on; plot(blob.xMid,blob.yMid,'r+')
        hold off
        
        % plot intersection of blob and circle
        subplot(1,2,2)
        imshow(cIM,[],'InitialMagnification','fit')
        hold on; plot(blob.xMid(i),blob.yMid(i),'r+')
        hold off
%         pause
    end
    
    % Erase circle from imSketch
    imSketch = imSketch & ~circFill;
    
    % Define current point as new mid point
    cX = blob.xMid(i); cY = blob.yMid(i);
    
    % Advance index
    i = i + 1;
    
    % Clear for next iteration
    clear j maxLen circ circFill cIM stats
end

% Orientation of body segments
for i = 1:(length(blob.sMid)-1)
   angl(i,1) = atan2(blob.yMid(i+1)-blob.yMid(i),blob.xMid(i+1)-blob.xMid(i)); 
end
angl = unwrap(angl);

% If last two points have large angular changes, trim them
if abs(angl(end-2)-angl(end-1))>Dangl_max
    idx = 1:(length(blob.sMid)-2);
    
% If last point has large change, trim that
elseif abs(angl(end-1)-angl(end))>Dangl_max
    idx = 1:(length(blob.sMid)-1);

% Otherwise, keep all points
else
    idx = 1:length(blob.sMid);
end

% Redefine points
blob.sMid = blob.sMid(idx);
blob.xMid = blob.xMid(idx);
blob.yMid = blob.yMid(idx);

% Reverse direction of body position
blob.sMid = max(blob.sMid) - blob.sMid;
[blob.sMid,idx] = sort(blob.sMid);
blob.xMid = blob.xMid(idx);
blob.yMid = blob.yMid(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize progress to this point
if 0
    figure
    imshow(blob.im,[],'InitialMagnification','fit')
    hold on
    plot(blob.xMid,blob.yMid,'go-')
    hold off
%     pause
end


%% FIND HEAD POINTS 

% Get peripheral shapes
bound   = bwboundaries(blob.BW,'noholes');
yBound  = bound{1}(:,1);
xBound  = bound{1}(:,2);

% % Angular position of radials 
% thetaR = [-extent :  extent/round(num_rad/2) : 0]';
% thetaL = [extent  : -extent/round(num_rad/2) : 0]';

% Indices for head and tail points
iHead = find((blob.sMid./max(blob.sMid))<pos_cran);

% Linear fit to head and tail points
cHx = polyfit(blob.sMid(iHead),blob.xMid(iHead),1);
cHy = polyfit(blob.sMid(iHead),blob.yMid(iHead),1);

% Index for origin (point in center of head)
iOrigin = round(length(iHead)*.75);

% Coordinates for origin
origin(1,1) = polyval(cHx,blob.sMid(iOrigin));
origin(1,2) = polyval(cHy,blob.sMid(iOrigin));

% Coordinates for rostrum
rost(1,1)   = polyval(cHx,blob.sMid(iHead(1)));
rost(1,2)   = polyval(cHy,blob.sMid(iHead(1)));

% Length of head
cran_len = hypot(origin(1)-rost(1),origin(2)-rost(2));

% Length of radials
%len = cran_len.*2;

% Define rotation matrix for local coord system
R = local_system(origin,rost);

% Coordinates for edge of right side of head
%[xEdgeR,yEdgeR,xLR,yLR,rEdgeR] = findEdge(R,origin,im,thetaR,len);

% Coordinates for edge of left side of head
%[xEdgeL,yEdgeL,xLL,yLL,rEdgeL] = findEdge(R,origin,im,thetaL,len);

% Local coordinates
[xBoundL,yBoundL] = global_to_local(origin,R,xBound,yBound);

% Chop off body
idx = xBoundL>0;
xBoundL = xBoundL(idx);
yBoundL = yBoundL(idx);

% Remove tail, if present
idx = abs(yBoundL) < 2*cran_len;
xBoundL = xBoundL(idx);
yBoundL = yBoundL(idx);

% Redefine rostrum point in local FOR by max x-value
iRost = find(xBoundL==max(xBoundL),1,'first');
rostL = [xBoundL(iRost) yBoundL(iRost)];

% Redefine origin with y-value shifted to the wider side of the head
originLnew = [0 mean(yBoundL)];

% Transform rostrum and origin coordinates into global FOR
[rostNew(1),rostNew(2)]       = local_to_global(origin,R,rostL(1),rostL(2));
[originNew(1),originNew(2)]   = local_to_global(origin,R,...
                                                originLnew(1),originLnew(2));

% Define new first and second midline points
blob.xMid(1) = rostNew(1);
blob.yMid(1) = rostNew(2);
blob.xMid(2) = originNew(1);
blob.yMid(2) = originNew(2);

% Fill out remaining points
idx     = max([3 iOrigin+1]):length(blob.xMid);
blob.xMid  = [blob.xMid(1:2); blob.xMid(idx)];
blob.yMid  = [blob.yMid(1:2); blob.yMid(idx)];

% Revise arclength values
blob = rmfield(blob,'sMid');
blob.sMid(1) = 0;
for i = 2:length(blob.xMid)
    blob.sMid(i,1) = blob.sMid(i-1) + ...
                   hypot(blob.yMid(i)-blob.yMid(i-1),blob.xMid(i)-blob.xMid(i-1));
end


%% FIND EYES (currently commented out) 

% % Redefine rotation matrix for local coord system with new points
% R = local_system(originNew,rostNew);
% 
% % Local coordinates
% [xBoundL,yBoundL] = global_to_local(originNew,R,xBound,yBound);
% 
% % Chop off body
% idx = xBoundL>0;
% xBoundL = xBoundL(idx);
% yBoundL = yBoundL(idx);
% 
% % Remove tail, if present
% idx = abs(yBoundL) < 2*cran_len;
% xBoundL = xBoundL(idx);
% yBoundL = yBoundL(idx);
% 
% % Separate left and right sides
% idx = yBoundL>=0;
% xLL = xBoundL(idx);
% yLL = yBoundL(idx);
% xLR = xBoundL(~idx);
% yLR = yBoundL(~idx);
% 
% % New edge values on right side in local FOR
% %[xEdgeR,yEdgeR,xLR,yLR,rEdgeR] = findEdge(R,originNew,im,thetaR,len);
% 
% % New edge values on left side in local FOR
% %[xEdgeL,yEdgeL,xLL,yLL,rEdgeL] = findEdge(R,originNew,im,thetaL,len);
% 
% % Number of points in the top quartile to define eye position
% num_pts = round(length(xLR)/5);
% 
% % LEFT EYE ---
% 
% % Sort y-values by x-values
% [xLL,iL] = sort(xLL,1,'descend');
% yLL = yLL(iL);
% 
% % Cut-off of y-value at eye
% Lcut = .75*max(yLL);
% 
% % Index of first point
% iFirst = find(yLL>=Lcut,1,'first');
% 
% % Index of y-values to include
% iL = iFirst:min([(iFirst+num_pts) length(yLL)]);
% 
% % Mean x-position
% xEyeL_l = mean(xLL(iL));
% 
% % Mean (& adjusted) y-value
% yEyeL_l = .75 * mean(yLL(iL));
% 
% 
% % RIGHT EYE ---
% 
% % Sort y-values by x-values
% [xLR,iR] = sort(xLR,1,'descend');
% yLR = yLR(iR);
% 
% % Cut-off of y-value at eye
% Rcut = .75*max(abs(yLR));
% 
% % Index of first point
% iFirst = find(abs(yLR)>=Rcut,1,'first');
% 
% % Index of y-values to include
% iR = iFirst:min([(iFirst+num_pts) length(yLR)]);
% 
% % Mean x-position
% xEyeR_l = mean(xLR(iR));
% 
% % Mean (& adjusted) y-value
% yEyeR_l = .75 * mean(yLR(iR));
% 
% clear iR iL
% 
% % % Identify the rise in width at the right eye
% % [xLR,iR] = sort(xLR,1,'descend');
% % yLR = yLR(iR);
% % %Rcut = quantile(abs(yLR),.75);
% % Rcut = .75*max(abs(yLR));
% % iR = find(abs(yLR)>=Rcut,num_pts);
% % xEyeR_l = mean(xLR(iR:(iL+num_pts)));
% % yEyeR_l = .75 * mean(yLR(iR:(iL+num_pts)));
% 
% % (x,y) coordinates for eye positions in local coordinates
% xMean = mean([xEyeL_l xEyeR_l]);
% yMean = mean([yEyeL_l -yEyeR_l]);
% 
% % Determine global coordinates for the two eyes
% [xEyeR,yEyeR] = local_to_global(originNew,R,xMean,-yMean);
% [xEyeL,yEyeL] = local_to_global(originNew,R,xMean,yMean);
% 
% % Store values for the eyes (global coordinates)
% blob.xEye = [xEyeR xEyeL];
% blob.yEye = [yEyeR yEyeL];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 0
%     % Get peripheral shapes
%     [b,l] = bwboundaries(blob.BWsmall,'noholes');
%     
%     subplot(2,1,1)
%     imshow(blob.im,[],'InitialMagnification','fit')
%     hold on; 
%     plot(blob.xMid,blob.yMid,'r.')
%     plot(b{1}(:,2),b{1}(:,1),'-',rostNew(1),rostNew(2),'rs')
%     
%     plot(blob.xMid(iHead),blob.yMid(iHead),'ro')
%     %plot(polyval(cHx,blob.sMid(iHead)),polyval(cHy,blob.sMid(iHead)),'g-')
%     plot(origin(1),origin(2),'g+')
%     plot(blob.xEye,blob.yEye,'yo')
% 
% 
%     hold off
%     
%     subplot(2,1,2)
%     plot(xLL,yLL,'-',xLR,yLR,'-',xMean,yMean,'o',xMean,-yMean,'o',...
%         [min(xLL) max(xLL)],[1 1].*Lcut,'--',...
%         [min(xLL) max(xLL)],[-1 -1].*Rcut,'--')
%     axis equal
%     pause
% end

%% Find Eyes & Eye orientation

% indicator for debugging plots
sane = 0;


% Redefine rotation matrix for local coord system with new points
R = local_system(originNew,rostNew);

% % Store values for eye position in eye structure (global coordinates)
% eye.xEye = zeros(1,2);
% eye.yEye = zeros(1,2);

% DEFINE ROI ----

% rostrum point global coordinates
xRost = rostNew(1);
yRost = rostNew(2);

% origin point global coordinates
xOrigin = originNew(1);
yOrigin = originNew(2);

% rostrum in local coordinates
[xRostL,yRostL] = global_to_local(originNew,R,xRost,yRost);

% origin in local coordinates (already the origin...should be (0,0))
xOriginL = 0; 
yOriginL = 0;

% 50% the distance between origin and rostrum
% roiL.radius = 0.25*hypot(xRostL-xOriginL,yRostL-yOriginL);
roiL.radius = 0.5 * xRostL;

% x-coordinates for ROI around eyes in local FOR
roiL.x = [xRostL + roiL.radius; xOriginL - roiL.radius;...
    xOriginL - roiL.radius; xRostL + roiL.radius]; 

% y-coordinates for ROI around eyes in local FOR
roiL.y = [yRostL + roiL.radius; yOriginL + roiL.radius;...
    yOriginL - roiL.radius; yRostL - roiL.radius];

% ROI in global coordinates
[roiG.x,roiG.y] = local_to_global(originNew,R,roiL.x,roiL.y);

% EYE BLOBS --------

% initialize blob number indicator
blobNum = 1;

% initialize threshold value indicator
thresh = 1;

% crop original image to blob ROI
im = imcrop(im,blob.roi_blob);

% create a binary mask based on the selected ROI in global coordinates
maskEyes = roipoly(im,roiG.x,roiG.y);

% set initial threshold value
tVal = 74;

% set the stepsize for decreasing tVal
tStep = 2;

% adaptively find threshold value s.t. there are only two blobs
while blobNum && thresh
    
    % Threshold image and include only ROI
    imBW1  = ~im2bw(im,tVal/255) & maskEyes;
    
    % square structuring element with radius=2
    seEye = strel('square',2);
    
    % Perform a morphological close (dilation + erosion) operation on the image.
    imBW1 = imopen(imBW1,seEye);
    % figure, imshow(imBW,'InitialMagnification','fit')
    % title('Threshold + Open')
    
    % find convex hulls of objects in imBW, with 4-connected neighborhood
    imBW2 = bwconvhull(imBW1,'objects',4);
    % figure, imshow(CH_objects,'InitialMagnification','fit');
    % title('Convex hulls')
    
    % find connected components in binary image within ROI
    cc = bwconncomp(imBW2);
    
    % get size of blobs
    ccProps = regionprops(cc,'Area');
    
    
    % if there are only two objects & their area is similar ...
    if cc.NumObjects==2 && (abs(ccProps(1).Area-ccProps(2).Area))<45
        
        % turn off indicator to exit loop
        blobNum = 0;
        
        % otherwise...
    else
        
        % decrease the threshold value
        tVal = tVal - tStep;
    end
    
    % check that threshold value doesn't get too low
    if tVal < 55
        % set threshold indicator to exit loop
        thresh = 0;
    end
        
    
end

% call regionprops with list of desired measurements 
eyeProps = regionprops(imBW2, 'Orientation', 'Centroid', 'Area');

% LEFT or RIGHT EYE ------

% y-position of first eye in local coordinates
[~,eye1Pos] = global_to_local(originNew,R,...
    eyeProps(1).Centroid(1),eyeProps(1).Centroid(2));

% NOTE: flip sign on eye1Pos for correct L/R eye 
% This is because image origin is top-left, and y-val increases as you go
% down the rows

% check if the first eye is the left or the right eye 
if -eye1Pos<0
    
    % save left eye position
    xEyeL = eyeProps(2).Centroid(1);
    yEyeL = eyeProps(2).Centroid(2);
    
    % save right eye position 
    xEyeR = eyeProps(1).Centroid(1);
    yEyeR = eyeProps(1).Centroid(2);
    
    % save orientation of left eye
    anglL = eyeProps(2).Orientation;
    
    % save orientation of right eye
    anglR = eyeProps(1).Orientation;
    
    % logical to identify which is right eye 
    % 1 corresponds to first element of eyeProps
    rightEye = 1;  
else
    % save left eye position 
    xEyeL = eyeProps(1).Centroid(1);
    yEyeL = eyeProps(1).Centroid(2);
    
    % save right eye position 
    xEyeR = eyeProps(2).Centroid(1);
    yEyeR = eyeProps(2).Centroid(2);
    
    % save orientation of left eye
    anglL = eyeProps(1).Orientation*pi/180;
    
    % save orientation of right eye
    anglR = eyeProps(2).Orientation*pi/180;  
    
    % logical to identify which is right eye
    % 1 corresponds to first element of eyeProps
    rightEye = 0;
end

% Store values for eye position (global coordinates)
blob.xEye = [xEyeR xEyeL];
blob.yEye = [yEyeR yEyeL];

% Store values for eye position in eye structure (global coordinates)
eye.xEye = [xEyeR xEyeL];
eye.yEye = [yEyeR yEyeL];

% store eye orientation (from ellipse fit)
eye.thetaR = anglR;
eye.thetaL = anglL; 

% body midline angle w.r.t. global x-axis (alternative calculation)
eye.anglBody = atan2(eye.xEye(1)-eye.xEye(2),eye.yEye(1)-eye.yEye(2));

% slope of line defining body midline
eye.m1 = (eye.xEye(1)-eye.xEye(2))/(eye.yEye(1)-eye.yEye(2));

%------------------------- Sanity Check ---------------------------------%
if sane
    close all
    
    figOrig = figure;
    % plot the binary mask for eye region
    subplot(2,1,1), imshow(imBW2)
    hold on; plot(xOrigin,yOrigin,'ro',xRost,yRost,'yo')
    axs1 = gca;

    % plot eyes 
    plot(axs1,blob.xEye,blob.yEye,'ow');
    
    % plot original image, rostrum & head points
    subplot(2,1,2), imshow(im)
    hold on; plot(xOrigin,yOrigin,'ro',xRost,yRost,'yo')
    axs2 = gca;

    % plot eyes
    plot(axs2,blob.xEye,blob.yEye,'ok');
% pause
end

% EYE ORIENTATION -------

% preallocate boundary points vectors
xPnts = zeros(2,2);
yPnts = zeros(2,2);

% preallocate boundary points vectors (local FOR)
xPntsL = zeros(2,2);
yPntsL = zeros(2,2);

% preallocate eye orientation vector
eyePhi = zeros(1,2);

% preallocate eye orientation vector (local FOR)
eyePhiL = zeros(1,2);

% preallocate eye orientation vector (local FOR, alternative calc.)
eyePhiL_alt = zeros(1,2);

% find object boundaries (points along the perimeter of eye blobs) 
eyeBound = bwboundaries(imBW2,'noholes'); 

for j=1:length(eyeBound)
    
    % extract x-coordinates of boundary
    xBoundary = eyeBound{j}(:,2);
    
    % extract y-coordinates of boundary
    yBoundary = eyeBound{j}(:,1);
    
    % compute distances between all points
    ptDists = dist([xBoundary'; yBoundary']);
    
    % find the max val for each column of ptDists & their index
    [Y,I]=max(ptDists);
    
    % find the largest of all maxima and its index in Y (one point)
    [~,ind1]=max(Y);
    
    % find the index of the other point
    ind2 = I(ind1);
    
    % store x-coordinates of points at far end of blobs
    xPnts(j,:) = [eyeBound{j}(ind1,2), eyeBound{j}(ind2,2)];
    
    % store y-coordinates of points at far end of blobs
    yPnts(j,:) = [eyeBound{j}(ind1,1), eyeBound{j}(ind2,1)];
    
    % x and y-coords in local FOR of same points
    [tmp_xPntsL, tmp_yPntsL] = global_to_local(originNew,R,...
        xPnts(j,:)',yPnts(j,:)');
    
    xPntsL(j,:) = tmp_xPntsL';
    
    yPntsL(j,:) = tmp_yPntsL';
  
    
    % NOTE: need to change order of difference in y-vals for correct
    % eye orientation angle
    
    % compute eye orientation (w.r.t. body axis) in radians
    eyePhiL(j) = atan((yPntsL(j,2)-yPntsL(j,1))/(xPntsL(j,1)-xPntsL(j,2)));
    
    % compute eye orientation (w.r.t. global x-axis) in radians 
%     eyePhi(j) = atan((yPnts(j,2)-yPnts(j,1))/(xPnts(j,1)-xPnts(j,2)));

    % slope of line defining eye orientation (global coordinates)
    m2 = (yPnts(j,2)-yPnts(j,1))/(xPnts(j,1)-xPnts(j,2));
    
    % alternative eye orientation calculation
%     eyePhiL_alt(j) = abs(eye.anglBody - eyePhi(j));
    eyePhiL_alt(j) = atan2((eye.m1-m2),(1+eye.m1*m2));
    
    %------------------------------- SANITY CHECK #2 --------------------%
    
    if sane
        % bring figure with original image to front
        figure(figOrig);
        
        % draw line between points
        line(xPnts(j,:),yPnts(j,:),'LineWidth',2);
        
        % include text of eye angle
        text(blob.xEye(1,j),blob.yEye(1,j),...
            num2str(eyePhiL_alt(j)*180/pi),'FontSize',12,'Color','r');
        
        pause
    end
    
%     plot(axs2,xBoundary,yBoundary,'m')
end

% if rightEye is true (i.e. if first element of eyeProps is right eye)...
if rightEye
    
    % do nothing...
else
    
    % flip order of xPnts
    xPnts = flip(xPnts);
    
    % flip order of yPnts
    yPnts = flip(yPnts);
    
    % flip order of eyePhi
    eyePhi = flip(eyePhi);
    
    % flip order of eyePhiL
    eyePhiL = flip(eyePhiL);
    
    % flip order of eyePhiL_alt
    eyePhiL_alt = flip(eyePhiL_alt);
end

% store original (cropped) image
blob.imOrig = im;

% store original (cropped) image
eye.im = im;

% store points used for eye orientation calculation
eye.xPnts = xPnts;
eye.yPnts = yPnts;

% store points used for eye orientation calculation (local FOR)
eye.xPntsL = xPntsL;
eye.yPntsL = yPntsL;

% store eye orientation (first element is right eye)
eye.Phi = eyePhi;

% store eye orientation (first element is right eye), local FOR
eye.PhiL = eyePhiL;

% store alternate eye orientation (first element is right eye), local FOR
eye.PhiL_alt = eyePhiL_alt;

% store rotation matrix
eye.R = R;

% store origin
eye.origin = originNew;

% store threshold value used
eye.tVal = tVal;




% function [imBW,imBWsmall] = giveBlobs(imBlob)
% 
% % Guess for initial threshold value
% tVal = min([0.95 graythresh(imBlob)+0.1]);
% 
% % Find largest blob
% [imBW,numBlob] = returnBlob(imBlob,tVal);
% 
% % Create tmp image (uses blob as mask)
% tmp = imBlob.*0+2.^16;
% tmp(imBW(:)) = imBlob(imBW(:));
% 
% % Starting low threshold value 
% tVal2 = tVal.*.75;
% 
% % Increment to increase threshold
% tIncr = .05;
% 
% % Loop to create image at lowest threshold that finds single blob
% while true
%  
%     % Find blob
%     [imBWsmall,numBlob] = returnBlob(tmp,tVal2);
%     
%     % Check number of blobs
%     if numBlob==1
%         break
%     else
%         tVal2 = tVal2 + tIncr;
%     end
%     
%     % Check for too high a threshold
%     if tVal2 >= tVal
%         tVal2 = tVal;
%         imBWsmall = imBW;
%         break
%     end
% end

% Find nano blob
%[imBWnano,numBlob] = returnBlob(tmp,0.6*tVal2);

%tmpBW = im2bw(tmp,tVal2);

% Close gaps with dilation and erosion
%se   = strel('disk',6);
%tmpBW = imdilate(tmpBW,se);
%tmpBW = imerode(tmpBW,se);


function [imBW,numBlob] = returnBlob(im,tVal)
% Returns binary image of blobs

% Apply threshold
imBW    = ~im2bw(im,tVal);

% Close gaps with dilation and erosion
se   = strel('disk',6);
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

% Return number of blobs
numBlob = max(LL(:));

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
% Assumes column vectors for coordinates

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
ptsT = (R \ pts')';
% ptsT = [inv(R) * pts']';

% Translate global coordinates wrt origin
ptsT(:,1) = ptsT(:,1) + origin(1);
ptsT(:,2) = ptsT(:,2) + origin(2);

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);