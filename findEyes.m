function eyeBlobs = findEyes(imPath)

close all

if nargin < 1
% imPath = '/Volumes/Backup/ZF_visuomotor/Raw video/2015-11-16/S01';
imPath = '/Volumes/VisualPred/ZF_visuomotor/Raw video/2015-11-16/S01';
imPath = [imPath filesep 'exp01_0006.jpg'];
end

% read in the image
imA = imread(imPath);

% initialize blob number indicator
blobNum = 1;

% initialize threshold value indicator
thresh = 1;

% use CLAHE for histogram equalization
% imA = adapthisteq(imA,'ClipLimit',0.03,'Distribution','uniform');

% Select region of interest from cropped image 
warning off
imshow(imA)
warning on
title('Choose ROI points, include both eyes.')

% Interactively find ROI
h = impoly;
roi_poly = wait(h);

% Store results
tmp = getPosition(h);
roi.x = tmp(:,1);
roi.y = tmp(:,2);

delete(h), close all;

% create a binary mask based on the selected ROI (select head region)
maskEyes = roipoly(imA,roi.x,roi.y);

% set intial threshold value (will be divided by 255 later)
tVal = 74;
% TO DO: figure out better way to find the threshold value

% set the stepsize for decreasing tVal
tStep = 2;

while blobNum && thresh
    
    % Threshold image and include only ROI
    imBW1  = ~im2bw(imA,tVal/255) & maskEyes;
    
    % square structuring element with radius=2
    se = strel('square',2);
    
    % Perform a morphological open (erosion + dilation) operation on the image.
    imBW1 = imopen(imBW1,se);
    % figure, imshow(imBW1,'InitialMagnification','fit')
    % title('Threshold + Open')
    
    % find convex hulls of objects in imBW, with 4-connected neighborhood
    imBW2 = bwconvhull(imBW1,'objects',4);
    % figure, imshow(imBW2,'InitialMagnification','fit');
    % title('Convex hulls')

    % find connected components in binary image within ROI
    cc = bwconncomp(imBW2);
    
    ccProps = regionprops(cc,'Area','Orientation');
    
    % if there are only two objects ...
    if cc.NumObjects==2 && (abs(ccProps(1).Area-ccProps(2).Area))<45
        
        % turn off indicator to exit loop
        blobNum = 0;
        
        % otherwise...
    else
        
        % decrease the threshold value
        tVal = tVal - tStep;
    end
    
        % check that threshold value doesn't get too low
    if tVal < 50
        % set threshold indicator
        thresh = 0;
    end
    
end
tVal

%% Find eye blobs and fit ellipses

% call regionprops to get area stats on connected components
% ccArea = regionprops(cc,'Area');

% call regionprops with list of desired measurements for fitting ellipses
eyeElips = regionprops(imBW2, 'Orientation', 'MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid', 'Area', 'ConvexImage');

% % find the largest blobs 
% idx = find([eyeElips.Area] > 45);
% 
% % define binary mask of eye blobs
% imBW2 = ismember(labelmatrix(cc), idx);
% 
% % redefine eyeElips so that only the large blobs are included
% eyeElips = eyeElips(idx);

%% Plots for Sanity Check
figure;
axs2 = axes;
imshow(imBW2,'InitialMagnification','fit');
hold on
title('Eye Blobs')

figOrig = figure;
axs1 = axes;
imshow(imA,'InitialMagnification','fit');
hold on
title('Original + Eye Outlines + Eye Orientation')

%% Find object boundaries and look for points with largest distance 

% turn on to use alternative method 1 
if 1
    
    % find oject boundaries
    B = bwboundaries(imBW2,'noholes');
    
    % preallocate eye orientation vector
    eyePhi = zeros(1,2);
    
    for j=1:length(B)
        
        % extract x-coordinates of boundary
        xBoundary = B{j}(:,2);
        
        % extract y-coordinates of boundary
        yBoundary = B{j}(:,1);
        
        % Plot boundary points over original image
%         plot(axs1,xBoundary,yBoundary,'m','LineWidth',2);
        
        % compute distances between all points
        ptDists = dist([xBoundary'; yBoundary']);
        
        % find the max val for each column of ptDists & their index 
        [Y,I]=max(ptDists);
        
        % find the largest of all maxima and its index in Y (one point)
        [~,ind1]=max(Y);
        
        % find the index of the other point
        ind2 = I(ind1);
        
        % store x-coordinates of points at far end of blobs
        xPnts = [B{j}(ind1,2), B{j}(ind2,2)];
        
        % store y-coordinates of points at far end of blobs
        yPnts = [B{j}(ind1,1), B{j}(ind2,1)];
        
        % plot centroid
        plot(axs2,eyeElips(j).Centroid(1),...
            eyeElips(j).Centroid(2),'ko');
        
        % bring figure with original image to front
        figure(figOrig);
        
        % draw line between points 
        line(xPnts,yPnts,'color','c','LineWidth',2);
        
        % NOTE: need to change order of difference in y-vals for correct
        % eye orientation angle 
        
        % compute eye orientation (w.r.t. global x-axis) in radians
%         eyePhi(j) = atan2(yPnts(2)-yPnts(1),xPnts(1)-xPnts(2));
        eyePhi(j) = atan((yPnts(2)-yPnts(1))/(xPnts(1)-xPnts(2)));
        
%         plot(axs2,xBoundary,yBoundary,'m','LineWidth',2);
    end

end

% flip for testing only
% eyePhi=flip(eyePhi);

%% parametrize ellipses and plot over original image 

phi = 0:pi/64:2*pi;
cosphi = cos(phi);
sinphi = sin(phi);

for k = 1:length(eyeElips)
    
    % center of ellipse
    xbar = eyeElips(k).Centroid(1);
    ybar = eyeElips(k).Centroid(2);
    
    % radius of major/minor axes 
    a = eyeElips(k).MajorAxisLength/2;
    b = (eyeElips(k).MinorAxisLength/2);

    % TO DO: Figure out if I need to adjust this angle in some way...
    % convert orientation angle to radians
    theta = pi*eyeElips(k).Orientation/180;
    
    % combine both computed values, see if this works better
    theta2 = mean([theta,eyePhi(k)]); 
    
    % rotation matrix for rotating axes (transpose of active rotation
    % matrix)
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];
     
    % define points of ellipse, centered at origin w/ orientation=0 (deg)
    xy = [a*cosphi; b*sinphi];
    
    % rotate axes about angle theta
    xy = R*xy;
    
    % define points of ellipse centered at (xbar, ybar)
    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    
    % define points to plot majorAxis
%     lineMajor = R*[xbar + a,  xbar - a; ybar, ybar];
    lineMajor = [x(1),x(65); y(1),y(65)];
    
    % plot ellipses over original image and eye blobs 
    plot(axs1,x,y,'r','LineWidth',2);
    
    % draw major axis line
    line(lineMajor(1,:),lineMajor(2,:),'color',[0 0.5 0],'Linewidth',2);
    
    plot(axs2,x,y,'r','LineWidth',2);
%     line(lineMajor(:,1),lineMajor(:,2),'color',[0 0.5 0],'Linewidth',2);
end
% hold off


%% skeltonize threshold image and find endpoints

% turn on to use alternative method 2
if false 
    
    % skeletonize threshold image
    skel = bwmorph(imBW2,'skel');
    
    % find extrema and area of eye blob skeleton
    s = regionprops(skel,'Area','Extrema');
    
    % find points in skeleton with largest distance
    for k=1:length(s)

        % extract x-coordinates of boundary
%         xPnts = s2(k).Extrema(:,1);
        
        % extract y-coordinates of boundary
%         yPnts = s2(k).Extrema(:,2);
        
        % Plot all extrema points over original image
%         plot(axs1,xPnts,yPnts,'m','LineWidth',2);
        
        % compute distances between all extrema points
%         ptDists = dist([xPnts'; yPnts']);
        ptDists = dist(s(k).Extrema');
        
        % find the max val for each column of ptDists & their index 
        [Y,I]=max(ptDists);
        
        % find the largest of all maxima and its index in Y (one point)
        [~,ind1]=max(Y);
        
        % find the index of the other point
        ind2 = I(ind1);
        
        % store x-coordinates of points at far end of blobs
%         xPnts = [B{j}(ind1,2), B{j}(ind2,2)];
        xPnts = [s(k).Extrema(ind1,1), s(k).Extrema(ind2,1)];
        
        % store y-coordinates of points at far end of blobs
%         yPnts = [B{j}(ind1,1), B{j}(ind2,1)];
        yPnts = [s(k).Extrema(ind1,2), s(k).Extrema(ind2,2)];
        
        % bring figure with original image to front
        figure(figOrig);
        
        % draw line between points 
        line(xPnts,yPnts,'color','m','LineWidth',2);
    end
end
    
    
    

