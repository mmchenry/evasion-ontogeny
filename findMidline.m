function  blob = findMidline(p,blob,pBlob,im)
% Finds the midline and eye of a fish

% Return the body as a binary
%[blob.BW,imBWsmall] = giveBlobs(im);
                
% Number of radials for head points
%num_rad = 200;

% The extent of radial positions toward the posterior of the head
%extent = 3*pi/4;

% Relative body position of the end of the cranium
pos_cran = 1/5;

% Maximum change in orientation between midline segments at the head
Dangl_max = pi/3;

% Get centroid of small blob
statsLarge = regionprops(blob.BW,'Area','Centroid');

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

clear skel endPts xCent yCent dists iVal


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
    %figure
    imshow(blob.im,[],'InitialMagnification','fit')
    hold on
    plot(blob.xMid,blob.yMid,'go-')
    hold off
%     pause
end

clear idx angl radius radius0 i gosmall imSketch cIM cX cY circ circFill
clear statsLarge xTail yTail Dangl_max cIMsmall pBlob xAll yAll xEnd yEnd


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

% Redefine origin in local FOR with y-value shifted to the wider side of the head
originLnew = [0 mean(yBoundL)];

% Transform rostrum and origin coordinates into global FOR
[rostNew(1),rostNew(2)]       = local_to_global(origin,R,rostL(1),rostL(2));
[originNew(1),originNew(2)]   = local_to_global(origin,R,...
                                                originLnew(1),originLnew(2));

% Redefine the rotation matrix with the new rostrum and origin
R = local_system(originNew,rostNew);

% Create coordinates for head pixel values
[iXleft,iYleft]   = meshgrid(-1.5*cran_len:1.5*cran_len,-1.5*cran_len:0);
[iXright,iYright] = meshgrid(-1.5*cran_len:1.5*cran_len,0:1.5*cran_len);

% Transform coordinates into global FOR
[iXleft,iYleft]   = local_to_global(originNew,R,iXleft,iYleft);
[iXright,iYright] = local_to_global(originNew,R,iXright,iYright);

% Capture pixel values at those coordinates for the 
imLeft  = cast(interp2(double(blob.im),iXleft,iYleft),class(blob.im));
imRight = cast(interp2(double(blob.im),iXright,iYright),class(blob.im));

% % Store blob stuff
% blob.R = R;
% blob.origin  = originNew;
% blob.imLeft  = imLeft;
% blob.imRight = imRight;

% Define new first and second midline points
blob.xMid(1) = rostNew(1);
blob.yMid(1) = rostNew(2);
blob.xMid(2) = originNew(1);
blob.yMid(2) = originNew(2);

% Fill out remaining points
idx        = max([3 iOrigin+1]):length(blob.xMid);
blob.xMid  = [blob.xMid(1:2); blob.xMid(idx)];
blob.yMid  = [blob.yMid(1:2); blob.yMid(idx)];

% Revise arclength values
blob = rmfield(blob,'sMid');
blob.sMid(1) = 0;
for i = 2:length(blob.xMid)
    blob.sMid(i,1) = blob.sMid(i-1) + ...
                   hypot(blob.yMid(i)-blob.yMid(i-1),blob.xMid(i)-blob.xMid(i-1));
end

clear idx i iRost rostL rostNew xBound yBound pos_cran R bound cHz cHy origin
clear rost originNew xBoundL yBoundL cHx cran_len iHead iOrigin

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

for i = 1:size(x,2)
    
    pts = [x(:,i) y(:,i)];
    
    % Rotate points
    ptsT = (R \ pts')';
    
    % Translate global coordinates wrt origin
    ptsT(:,1) = ptsT(:,1) + origin(1);
    ptsT(:,2) = ptsT(:,2) + origin(2);
    
    % Extract columns of points
    xT(:,i) = ptsT(:,1);
    yT(:,i) = ptsT(:,2);
    
    clear ptsT pts
end
