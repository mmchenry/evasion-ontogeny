function pred_prey_acq(vPath,skip_prey)

% Acquire kinematics of predator and prey fish



%% Parameter values

% Extension for image files
nameSuffix = 'tif';

% Number of digits for frame number in filename
num_digit = 6;

% Initial visualizing predator acquisition
visSteps = 1;

% Number of predtaor frames to visualize at start of acquisition
numVis = 50;

% Max number of frames for creating the mean image
maxFrames = 1000;

% Whether to invert the image
invert = 1;

% Radius of prey roi (in pixels)
py_roi = 30;

if nargin < 2
    skip_prey = 0;
end

% Frame rate (fps)
p.framerate = 1000;



%% Get path of data file, load data

% if nargin < 1
%     vPath = uigetdir(pwd,'Select directory');
%     if vPath==0
%         return
%     end
% end

% Path to video
vPath = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Video for Matt';

% Path to data
dPath = vPath;

% Load filenames for frames
a = dir([vPath  filesep '*.' nameSuffix]);

if isempty(a)
    error('No video frames found');
end

% Get indicies for video frames
for i = 1:length(a)
    
    % Read filename
    frNum = str2num(a(i).name(end-num_digit-length(nameSuffix):...
                    end-length(nameSuffix)-1));
                
    % Check ordering    
    if (i>1) && (frNum ~= p.frNums(i-1)+1)
        error('Frame numbers not consecutive')
    end
    
    % Store
    p.frNums(i,1) = frNum;
    p.filename{i} = a(i).name;
end


%% Define roi

% Look for roi data
a2 = dir([vPath filesep 'roi.mat']);

if isempty(a2)  
    % Read first frame
    im = imread([vPath  filesep a(1).name]);
    
    % Select dimensions of roi
    txt = 'Clockwise, starting with upper-left corner';
    figure;
    [p.roi_x,p.roi_y]   = choosePoints(im,1,txt);
    close 
    
    % Trim extraneous points
    p.roi_x = p.roi_x(1:4);
    p.roi_y = p.roi_y(1:4);      
    
    % Save roi data
    save([vPath filesep 'roi.mat'],'p')  
    
else
    load([vPath filesep 'roi.mat'])
    
end

% Clear variables
clear im txt a2


%% Prompt for sequence info

% % Look for mean image
% a2 = dir([vPath filesep 'start_point.mat']);
% 
% if isempty(a2)
%     
%     % Load p for defined for roi
%     load([vPath filesep 'roi.mat'])
%     
%     warning off all
%     
%     % Get indicies for video frames
%     idx = 1;
%     for i = 1:length(a)
%         frNum = str2num(a(i).name(end-num_digit-length(nameSuffix):...
%             end-length(nameSuffix)-1));     
%         if (frNum >= startFrame) && (frNum <= endFrame) 
%             p.frNums(idx) = frNum;
%             p.filename{idx} = a(i).name;
%             idx = idx + 1;
%         end 
%     end
%     
%     % Measure body length, initial position & orientation
%     txt = 'Select nose, then caudal peduncle';
%     img = imread([vPath filesep p.filename{1}]);
%     [xT,yT]   = choosePoints(img,1,txt);
%     p.bLength = ((xT(2)-xT(1))^2 + (yT(2)-yT(1))^2)^0.5;
%     p.xHead = xT(1);
%     p.yHead = yT(1);
%     p.xTail = xT(2);
%     p.yTail = yT(2);
%     p.x = mean(xT);
%     p.y = mean(yT);
%     clear xT yT txt
%     close
% 
%     warning on all
% 
%     save([vPath filesep 'start_point.mat'],'p');
%     
%     clear startFrame endFrame prompt name numlines defaultanswer a
%     
% else % if seq_param exists, load
% 
%     disp(' '); disp('Loading existing starting point data . . .'); 
%     load([vPath filesep 'start_point.mat'])
% 
% end
% 
% clear img


%% Create or load mean image

% Look for mean image
a2 = dir([vPath filesep 'meanImage.tif']);

% Calculate mean image does not exist
if isempty(a2)   
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if length(p.frNums) > maxFrames
        dframe = floor(length(p.frNums)/maxFrames);
        frIdx = 1:dframe:length(p.frNums);
        clear dframe
    else
        frIdx = 1:length(p.frNums);
    end
    
    % Create waitbar
    h = waitbar(0,...
            ['Mean image: ' num2str(1)],...
             'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vPath  filesep p.filename{1}]);
    imSum = double(imCurr);
    clear imCurr tmp
      
    % Loop through frames 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vPath  filesep p.filename{frIdx(i)}]);
        imSum        = imSum + double(imCurr);
        clear tmp imCurr
        
        % Update status bar
        h = waitbar(i/length(frIdx),h,...
            ['Mean image: ' num2str(i) ' of ' num2str(length(frIdx)) ' frames']);
        
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
        
    end
    
    % Calculate mean from sum image
    imMean = uint8(round(imSum./length(frIdx)));
    
    imMean = imMean(:,:,1);
    
    % Write image to movie dir
    imwrite(imMean,[vPath filesep 'meanImage.tif'],'tif',...
            'Compression','none');
    
    close force
    clear frIdx h i imSum
        
    %imMean = rgb2gray(imMean);
      
    
% Load mean image, if present
else
    
    disp(' ')
    disp('Loading mean image . . .');
    imMean = imread([vPath filesep 'meanImage.tif']);
    
end


%% Select thresholds

% Look for mean image
a2 = dir([vPath filesep 'seq_params.mat']);

if 0 % isempty(a2)
    
    % Grab frames for threshold finding
    img = grabFrame(vPath,p.filename{1},invert,p.roi_x,p.roi_y);
    
    % Matlab guesses a threshold value
    p.tVal = graythresh(img)+0.1;
    
    % Store path info in p
    p.path   = [vPath ];

    % Run threshFinder to find threshold values for predator and prey
    % note: threshFinder saves p in seq_params.mat
    disp(' ')
    disp('Choose threshold for the predator')
    
    %waitfor(threshFinder(img,p))
    
    load([vPath filesep 'seq_params.mat'])
 
    
%     % Locate new postion of larva
%     [p.xPrey,p.yPrey,p.areaPrey,x_roipy,y_roipy,imBW] = findLarva(img,p.xPrey, ...
%                                 p.yPrey,50,py_roi,p.tVal_py,'dist');
    
    
%     % Define prey roi 
%     x_py = p.xPrey + py_roi.*cos(linspace(0,2*pi,200));
%     y_py = p.yPrey + py_roi.*sin(linspace(0,2*pi,200));
%     pyROI = roipoly(img,x_py,y_py);
%     imBW  = ~im2bw(img,p.tVal_py) & pyROI;
%     
%     % Dilate im & get properties
%     se    = strel('disk',4,4);
%     imBW = imdilate(imBW,se);
%     LL    = bwlabel(imBW);
%     props = regionprops(LL,'Centroid','Area');
%     
%     % Select blob closest to selected point
%     if length(props)>1
%         dist = 10^10;
%        for i = 1:length(props)
%            cDist = sqrt((props(i).Centroid(1)-p.xPrey)^2 + ...
%                         (props(i).Centroid(2)-p.yPrey)^2);
%            if cDist < dist
%                dist = cDist;
%                py.areaPrey = props(i).Area;
%            end
%        end
%     else
%         p.areaPrey = props.Area;
%     end
    
    % Save 'p' structure
    save([vPath filesep 'seq_params'],'p')
    
    clear im img x_roi y_roi dist props LL imBW pyROI x_py y_py se i
    
else % if seq_param exists, load

%     disp(' '); disp('Loading existing sequence parameters . . .'); 
%     load([vPath filesep 'seq_params.mat'])
%     disp(' ');

end

clear img


%% Step through frames for position of predator

a3 = dir([vPath filesep 'pred_coords.mat']);

if isempty(a3)
    
    f = figure;
    set(f,'DoubleBuffer','on')
    
    % Get image
    img = grabFrame(vPath,p.filename{1},invert,p.roi_x,p.roi_y);
    
    % Gues for initial threshold value
    p.tVal = graythresh(img)+0.1;
    
    % Loop through frames
    for i = 1:length(p.frNums)
        
        % Define roi coordinates
        %[x_roi,y_roi] = roiCoords(p);
        
        % Read full frame
        im = imread([vPath  filesep p.filename{i}]);
        
        % Grab frame, threshold image & choose roi
        img = grabFrame(vPath,p.filename{i},invert,p.roi_x,p.roi_y);
        imROI   = roipoly(img,p.roi_x,p.roi_y);
        %& ...
         %         roipoly(img,pd.xPerim{i},pd.yPerim{i});
        img2 = img.*0;
        img2(imROI) = img(imROI);
        imBW    = ~im2bw(img,p.tVal);
        imBW    = imBW & imROI;
        
        % Dilate im & get properties
        se    = strel('disk',4,4);
        imBW = imdilate(imBW,se);
        %imBW = imerode(imBW,se);
        LL    = bwlabel(imBW);
        props = regionprops(LL,'Centroid','Area');
        
        % Get peripheral shapes
        [B,L] = bwboundaries(imBW,'noholes');
        
        % Select blob with greatest periphery
        maxB = 0;
        for j = 1:length(B)
            if length(B{j}) > maxB
                maxB = length(B{j});
                perim = B{j};
            end
        end
        
        imBlob{i} = imcrop(img,[min(perim(:,2)) min(perim(:,1)) ...
                            range(perim(:,2)) range(perim(:,1))]);
        
        %imBlob2 = adapthisteq(imBlob,'clipLimit',0.02,'Distribution','rayleigh');                
        %imBlob = imadjust(imBlob);     
        
        % Store away data
        pd(i).frame = p.frNums(i);
        pd(i).filename = p.filename{i};
        pd(i).xPerim = perim(:,1);
        pd(i).yPerim = perim(:,2);
        
        
        
        % Visualize frames
        if visSteps
  
            figure(f)
            warning off
            subplot(3,1,1:2)
            imshow(im)
            hold on
            plot(pd(i).yPerim,pd(i).xPerim,'r-')
            title(['Frame ' num2str(p.frNums(i)) ' of ' ...
                num2str(p.frNums(end))])
            hold off
            subplot(3,1,3)
            imshow(imBlob)
            pause(.2)
            warning on
            clear im
        else
            disp(['Predator acquire: Frame ' num2str(p.frNums(i)) ' of ' ...
                num2str(p.frNums(end))]);
        end
        
        % Clear variables for next loop
        clear im img imBW imBW2 props imROI se x_roi y_roi maxB
%         
%         if i > numVis
%             close
%             visSteps = 0;
%         end
    end
    
    % Save data
    save([vPath filesep 'pred_coords'],'pd')
    
    % Save blob images
    save([vPath filesep 'Blobs'],'imBlob')
    
else
    
    % Load 'pd' structure of predator coordinates
    disp('Loading predator data . . .')
    load([vPath filesep 'pred_coords.mat'])
    
    % Load blob images ('imBlob')
    load([vPath filesep 'Blobs'])

end


%% Analyze blobs

% Define start and end frames
first_frame = 1;
last_frame = length(pd);

% Step though frames
for i = first_frame:last_frame
    
    % Extract current frame, with enhanced contrast
    imBlob2 = imadjust(imBlob{i}); 
    
    % Guess for initial threshold value
    tVal = graythresh(imBlob2)+0.1;
    
    % Apply threshold
    imBW    = ~im2bw(imBlob2,tVal);
    
    % Close gaps with dilation and erosion
    se    = strel('disk',4,4);
    imBW = imdilate(imBW,se);
    imBW = imerode(imBW,se);
    
    % Identify blobs
    LL    = bwlabel(imBW);
    props = regionprops(LL,'Centroid','Area');
        
    % Get peripheral shapes
    [B,L] = bwboundaries(imBW,'noholes');
        
    % Select blob with greatest periphery
    maxB = 0;
    idx = [];
    for j = 1:length(B)
        if length(B{j}) > maxB
            maxB = length(B{j});
            perim = B{j};
            xPerim = perim(:,1);
            yPerim = perim(:,2);
            idx = j;
        end
    end
    
    % Define image as having only largest blob
    imBW = LL==idx;
    
    % Skeletonize
    %imSkel = bwmorph(imBW,'skel',Inf);
    %[ySkel,xSkel] = find(imSkel);
        
    % Skeletonize blob
    skel= bwmorph(imBW,'skel',Inf);
    
    % Find branch and endpoints
    B = bwmorph(skel, 'branchpoints');
    E = bwmorph(skel, 'endpoints');
    
    % Get coordinates for blob, branch & end points
    [yEnd,xEnd]    = find(E);
    [yBran,xBran]  = find(B);
    [ySkel,xSkel]  = find(skel);
    
    % FIND FIRST HEAD & TAIL POINTS ---------------------------------------
    
    % Identify tail on first run thru code (works best for straight body)
    if i==first_frame
        
        % Compute distances between points
        for j = 1:length(xEnd)
            dists(:,j) = sqrt((xEnd-xEnd(j)).^2 + (yEnd-yEnd(j)).^2);
        end
        
        % Get indices of points furthest apart
        [iR,iC] = find(dists==max(dists(:)));
        
        % Theta & radius for ROI
        theta = linspace(0,2*pi,100);
        r = max(size(B))/8;
        
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
    end
          
    % FIND NEW HEAD & TAIL POINTS -----------------------------------------
    
    % Distance between prior and tail point and endpoints
    xTail = repmat(xTail,length(xEnd),1);
    yTail = repmat(yTail,length(yEnd),1);
    dists = sqrt((xEnd-xTail).^2 + (yEnd-yTail).^2);
    
    % Index for new tail point as min distance to endpoint
    iVal = find(dists==min(dists),1,'first');
    xTail = xEnd(iVal);
    yTail = yEnd(iVal);
    
    % Distance between prior and tail point and endpoints
    xHead = repmat(xHead,length(xEnd),1);
    yHead = repmat(yHead,length(yEnd),1);
    dists = sqrt((xEnd-xHead).^2 + (yEnd-yHead).^2);
    
    % Index for new tail point as min distance to endpoint
    iVal = find(dists==min(dists),1,'first');
    xHead = xEnd(iVal);
    yHead = yEnd(iVal);

    % FIND MIDLINE POINTS -------------------------------------------------
 
    % Set last point (on repeat)
    %xLast = repmat(xTail,length(xSkel),1);
    %yLast = repmat(yTail,length(xSkel),1);
    %xLast = xTail;
    %yLast = yTail;
    
    % Prune branches from skeleton
    skel = prune_branches(skel,xHead,yHead,xTail,yTail);
        
    % Remove branches from binary image
    %skel = remove_branches(skel,B,E,xHead,yHead,xTail,yTail);
    
    [sMid,xMid,yMid] = find_mid(skel,xHead,yHead,xTail,yTail);
    
    
    % Find midline points
    %[yMid,xMid] = find(skel);
   
     % Visualize frames
     if visSteps
        
        warning off
       % subplot(1,2,1)
        imshow(imBlob2)
        hold on
        plot(yPerim,xPerim,'r-')
        plot(xMid,yMid,'g+')
        title(['Frame ' num2str(i)])
        hold off
        warning on
        pause(0.01)
     end
    
     clear imBlob2 imBW imBW2 LL props B L maxB idx dists xEnd yEnd
end

return


%TODO: postprocessing on predator to determine head


% function  [skel,xLast,yLast,isBranch,isEnd] = chooseNext(skel,...
%                  xLast,yLast,xBran,yBran,xEnd,yEnd)
% % Selects next-closest point among skeleton points and then removes it
% 
% [ySkel,xSkel]  = find(skel);
% 
% % Repeat values
% xLast = repmat(xLast,length(xSkel),1);
% yLast = repmat(yLast,length(xSkel),1);
% 
% % Distance between all skeletal points and last tail point
% dists   = sqrt((xSkel-xLast).^2 + (ySkel-yLast).^2);
%     
% % Index of closest next point
% iClose = find(dists==min(dists),1,'first');        
% 
% % Determine if branch point
% if max( (xSkel(iClose)==xBran) & (ySkel(iClose)==yBran) )
%     isBranch = 1;
% else
%     isBranch = 0;
% end
% 
% % Determine if end point
% if max( (xSkel(iClose)==xEnd) & (ySkel(iClose)==yEnd) )
%     isEnd = 1;
% else
%     isEnd = 0;
% end
% 
% % Set next-closest as new last point (on repeat)
% xLast = xSkel(iClose);
% yLast = ySkel(iClose);
% 
% % Remove the next point from the skeleton
% skel(yLast,xLast) = false;

function  [s,xSm,ySm] = find_mid(skel,xHead,yHead,xTail,yTail)
% Finds a smooth midline from skeletonized image

% Number of pixels between midline points
mid_space = 5;

tol = 10;

% Store iniital image
skel_start = skel;

% Pad image
skel = pad_im(skel);



% Pad coordinates
xHead = xHead + 1;
yHead = yHead + 1;
xTail = xTail + 1;
yTail = yTail + 1;

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

% Spline smoothed coordinates
[xSp,xSm] = spaps(s,xMid,tol);
[ySp,ySm] = spaps(s,yMid,tol);

% Body positions that define margins of cranium and tail
pos_cran = 1/3;
pos_tail = 2/3;

% Indices for head and tail points
iHead = (s./max(s))<pos_cran;
iTail = (s./max(s))>pos_tail;

% % Linear fit to head and tail points
% cH = polyfit(xMid(iHead),yMid(iHead),1);
% cT = polyfit(xMid(iTail),yMid(iTail),1);
% %cH = polyfit(s(iHead),yMid(iHead),1);
% %cT = polyfit(s(iTail),yMid(iTail),1);
% 
% 
% vect1 = [1 cH(1)]; % create a vector based on the line equation
% vect2 = [1 cT(1)];
% dp = dot(vect1, vect2);
% 
% % compute vector lengths
% length1 = sqrt(sum(vect1.^2));
% length2 = sqrt(sum(vect2.^2));
% 
% % obtain the larger angle of intersection in degrees
% angl = acos(dp/(length1*length2));
% 
% % Intersection between two lines
% xInt = (cH(2)-cT(2))/(cT(1)-cH(1));
% yInt = polyval(cT,xInt);
% 
% 
% xPts = 1:size(skel,2);
% 
% if 1 
%     warning off
%     imshow(skel_start)
%     hold on
%     plot(xSm,ySm,'r+',xPts,polyval(cH,xPts),'g-',xPts,polyval(cT,xPts),'y-',...
%         xInt,yInt,'go')
%     title(['Angle = ' num2str(angl/pi*180) ' deg'])
%     hold off
%     warning on
%     pause(0.1)
% end

%[xSp,xSm] = spaps(s,xMid,tol);
%[ySp,ySm] = spaps(s,yMid,tol);

%xMid = xSm
% Visual check
if 0
    warning off
    imshow(skel_start);hold on
    plot(xSm,ySm,'r-');hold off
    warning on
end


function [py_x,py_y,py_a,x_roipy,y_roipy,imBW] = ...
                findLarva(img,cX,cY,cA,py_roi,tVal,mthd)

% mthd - ('area' or 'dist') criterion for selecting blob            
            
% Define prey roi coordinates
x_roipy = cX + py_roi.*cos(linspace(0,2*pi,200));
y_roipy = cY + py_roi.*sin(linspace(0,2*pi,200));

% Binary image of the roi around the prey
pyROI = roipoly(img,x_roipy,y_roipy);

% Slice up img by the rois
imBW    = ~im2bw(img,tVal);
imBW    = imBW & pyROI;

clear img

% Dilate im & get properties
se    = strel('disk',ceil(sqrt(cA/pi)),4);
imBW = imdilate(imBW,se);
imBW = imerode(imBW,se);
LL    = bwlabel(imBW);
props = regionprops(LL,'Centroid','Area');

% Halt, if no blob
if isempty(props)
    py_x = nan;
    py_y = nan;
    py_a = nan;
    warning(['Lost larva -- try expanding the roi and/or ' ...
        'adjusting the threshold']);
    
% Store, if one blob
elseif length(props)==1
    py_x = props.Centroid(1);
    py_y = props.Centroid(2);
    py_a = props.Area;
    
    
else
    
    % Select blob with area closest to last
    if strcmp(mthd,'area')
        
        tmp = 10^10;
        for j = 1:length(props)
            if abs(cA - props(j).Area) < tmp
                tmp = abs(cA - props(j).Area);
                py_x = props(j).Centroid(1);
                py_y = props(j).Centroid(2);
                py_a = props(j).Area;
            end
        end
        
    % Select closest distance from last   
    elseif strcmp(mthd,'dist')
        
        dist = 10^10;
        for i = 1:length(props)
            cDist = sqrt((props(i).Centroid(1)-cX)^2 + ...
                (props(i).Centroid(2)-cY)^2);
            if cDist < dist
                dist = cDist;
                py_x = props(i).Centroid(1);
                py_y = props(i).Centroid(2);
                py_a = props(i).Area;
            end
        end
        
        
    else
        error('invalid entry for mthd');
        
    end
        
end



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
    



function img = grabFrame(dirPath,filename,invert,x_roi,y_roi)

% Load image
img = imread([dirPath filesep filename]);

% Deinterlace image
%img = deinterlace(img);

%img = adapthisteq(img,'clipLimit',0.02,'Distribution','rayleigh');

% Load subtraction image
imSub  = imread([dirPath filesep,'meanImage.tif']);   

% Adjust grayscale values and convert to double
im     = (imadjust(img));
imSub  = (imadjust(imSub));

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

%im(find(im>255))  = 255;

if invert
    im = imcomplement(im);
end

% Use roi to crop image
if nargin > 3
    roiI = roipoly(im,x_roi,y_roi);
    img = uint8(255.*ones(size(im,1),size(im,2)));
    img(roiI(:)) = im(roiI(:));
else
    img = uint8(255.*ones(size(im,1),size(im,2)));
end


function img = deinterlace(img)
% This version uses a single field from an interlaced video frame
% Note: code could be modified to double temporal resolution

% Get coordinates for whole frame and individual fields
[X,Y] = meshgrid(1:size(img,2), 1:size(img,1));
[X1,Y1] = meshgrid(1:size(img,2), 1:2:size(img,1));
[X2,Y2] = meshgrid(1:size(img,2), 2:2:size(img,1)-2);

% Extract fields
fld1 = img(1:2:size(img,1),:);
%fld2 = img(2:2:size(img,1),:);

% Interpolate between scan lines of field 1
warning off
fr2 = uint8(interp2(X1,Y1,double(fld1),X2,Y2));
warning on

% Replace field 2 with interpolated values
for i=1:size(fr2,1)
    img(2*i,:) = fr2(i,:);
end


function [x_roi,y_roi] = roiCoords(p)
%Provides coordinates for an elliptical region of interest

numPts  = 400;
x_h     = p.roi_h.x(1:2);
y_v     = p.roi_v.y(1:2);
r_h     = abs(x_h(1)-x_h(2))/2;
r_v     = abs(y_v(1)-y_v(2))/2;
x_roi   = [];
y_roi   = [];

theta   = linspace(0,pi/2,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi/2,pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi,1.5*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(1.5*pi,2*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];


function [x,y] = choosePoints(img,link,txt)
%Used for finding coordinate points on a static image 'img'.
warning off all
imshow(img);
title(txt)
hold on;
set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')
n = 0;
but = 1;
while 1
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
        if link
            h = plot(x,y,'ro-');
        else
            h = plot(x,y,'ro');
        end
    elseif but==3
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        hold off
        imshow(img);
        title(txt)
        hold on
        if link
            h = plot(x,y,'ro-');
        else
            h = plot(x,y,'ro');
        end
    end
end

delete(h)

x = x'; y = y';
warning on all