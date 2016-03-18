function anaFrames(dPath,vPath,tPath,cPath,p,roi,...
                   includeCalibration,startFrame,newMean,endFrame,skipFrame)
% Steps thru frames, creates a thumbnail image and then extracts the
% midline and eye positions


%% Parameters

% Report progress at command line
reportStatus = 1;


%% Initialize data structure 'B'

% Load filenames of frames
% a = dir([vPath  filesep '*.' p.nameSuffix]);

% Alternate to load filenames (use when 'ghost' files present) 
a = dir([vPath filesep 'exp*']);

if isempty(a)
    error('No video frames found');
end

% Look for data file
a2 = dir([dPath filesep 'blob data.mat']);

% If data exist . . .
if ~isempty(a2)
    
    % Load data 'B'
    load([dPath filesep 'blob data.mat']);
    
% Otherwise . . .
else
    
    % find total number of video frames
    frTot = length(a);
    
    % initialize data structure 'B' with fields 'fr_num' and 'filename'
    B(frTot).fr_num = frTot;
    B(frTot).filename = [];
    
    % Loop thru video frames
    for i = 1:frTot
        
        % Read frame number
        frNum = str2num(a(i).name(end-p.num_digit_frame-length(p.nameSuffix):...
            end-length(p.nameSuffix)-1));
        
        % Store
        B(i).fr_num     = frNum;
        B(i).filename   = a(i).name;
        
    end
    
    % Save data 'B' to file 'blob data.mat'
    save([dPath filesep 'blob data.mat'],'B')
end

% Set default end frame and skipFrame
if nargin < 11
    skipFrame = 0;
    if nargin < 10
        endFrame = length(B);
    end
end

if includeCalibration
    % Load calibration data ('cal')
    tmp = load([cPath filesep 'calibration data.mat']);
    cal = tmp.cal;
    clear tmp
else
    cal = [];
end

% Make frame list for analysis
frames = startFrame:(1+skipFrame):endFrame;


%% Make Mean image

imMean = makeMeanImage(dPath,vPath,cal,B,newMean);


%% Loop thru frames sequentially, finding the body midline

% Update status
disp('      Analyzing frames . . .')

% Start timer
tstart = tic;

% Check for consecutive numbers
if max(diff(frames)~=1)
    error('Frame numbers not consecutive');
end

% Use a placeholder for initial iteration
blob = nan;
pBlob = blob;

% Use a placeholder for initial iteration
eye = nan;

% If not all frames analyzed . . .
% if length(dir([tPath filesep '*.mat']))<length(frames) 
    % Loop thru frames
    for i =  1:length(frames)
        %warning('TODO: reset starting frame to "1"')
        
        % Current frame
        cFrame = frames(i);
        
        
        % Read frame
        [im,~] = imread([vPath filesep B(cFrame).filename]);
        
        if ~isempty(cal)
            % Apply calibration to undistort image
            im = undistortImage(im, cal.cameraParams,'OutputView','full');
        end
        
        % Find blobs that define the fish in frame
        blob = findBlobs(im,imMean,roi);
        
        % If blob is not a nan . . .
        if ~isnan(blob.im(1))
            
            % Attempt to find landmarks
            try
                % Code that finds landmarks

                blob = findMidline(p,blob,pBlob,im);
                
                % If error above . . .
            catch
                % Write data
                blob.xMid = nan;
                
            end
            
            if 0
                figure(10)
                imshow(blob.im,[],'InitialMagnification','fit')
                hold on
                plot(blob.xMid,blob.yMid,'go-')
                hold off
                pause%(0.1)
            end
            
            % Visualize frame for debugging (switch "parfor" loop to "for")
            % (Must be commented if using parfor)
%             visData(blob,im,['Frame ' num2str(cFrame)])
            
            % Store time
            mid.t(i,1) = i./p.framerate;
            
            % Store predator blob roi
            mid.roi_blob(i,:) = blob.roi_blob;
         
            % Store coordinates
            if ~isnan(blob.xMid)
                % Store away data in global FOR
                mid.xRost(i,1)  = blob.xMid(1) + blob.roi_blob(1);
                mid.yRost(i,1)  = blob.yMid(1) + blob.roi_blob(2);
                mid.xHead(i,1)  = blob.xMid(2) + blob.roi_blob(1);
                mid.yHead(i,1)  = blob.yMid(2) + blob.roi_blob(2);
            else
                mid.xRost(i,1)  = nan;
                mid.yRost(i,1)  = nan;
                mid.xHead(i,1)  = nan;
                mid.yHead(i,1)  = nan;
            end
            
            % Store prey coordinates
            if blob.preyOn
                prey.t(i,1)         = mid.t(i,1);
                prey.xPrey(i,1)     = blob.xPrey;
                prey.yPrey(i,1)     = blob.yPrey;
                prey.thetaPrey(i,1) = blob.thetaPrey;
                prey.size(i,1)      = blob.sizePrey;
            else
                prey.xPrey(i,1)     = nan;
                prey.yPrey(i,1)     = nan;
                prey.thetaPrey(i,1) = nan;
                prey.size(i,1)      = nan;
            end
                
        end
  % ------ NOTE: pBlob is not used in findLandmarks        
        % Store current blob for next iteration
        pBlob = blob;
        
        % Run analysis code
        %blob = ana_frame(action,cFrame,frames,tPath,vPath,dPath,B,p,...
        %    cal,roi,imMean,blob);
        
        % Update status
        disp(['          Midline done for ' num2str(i) ' of ' ...
            num2str(length(frames))])
        
    end
% end


%% Close execution
    
% Report duration
telap = toc(tstart)/60;
disp(['     . . . completed in ' num2str(telap) ' min'])

% Save data 'B'
save([dPath filesep 'blob data.mat'],'B')

% Save data 'mid'
save([dPath filesep 'midline data.mat'],'mid')

% Save data 'prey'
save([dPath filesep 'prey data.mat'],'prey')


function B = setB(B,Bc,Btmp,cFrame)     
% Store results
B(cFrame).sMid       = Bc.sMid;
B(cFrame).xMid       = Bc.xMid;
B(cFrame).yMid       = Bc.yMid;
B(cFrame).xEye       = Bc.xEye;
B(cFrame).yEye       = Bc.yEye;
B(cFrame).xPerim     = Btmp.xPerim;
B(cFrame).yPerim     = Btmp.yPerim;
B(cFrame).roi_blob   = Btmp.roi_blob;

function [im,cmap] = readFrame(vPath,B,cFrame)
% reads current frame, 'im' is an MxN array containing image data
[im,cmap] = imread([vPath filesep B(cFrame).filename]);

function lg = setLog(val)
lg.success = val;

function saveBlob(tPath,blob,fileName)

% Save blob data
save([tPath filesep fileName(1:end-4)],'blob')

function saveEyes(tPath,eye,fileName)

% Save eye data
save([tPath filesep fileName(1:end-4)],'eye','-append')

function writeFile(pathh,lg)
% Allows writing a file within 'parfor' loop
save(pathh,'lg')

function visData(blob,imWholecl,ttxt)

%TODO:fix

% Read full frame
%[imWholecl,cmap] = imread([vPath filesep ...
%p.filename{cFrame}(1:(end-3)) ...
%p.nameSuffix]);
%warning off
subplot(1,2,1)
imshow(imWholecl,[],'InitialMagnification','fit')
hold on

subplot(1,2,2)
imshow(blob.im,[],'InitialMagnification','fit')
hold on
h = plot(blob.xMid,blob.yMid,'g-o',...
         blob.xEye,blob.yEye,'yo');
set(h(2),'MarkerFaceColor','y','MarkerSize',3)
title(ttxt)
hold off
%warning on
pause(0.01)

function blob = findBlobs(imStart,imMean,roi)

% Adjust grayscale values  
 
%---SKIP imadjust steps if photos have been preprocessed in Photoshop----
im     = (imadjust(imStart));
imSub  = (imadjust(imMean));

% im = imStart;
% imSub = imMean;

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

% Get inverse of image
im = imcomplement(im);

% Use roi to crop image
%---------------NOTE: if using full frame, comment out roiI
% roiI = roipoly(im,roi.x,roi.y);

% Find threshold
tVal = min([0.95 graythresh(im)+0.1]);

% Threshold image
imBW    = ~im2bw(im,tVal); % & roiI;

% % Dilate im & get properties
% se    = strel('disk',6);
% imBW = imdilate(imBW,se);
% %imBW = imerode(imBW,se);
% LL    = bwlabel(imBW);
% props = regionprops(LL,'Centroid','Area');

% Image showing intersection of roi boundary & fish
% imOverlap = ~roiI & imBW;

% Get peripheral shapes
% By = bwboundaries(imBW,'noholes');

% Get connected components 
% ccBlobs = bwconncomp(imBW);

% Get peripheral shapes, query size and boundaries
By2 = regionprops(imBW,'ConvexHull','ConvexArea','ConvexImage',...
    'BoundingBox','Orientation','Centroid');

% If no fish or touching a wall . . .
if isempty(By2) %|| (max(imOverlap(:))==1)
%     if isempty(By)
%         % Declare warning
%         warning(['No fish found in ' p.filename{i}]);
%     else
%         % Declare warning
%         warning(['Fish touching wall in ' p.filename{i}]);
%     end
    
    % nan blob
    imBlob = nan;
    
    % Store nans
    blob.xPerim   = nan;
    blob.yPerim   = nan;
    blob.roi_blob = nan;
    
% If fish . . .
else
    
    %----NOTE----
    % This may be a good place to use 'regionprops' instead of
    % 'bwboundaries' to get the largest AND second largest blob 
    
    % Filter out small objects (<15 px) found by 'regionprops'
    By2 = By2([By2.ConvexArea] > 15);
    
    % find largest object and get is bounding box
    [~,ind2] = max([By2.ConvexArea]);
    perim2 = By2(ind2).BoundingBox;
    
    % find smaller object
    [~,ind3] = min([By2.ConvexArea]);
    
    % check to see if there is a prey in this frame
    if length(By2) > 1 
        % turn on indicator for prey
        blob.preyOn = 1;
        
        % save relevent prey data
        blob.xPrey      = By2(ind3).Centroid(1);
        blob.yPrey      = By2(ind3).Centroid(2);
        blob.thetaPrey  = By2(ind3).Orientation;
        blob.sizePrey   = By2(ind3).ConvexArea;
    else
        % turn off indicator for prey
        blob.preyOn = 0;
    end
    
    % Number of pixels that pad the blob
    pad_val = 10;
    
    % Rectangle for blob roi: [XMIN YMIN WIDTH HEIGHT]
    rect = [perim2(1)-pad_val, perim2(2)-pad_val,...
                         perim2(3)+ 2*pad_val, perim2(4)+2*pad_val];
                     
    % Crop down image to predator blob using 'BoundingBox' output
    % crop image syntax: imcrop(im,[XMIN YMIN WIDTH HEIGHT])
    imBlob = imcrop(im,rect);
    
    imBlobBW = bwconvhull(imBlob<(tVal*255),'objects');
    
    % Store predator data
%     blob.xPerim    = perim(:,2);
%     blob.yPerim    = perim(:,1);
    blob.roi_blob  = rect;
    
end

% If there is a blob, find binary blobs (big and small)
if ~isnan(blob.roi_blob(1))
    
    
    % Find largest blob
%     [imBW,numBlob] = returnBlob(~im2bw(imBlob,tVal));
    [imBW,numBlob] = returnBlob(imBlobBW);
    
    % Create tmp image (uses blob as mask)
    tmp = imBlob.*0+2.^16;
    tmp(imBW(:)) = imBlob(imBW(:));
    
    % Starting low threshold value (for small blob)
    tVal2 = tVal.*.75;
    
    % Increment to increase threshold (for small blob)
    tIncr = .05;
    
    % Loop to create image at lowest threshold that finds single blob
    while true
        
        % Find blob
        [imBWsmall,numBlob] = returnBlob(~im2bw(tmp,tVal2));
        
        % Check number of blobs
        if numBlob==1
            break
        else
            tVal2 = tVal2 + tIncr;
        end
        
        % Check for too high a threshold
        if tVal2 >= tVal
            tVal2 = tVal;
            imBWsmall = imBW;
            break
        end
    end

% Otherwise, return nans
else
    imBW      = nan;
    imBWsmall = nan;
end

% Store images
blob.im      = imBlob;
blob.BW      = imBW;
blob.BWsmall = imBWsmall;


function [imBW,numBlob] = returnBlob(imBW)
% Returns binary image of blobs

% Close gaps with dilation and erosion (should be able to use 'imclose')
se   = strel('disk',5);
imBW = imdilate(imBW,se);
imBW = imerode(imBW,se);

% Identify blobs (BWCONNCOMP is more memory efficient) 
CC = bwconncomp(imBW);
LL = labelmatrix(CC);
% LL    = bwlabel(imBW);

% Get peripheral shapes
[bb,~] = bwboundaries(imBW,'noholes');

% Select blob with greatest periphery
maxB = 0;
idx = [];
for j = 1:length(bb)
    if length(bb{j}) > maxB
        maxB = length(bb{j});
        idx = j;
    end
end

if isempty(idx)
    imBW = nan;
    numBlob = 0;
else
    % Define image as having only largest blob
    imBW = LL==idx;
    
    % Return number of blobs
    numBlob = max(LL(:));
end




