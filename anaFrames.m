function motion = anaFrames(action,dPath,vPath,tPath,cPath,p,roi,...
                            skipFrame,startFrame,endFrame)
% Steps thru frames, creates a thumbnail image and then extracts the
% midline and eye positions


%% Parameters

% Report progress at command line
reportStatus = 1;


%% Initialize data structure 'B'

% Load filenames of frames
a = dir([vPath  filesep '*.' p.nameSuffix]);

if isempty(a)
    error('No video frames found');
end

% Look for data file
a2 = dir([dPath filesep 'blob data.mat']);

% If data exist . . .
if ~isempty(dir([dPath filesep 'blob data.mat']))
    
    % Load data 'B'
    load([dPath filesep 'blob data.mat']);

% Otherwise . . .
else
  
    % Loop thru video frames
    for i = 1:length(a)
        
        % Read frame number
        frNum = str2num(a(i).name(end-p.num_digit_frame-length(p.nameSuffix):...
                        end-length(p.nameSuffix)-1));

        % Store
        B(i).fr_num     = frNum;
        B(i).filename   = a(i).name;
            
    end
    
    % Save data 'B'
    save([dPath filesep 'blob data.mat'],'B')
end

% Set default start and end frames
if nargin < 10
    endFrame = length(B);
    if nargin < 9
        startFrame = 1;
    end
end

% Load calibration data ('cal')
tmp = load([cPath filesep 'calibration data.mat']);
cal = tmp.cal;
clear tmp

% Make list for analysis
frames = startFrame:(1+skipFrame):endFrame;


%% Analyze mean image for motion

% Get mean image
[imMin,imMean] = makeMeanImage(dPath,vPath,cal,B);

blobMin = findBlobs(imMin,imMean,roi);
propsMin = regionprops(blobMin.BW,'Area');

if propsMin.Area < p.minBlobArea    
    motion = 0;
    disp('Fish did not move during the recording')
    return
else 
    motion = 1;
end

clear blobMin propsMin imMin
    

%% Initialize execution

% Update status
disp('      Analyzing frames . . .')

% Start timer
tstart = tic;


%% Loop thru frames in parallel

% Step though frames (Parallel for preliminary run)
if strcmp(action,'parallel')
    
    % Set up parallel pool (2 cores)
    %parpool(2)
    
    % Set up log directory
    if isempty(dir([dPath filesep 'log']))
        mkdir([dPath filesep 'log'])
    else
        delete([dPath filesep 'log' filesep '*']);
    end

    % Loop thru frames in parallel
    %parfor i =  1:length(frames); 
    for i =  1:length(frames); 
        
        % Current frame
        cFrame = frames(i);

        % Run analysis code
        blob = ana_frame(action,cFrame,frames,tPath,vPath,dPath,B,p,...
                         cal,roi,imMean);
        
        % If it failed . . .            
        if isnan(blob.xMid(1))
             % Log result
            lg = setLog(0);
        else
             % Log result
            lg = setLog(1);
        end
        
        % Write to log (to track progress)
        writeFile([dPath filesep 'log' filesep B(cFrame).filename(1:end-4)],lg);
        
        % Look at log
        aaa = dir([dPath filesep 'log']);
        
        % Update status
        if reportStatus
            if ~isnan(blob.xMid(1))
                disp(['           Frame succeeded :)  ' ...
                    num2str(floor(100*length(aaa)/length(frames))) ...
                    '% of blob analysis completed'])
            else
                disp(['           Frame failed :(  ' ...
                    num2str(floor(100*length(aaa)/length(frames))) ...
                    '% of blob analysis completed'])
            end
        end

    end
    
    % Remove log files and directory
    if ~isempty(dir([dPath filesep 'log']))
        pause(0.1)
        delete([dPath filesep 'log' filesep '*']);
        rmdir([dPath filesep 'log']);
    end
end  
 

%% Loop thru frames sequentially

% Run, not parallel, for fast start
if strcmp(action,'sequential')
    
    % Check for consecutive numbers
    if max(diff(frames)~=1)
        error('Frame numbers not consecutive');
    end
    
    % Use a placeholder for initial iteration
    blob = nan;
    
    % Loop thru frames
    for i =  1:length(frames); 
        
        % Current frame
        cFrame = frames(i);
        
        % Run analysis code
        blob = ana_frame(action,cFrame,frames,tPath,vPath,dPath,B,p,...
                         cal,roi,imMean,blob);
        
        % Look at log
        aaa = dir([dPath filesep 'log']);
        
        % Update status        
        if reportStatus
            if ~isnan(blob.xMid(1))
                disp(['           Frame succeeded :)  ' ...
                    num2str(round(100*i/length(frames))) ...
                    '% of blob analysis completed'])
            else
                disp(['           Frame failed :(  ' ...
                    num2str(round(100*i/length(frames))) ...
                    '% of blob analysis completed'])
            end
        end   
    end
end


%% Close execution
    
% Report duration
telap = toc(tstart)/60;
disp(['     . . . completed in ' num2str(telap) ' min'])

% Save data 'B'
save([dPath filesep 'blob data.mat'],'B')




function blob = ana_frame(action,cFrame,frames,tPath,vPath,...
                          dPath,B,p,cal,roi,imMean,pBlob)

% Set default
if (nargin < 9) 
    bBlob = nan;
end

% Read frame
[im,cmap] = readFrame(vPath,B,cFrame);
[im,cmap] = imread([vPath filesep B(cFrame).filename]);

% Apply calibration to undistort image
im = undistortImage(im, cal.cameraParams,'OutputView','full');

% Find blobs that define the fish in frame
blob = findBlobs(im,imMean,roi);

% If blob is not a nan . . .
if ~isnan(blob.im(1))
    
    % Attempt to find landmarks
    try
        % If looking at fast start (beyond first frame) . . .
        if strcmp(action,'sequential') && (cFrame>frames(1))
            % Code that finds landmarks
            blob = findLandmarks(blob,pBlob);
            
        % Otherwise . . .
        else
            % Code that finds landmarks
            blob = findLandmarks(blob);
        end

        
    % If error above . . .
    catch   
        
        % Write data
        blob.xMid = nan;
        
    end
    
    % Visualize frame for debugging (switch "parfor" loop to "for")
    % (Must be commented if using parfor)
    %visData(blob,im,['Frame ' num2str(cFrame)])

    % Save thumbnail file
    saveBlob(tPath,blob,B(cFrame).filename)
    
end


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
[im,cmap] = imread([vPath filesep B(cFrame).filename]);


function lg = setLog(val)
lg.success = val;


function saveBlob(tPath,blob,fileName)

% Save blob data
save([tPath filesep fileName(1:end-4)],'blob')

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

% Adjust grayscale values and convert to double
im     = (imadjust(imStart));
imSub  = (imadjust(imMean));

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

% Get inverse of image
im = imcomplement(im);

% Use roi to crop image
roiI = roipoly(im,roi.x,roi.y);

% Find threshold
tVal = min([0.95 graythresh(im)+0.1]);

% Threshold image
imBW    = ~im2bw(im,tVal) & roiI;

% % Threshold value
% %tVal = min([0.95 graythresh(im)+0.1]);
% %imBW    = ~im2bw(im,tVal) & roiI;
% 
% % Dilate im & get properties
% se    = strel('disk',6);
% imBW = imdilate(imBW,se);
% %imBW = imerode(imBW,se);
% LL    = bwlabel(imBW);
% props = regionprops(LL,'Centroid','Area');

% Get peripheral shapes
By = bwboundaries(imBW,'noholes');

% Image showing intersection of roi boundary & fish
imOverlap = ~roiI & imBW;

% Store file info
%pd(i).frame = p.frNums(i);
%pd(i).filename = p.filename{i};

% If no fish or touching a wall . . .
if isempty(By) || (max(imOverlap(:))==1)
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
    
    % Select blob with greatest periphery
    maxB = 0;
    for j = 1:length(By)
        if length(By{j}) > maxB
            maxB = length(By{j});
            perim = By{j};
        end
    end
    
    % Number of pixels that pad the blob
    pad_val = 10;
    
    % Min and max range of roi
    Croi_min = max([1 min(perim(:,2))-pad_val]);
    Rroi_min = max([1 min(perim(:,1))-pad_val]);
    Croi_max = min([size(im,2) max(perim(:,2)+pad_val)]);
    Rroi_max = min([size(im,1) max(perim(:,1)+pad_val)]);
    
    % Crop down to blob
    imBlob = imcrop(im,[Croi_min Rroi_min ...
                        Croi_max-Croi_min Rroi_max-Rroi_min]);
    
    % Enhance blob contrast
    %imBlob = adapthisteq(imBlob,'clipLimit',0.02,'Distribution','rayleigh');
    %imBlob = imadjust(imBlob);
    
    % Store data
    blob.xPerim    = perim(:,2);
    blob.yPerim    = perim(:,1);
    blob.roi_blob  = [Croi_min Rroi_min Croi_max-Croi_min Rroi_max-Rroi_min];
    
    %clear perim Croi_min Rroi_min Croi_max Rroi_max
end

% If there is a blob, find binary blobs (big and small)
if ~isnan(blob.xPerim(1))
    
    
    % Find largest blob
    [imBW,numBlob] = returnBlob(~im2bw(imBlob,tVal));
    
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

if isempty(idx)
    imBW = nan;
    numBlob = 0;
else
    % Define image as having only largest blob
    imBW = LL==idx;
    
    % Return number of blobs
    numBlob = max(LL(:));
end




