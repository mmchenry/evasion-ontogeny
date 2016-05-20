function findPrey(dPath,vPath,startFrame,redoPrey,includeCalibration,p)
% Steps thru frames, creates a thumbnail image and then extracts the
% position of the prey and its orientation


%% Parameters

% Report progress at command line
% reportStatus    = 1;

% Indicator for visualising progress
% figOn           = 0;

% Indicator for adjusting image contrast
adjustON        = 0;


%% Load data structure 'B' & midline data

% Load midline data
load([dPath filesep 'Midline data.mat'])

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
    error('Run anaFrames on this video first.')
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

%% Initial position of prey

if isempty(dir([dPath filesep 'prey data.mat'])) || redoPrey
    
    close all;
    
    if 0% redoPrey
        % Load prey data
        load([dPath filesep 'prey data.mat'])
        
        % Get initial position of prey
        pPrey.x = prey.x0;
        pPrey.y = prey.y0;
    else
        
        % Load image of first video frame
        im = imread([vPath filesep a(startFrame).name]);
        
        f = figure;
        imshow(im,'InitialMagnification','fit')
        title('Zoom into prey region and press return')
        zoom on;
        
        % Wait for the most recent key to become the return/enter key
        waitfor(f,'CurrentKey','return');
        zoom reset;
        zoom off;
        
        title('Select rostrum and COM')
        [x,y,~] = ginput(2);
        hold on
        plot(x(1),y(1),'+r',x,y,'r-')
        
        % Save initial position of prey
        pPrey.x = x(2);
        pPrey.y = y(2);
        
        % Store initial position of prey 
        prey.x0 = x(2);
        prey.y0 = y(2);
        
        % Set initial prey size in pixels (rough estimate)
        pPrey.size = 20;
        
        % close figure
        close(f)
    
    end
    
    clear x y

%% Make Mean image

newMean2 = 1;

imMean = makeMeanImage(dPath,vPath,cal,B,0,newMean2);


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
    
    % Loop thru frames
    for k =  1:length(frames)
        
        % Current frame
        cFrame = frames(k);
        
        % Read frame
        [im,~] = imread([vPath filesep B(cFrame).filename]);
        
        if ~isempty(cal)
            % Apply calibration to undistort image
            im = undistortImage(im, cal.cameraParams,'OutputView','full');
        end
        
        % Get predator roi
        xMin    = mid.roi_blob(k,1);
        yMin    = mid.roi_blob(k,2);
        w       = mid.roi_blob(k,3);
        h       = mid.roi_blob(k,4);
        
        roi = [xMin yMin; xMin+w yMin; xMin+w yMin+h; xMin yMin+h;...
               xMin yMin];
        
        % Find blobs that define the fish in frame
        blob = findBlobs(im,imMean,pPrey,adjustON,roi);
        
        % If blob is not a nan . . .
        if ~isnan(blob.im(1))

            % Store time
            prey.t(k,1) = k./p.framerate;
            
            % Store prey blob roi
            prey.roi_blob(k,:) = blob.roi_blob;
            
            % Store prey coordinates
            prey.xPrey(k,1)     = blob.xPrey;
            prey.yPrey(k,1)     = blob.yPrey;
            prey.thetaPrey(k,1) = blob.thetaPrey;
            
            % Store prey size parameters
            prey.size(k,1)      = blob.sizePrey;
            prey.MajorAxis(k,1) = blob.MajorAxis;
            prey.MinorAxis(k,1) = blob.MinorAxis;
            
            % Update blob pixel positions for next iteration
            pPrey.x = blob.xVals;
            pPrey.y = blob.yVals;
            
            % Set size expection for next iteration
            pPrey.size = mean(prey.size);
            
            % otherwise ... 
        else
            % Display warning
            warning(['No fish found in ' B(cFrame).filename]);
            
            % Store time
            prey.t(k,1) = k./p.framerate;

            % Store prey blob roi
            prey.roi_blob(k,:) = blob.roi_blob;
            
            % Store nans
            prey.xPrey(k,1)     = nan;
            prey.yPrey(k,1)     = nan;
            prey.thetaPrey(k,1) = nan;
            prey.size(k,1)      = nan;
            
        end
      
        % ------ NOTE: pBlob is not used in anything
        % Store current blob for next iteration
        pBlob = blob;
        
        % Update status
        disp(['          Prey done for ' num2str(k) ' of ' ...
            num2str(length(frames))])
    end    

%% Close execution

% Report duration
telap = toc(tstart)/60;
disp(['     . . . completed in ' num2str(telap) ' min'])


% Save data 'prey'
save([dPath filesep 'prey data.mat'],'prey')

else
    disp('            Loading prey data')
    load([dPath filesep 'prey data.mat'])
end

%% Plot results
figure
subplot(2,1,1)
plot(prey.xPrey,prey.yPrey)
hold on;
plot(prey.xPrey(1),prey.yPrey(1),'ok')
xlabel('x-position')
ylabel('y-position')

subplot(2,1,2)
plot(prey.t,unwrap(prey.thetaPrey))
xlabel('time (sec)')
ylabel('Orientation (deg)')


function blob = findBlobs(imStart,imMean,pPrey,adjustON,roi)

% Adjust grayscale values if adjustON indicator is set to 1

if adjustON
    im     = (imadjust(imStart));
    imSub  = (imadjust(imMean));
else
    im = imStart;
    imSub = imMean;
end

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

% Get inverse of image
im = imcomplement(im);

% Use predator roi to exclude predator region
roiP = roipoly(im,roi(:,1),roi(:,2));

% Find threshold
tVal = min([0.95 graythresh(im)+0.1]);

% Threshold image
imBW    = ~im2bw(im,tVal); %& ~roiP;

% Close image to get rid of small contamination
se    = strel('diamond',5);
imBW = imclose(imBW,se);

% Image showing intersection of roi boundary & fish
% imOverlap = ~roiI & imBW;

% Find all blobs in binary image
% props = regionprops(imBW,'Centroid','Area','Orientation',...
%     'BoundingBox','PixelList');

% Select blob that includes any of the previous positions
imBW1 = bwselect(imBW,pPrey.x,pPrey.y,8);

% Measure its properties
props2 = regionprops(imBW1,'Centroid','Area','Orientation',...
    'BoundingBox','PixelList','MajorAxisLength','MinorAxisLength');

% Filter out small objects (<0.33 px of preySize) found by 'regionprops'
props2 = props2([props2.Area] > pPrey.size/3);

% If no fish  . . .
if isempty(props2)

    % Store nans
    blob.im         = nan;
    blob.roi_blob   = nan;

% If more than one blob ... 
elseif length(props2)>1
    
    % Interactively select the prey position (press return after selection)
    imBW1 = bwselect(imBW);
    
    props2 = regionprops(imBW1,'Centroid','Area','Orientation',...
    'BoundingBox', 'PixelList','MajorAxisLength','MinorAxisLength');

% Get bounding box for prey
    perim = props2.BoundingBox;

    % Number of pixels that pad the blob
    pad_val = 5;
    
    % Rectangle for blob roi: [XMIN YMIN WIDTH HEIGHT]
    rect = [perim(1)-pad_val, perim(2)-pad_val,...
        perim(3)+ 2*pad_val, perim(4)+2*pad_val];
    
    % Crop down image to prey region using 'BoundingBox' output
    % crop image syntax: imcrop(im,[XMIN YMIN WIDTH HEIGHT])
    imBlob = imcrop(im,rect);
    
    % Crop down binary image to prey region
    imBlobBW = imcrop(imBW,rect);
    
    % Store prey data
    blob.sizePrey    = props2.Area;
    blob.xPrey       = props2.Centroid(1);
    blob.yPrey       = props2.Centroid(2);
    blob.thetaPrey   = props2.Orientation;
    blob.roi_blob    = rect;
    blob.xVals       = props2.PixelList(:,1);
    blob.yVals       = props2.PixelList(:,2);
    blob.MajorAxis   = props2.MajorAxisLength;
    blob.MinorAxis   = props2.MinorAxisLength;
    
    % Store cropped images
    blob.im      = imBlob;
    blob.BW      = imBlobBW;

% If exactly one blob ...
else 
    
% Get bounding box for prey
    perim = props2.BoundingBox;

    % Number of pixels that pad the blob
    pad_val = 5;
    
    % Rectangle for blob roi: [XMIN YMIN WIDTH HEIGHT]
    rect = [perim(1)-pad_val, perim(2)-pad_val,...
        perim(3)+ 2*pad_val, perim(4)+2*pad_val];
    
    % Crop down image to prey region using 'BoundingBox' output
    % crop image syntax: imcrop(im,[XMIN YMIN WIDTH HEIGHT])
    imBlob = imcrop(im,rect);
    
    % Crop down binary image to prey region
    imBlobBW = imcrop(imBW,rect);
    
    % Store prey data
    blob.sizePrey    = props2.Area;
    blob.xPrey       = props2.Centroid(1);
    blob.yPrey       = props2.Centroid(2);
    blob.thetaPrey   = props2.Orientation;
    blob.roi_blob    = rect;
    blob.xVals       = props2.PixelList(:,1);
    blob.yVals       = props2.PixelList(:,2);
    blob.MajorAxis   = props2.MajorAxisLength;
    blob.MinorAxis   = props2.MinorAxisLength;
    
    % Store cropped images
    blob.im      = imBlob;
    blob.BW      = imBlobBW;
    
end
