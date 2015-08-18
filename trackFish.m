function trackFish(vPath,tPath,dPath,cPath,roi,p,overwrite)

% Acquire kinematics of predator and prey fish
% A previous version of this function was known as pred_prey_acq


%% Parameter values

% Extension for image files
%nameSuffix = 'tif';

% Number of digits for frame number in filename
%num_digit = 6;

% Initial visualizing predator acquisition
visSteps = 0;

% Number of predtaor frames to visualize at start of acquisition
%numVis = 50;

% Max number of frames for creating the mean image
maxFrames = 1000;

% Whether to invert the image
%invert = 1;


%% Get path of data file, load data

% Load filenames for frames
a = dir([vPath  filesep '*.' p.nameSuffix]);

if isempty(a)
    warning('No video frames found');
    return
end

% Get indicies for video frames
for i = 1:length(a)
    
    % Read filename
    frNum = str2num(a(i).name(end-p.num_digit_frame-length(p.nameSuffix):...
                    end-length(p.nameSuffix)-1));
                
    % Check ordering    
    if (i>1) && (frNum ~= p.frNums(i-1)+1)
        error('Frame numbers not consecutive')
    end
    
    % Store
    p.frNums(i,1) = frNum;
    p.filename{i} = a(i).name;
end

% Load calibration data ('cal')
load([cPath filesep 'calibration data.mat'])


%% Create or load mean image

% Look for mean image
a2 = dir([dPath filesep 'meanImage.tif']);

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
    h = waitbar(0,['Mean image: ' num2str(1)],'Name','Mean image',...
                     'CreateCancelBtn',...
                     'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vPath filesep p.filename{1}]);
    
    % Apply calibration to undistort image
    imCurr = undistortImage(imCurr, cal.cameraParams,'OutputView','full');
    
    imSum = double(imCurr);
    clear imCurr tmp
      
    % Loop through frames 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vPath  filesep p.filename{frIdx(i)}]);
        
         % Apply calibration to undistort image
         imCurr = undistortImage(imCurr, cal.cameraParams,'OutputView','full');
    
        % Sum pixel values
        imSum        = imSum + double(imCurr);
        clear tmp 
        
        % Update status bar
%         h = waitbar(i/length(frIdx),h,...
%             [num2str(i) ' of ' num2str(length(frIdx)) ' frames']);
        
        % Check for Cancel button press
        if getappdata(h,'canceling')
            close(h,'force')
            return
            
            % Otherwise, update status
        else
            waitbar(i/length(frIdx),h,...
                ['Mean image: ' num2str(i) ' of ' ...
                num2str(length(frIdx)) ' frames'])
        end
            
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
        
    end
    
    if strcmp(class(imCurr),'uint8')   
        % Calculate mean from sum image
        imMean = uint8(round(imSum./length(frIdx)));        
    elseif strcmp(class(imCurr),'uint16')
        % Calculate mean from sum image
        imMean = uint16(round(imSum./length(frIdx)));     
    end
    
    %imMean = imMean(:,:,1);
    
    % Write image to movie dir
    imwrite(imMean,[dPath filesep 'meanImage.tif'],'tif',...
            'Compression','none');
        
    % Close status bar
    close(h,'force')
    
    clear frIdx h i imSum imCurr
 
% Load mean image, if present
else
    imMean = imread([dPath filesep 'meanImage.tif']);  
end


%% Step through frames for position of predator

% Look for data file
a3 = dir([dPath filesep 'fish_coords.mat']);

% If data file not present . . .
if isempty(a3) || overwrite
    
    % Update status
    disp('      Creating thumbnails . . .')   
    
    % Status update
    if visSteps
        % Figure window
        f = figure;
        set(f,'DoubleBuffer','on')
    else
        hW = waitbar(0,'1','Name','Creating thumbnails',...
                     'CreateCancelBtn',...
                     'setappdata(gcbf,''canceling'',1)');
        setappdata(hW,'canceling',0)
    end
 
    % Loop through frames
    for i = 1:length(p.frNums)
        
        % Read full frame
        [imStart,cmap] = imread([vPath  filesep p.filename{i}]);
    
        % Apply calibration to undistort image
        imStart = undistortImage(imStart, cal.cameraParams,'OutputView','full');
        
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
        
        % Threshold value
        p.tVal = min([0.95 graythresh(im)+0.1]);
        imBW    = ~im2bw(im,p.tVal) & roiI;
        
        % Dilate im & get properties
        se    = strel('disk',6);
        imBW = imdilate(imBW,se);
        %imBW = imerode(imBW,se);
        LL    = bwlabel(imBW);
        props = regionprops(LL,'Centroid','Area');
        
        % Get peripheral shapes
        [B,L] = bwboundaries(imBW,'noholes');
        
        % Image showing intersection of roi boundary & fish
        imOverlap = ~roiI & imBW;
        
        % Store file info
        pd(i).frame = p.frNums(i);
        pd(i).filename = p.filename{i};
        
        % If no fish or touching a wall . . .
        if isempty(B) || (max(imOverlap(:))==1)
            if isempty(B)
                % Declare warning
                warning(['No fish found in ' p.filename{i}]);
            else
                % Declare warning
                warning(['Fish touching wall in ' p.filename{i}]);
            end
            
            % nan blob
            imBlob = nan;
            
            % Store nans          
            pd(i).xPerim = nan;
            pd(i).yPerim = nan;
            pd(i).roi_blob = nan;
          
        % If fish . . .
        else
            
            % Select blob with greatest periphery
            maxB = 0;
            for j = 1:length(B)
                if length(B{j}) > maxB
                    maxB = length(B{j});
                    perim = B{j};
                end
            end
            
            % Crop down to blob
            imBlob = imcrop(im,[min(perim(:,2)) min(perim(:,1)) ...
                range(perim(:,2)) range(perim(:,1))]);
                    
            % Enhance blob contrast
            %imBlob = adapthisteq(imBlob,'clipLimit',0.02,'Distribution','rayleigh');
            %imBlob = imadjust(imBlob);
            
            % Store data
            pd(i).xPerim = perim(:,2);
            pd(i).yPerim = perim(:,1);
            pd(i).roi_blob = [min(perim(:,2)) min(perim(:,1)) ...
                range(perim(:,2)) range(perim(:,1))];
        end
        
        % Save blob image
        save([tPath filesep p.filename{i}(1:(end-length(p.nameSuffix)-1))],'imBlob')
        
        % Visualize frames
        if visSteps
  
            figure(f)
            warning off
            subplot(3,1,1:2)
            imshow(imStart,cmap)
            hold on
            plot(pd(i).yPerim,pd(i).xPerim,'r-')
            plot([roi.x; roi.x(1)],[roi.y; roi.y(1)],'g-')
            title(['Frame ' num2str(p.frNums(i)) ' of ' ...
                num2str(p.frNums(end))])
            hold off
            subplot(3,1,3)
            imshow(imBlob)
            pause(.001)
            warning on
            clear im
            
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
        
        % Clear variables for next loop
        clear im  imBW  props imROI se x_roi y_roi maxB imBlob props LL       
    end
    
    % Save data
    save([dPath filesep 'fish_coords'],'pd')
    
    % Close status bar
    close(hW,'force')
    
else % If data file present . . .
    
    % Load 'pd' structure of predator coordinates
    disp('      Thumbnails previously acquired')
    
end
    


