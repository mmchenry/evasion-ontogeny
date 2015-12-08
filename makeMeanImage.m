function [imMin,imMean] = makeMeanImage(dPath,vPath,cal,B)


%% Parameters

% Max number of frames for creating the mean image
maxFrames = 1000;

% total number of frames in 'B'
frTot = length(B);


%% Create or load mean image

% Look for mean image
a2 = dir([dPath filesep 'meanImage.tif']);

% Calculate mean image if it does not exist
if isempty(a2)
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if frTot > maxFrames
        dframe = floor(frTot/maxFrames);
        frIdx = 1:dframe:frTot;
        clear dframe
    else
        frIdx = 1:frTot;
    end
    
    % Create waitbar
    h = waitbar(0,['Mean image: ' num2str(1)],'Name','Mean image',...
                     'CreateCancelBtn',...
                     'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vPath filesep B(1).filename]);
    
    % Start sum image
    imSum = double(imCurr.*0);
    
    % Min image (hilights?) 
    if strcmp(class(imCurr),'uint8')
        imMin = imCurr.*0 + 255; 
        
    elseif strcmp(class(imCurr),'uint16')
        imMin = imCurr.*0 + 2.^16;
        
    else
        error('Code only supports 16 and 8 bit images')
        
    end
    
    clear tmp
      
    % Loop through frames (iteratively updates 'imSum' and 'imMin') 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vPath  filesep B(i).filename]);
        
         % Apply calibration to undistort image
         %imCurr = undistortImage(imCurr, cal.cameraParams,'OutputView','full');
    
        % Sum pixel values (cumulative sum)
        imSum  = imSum + double(imCurr);
       
        % Update min image
        tMin(:,:,1) = imMin;
        tMin(:,:,2) = imCurr;
        imMin = min(tMin,[],3);
        
        % Clear for next
        clear tmp 
        
        % Check for Cancel button press
        if getappdata(h,'canceling')
            close(h,'force')
            error('Execution stopped by user');
            
            % Otherwise, update status
        else
            waitbar(i/length(frIdx),h,...
                ['Mean image: ' num2str(i) ' of ' ...
                num2str(length(frIdx)) ' frames'])
        end       
    end
    
    if strcmp(class(imCurr),'uint8')
        % Calculate mean from sum image
        imMean = uint8(round(imSum./length(frIdx)));
    elseif strcmp(class(imCurr),'uint16')
        % Calculate mean from sum image
        imMean = uint16(round(imSum./length(frIdx)));
    end
    
    if ~isempty(cal)
        % Apply calibration to undistort images
        imMean = undistortImage(imMean, cal.cameraParams,'OutputView','full');
        imMin  = undistortImage(imMin, cal.cameraParams,'OutputView','full');
    end
    
    % Write mean image to movie dir
    imwrite(imMean,[dPath filesep 'meanImage.tif'],'tif','Compression','none');
        
    % Write min image to movie dir
    imwrite(imMin,[dPath filesep 'minImage.tif'],'tif','Compression','none');
    
    % Close status bar
    close(h,'force')
    
    clear frIdx h i imSum imCurr
 
% Load images, if present
else
    imMean = imread([dPath filesep 'meanImage.tif']);  
    imMin  = imread([dPath filesep 'minImage.tif']); 
end




