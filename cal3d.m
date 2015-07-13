function cal3d(root)


%% Code execution

% Plot results of each step
show_steps = 1;


%% Parameters

% Starting number of calibration images
num_im = 20;

num_increment = 10;

% Size of checkerboard square (m)
sqr_size = 8.261e-3;


%% Define paths

% Root
%root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/3D Calibration';
%root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Sequences/Calibration/4mm';
if nargin<1
root = '/Volumes/SP UFD U3/Sequences/4mm';
end

% Path to Right cam
rcam_path = 'dorsal';

% Path to Left cam
lcam_path = 'lateral';


%% Get information about tiffs

% Gather info about the tiff files from right cam, save
if isempty(dir([root filesep 'Rfileinfo.mat']))
    R = info_tiffs(root,rcam_path);
    save([root filesep 'Rfileinfo.mat'],'R');
else
    % Load 'R'
    load([root filesep 'Rfileinfo.mat'])
end

% Gather info about the tiff files from left cam, save
if isempty(dir([root filesep 'Lfileinfo.mat']))
    L = info_tiffs(root,lcam_path);
    save([root filesep 'Lfileinfo.mat'],'L');
else
    % Load 'L'
    load([root filesep 'Lfileinfo.mat'])
end



%% 'Create calibration images' %------------------------------------

% Status variable
trying = 1;


% Report on start
disp(['Trying calibration with ' num2str(num_im) ' frames . . .'])

while trying


% Index for frames to use in calibration
iFiles = round(linspace(length(R.time).*.1,length(R.time).*.8,num_im));

% Gather listing of calibration images
for i = 1:length(iFiles)
   % Right cam
   Rfiles{i} = [root filesep rcam_path filesep R.name{i}];
   
   % Left cam
   Lfiles{i} = [root filesep lcam_path filesep L.name{i}];
    
end

% Detect the checkerboard corners in the images.
[imPoints,boardSize,imUsed] = detectCheckerboardPoints(Rfiles,Lfiles);

% Visualize this step
if show_steps
    figure
    subplot(1,2,1)
    warning off
    imshow(imread(Lfiles{1}))
    warning on
    if ~isempty(imPoints)
        hold on
        plot(imPoints(:, 1, 1, 2), imPoints(:, 2, 1, 2), '+g',...
            imPoints(1, 1, 1, 2), imPoints(1, 2, 1, 2), 'sy');
    end
    title('Left cam')
    
    subplot(1,2,2)
    warning off
    imshow(imread(Rfiles{1}))
    warning on
    if ~isempty(imPoints)
        hold on
        plot(imPoints(:, 1, 1, 1), imPoints(:, 2, 1, 1), '+g',...
            imPoints(1, 1, 1, 1), imPoints(1, 2, 1, 1), 'sy');
    end
    title('Right cam')
end

if isempty(imPoints)
    error('Failed to detect checkboards');
end

% Report on images
disp(['   ' num2str(sum(imUsed)) ' of ' num2str(length(imUsed)) ...
    ' images were usable'])

% Generate the world coordinates of the checkerboard corners in the
% pattern-centric coordinate system, with the upper-left corner at (0,0).
worldPoints = generateCheckerboardPoints(boardSize, sqr_size);

% Calibrate the camera
try
    stereoParams = estimateCameraParameters(imPoints,worldPoints, ...
        'EstimateSkew', true, 'EstimateTangentialDistortion', true, ...
        'NumRadialDistortionCoefficients', 3, 'WorldUnits', 'm');
    
    
    trying = 0;
   
catch
    % Report results
    beep;beep;
    disp('Calibration failed'); disp(' ')
    
    num_im = num_im + num_increment;
    
    if num_im > 50
        disp(' ')
        disp('All attempts failed');
        disp(' ')
        return
    else
        disp(['Trying again with ' num2str(num_im) ' frames . . .'])
    end
    stereoParams = [];
    
end


end



if show_steps
    figure;
    
    showReprojectionErrors(stereoParams);
end


% Report results
disp('Calibration completed !'); disp(' ')



% Save data
cal.Rfiles = Rfiles;
cal.Lfiles = Lfiles;
cal.stereoParam = stereoParams;
cal.imPoints = imPoints;
cal.boardSize = boardSize;
cal.imUsed = imUsed;
cal.worldPoints = worldPoints;

save([root filesep 'stereo calibration.mat'],'cal');
    


return



    if ~isempty(cameraParams)
        
        % View reprojection errors
        figure;
        subplot(2,2,1)
        showReprojectionErrors(cameraParams, 'BarGraph');
        title(cal_name)
        
        % Visualize pattern locations
        subplot(2,2,2)
        showExtrinsics(cameraParams, 'CameraCentric');
        view([0 -45])
        
        % Use the calibration data to remove effects of lens distortion.
        originalImage = imread(files{1});
        I = rgb2gray(originalImage);
        undistortedImage = undistortImage(I, cameraParams);
        
        subplot(2,2,3)
        %I = insertText(I, imagePoints(:,:,1), 1:size(imagePoints, 1));
        I = insertMarker(I, imagePoints(:,:,1), 'o', 'Color', 'red', 'Size', 5);
        warning off
        imshow(I)
        warning on
        title('Original image')
        
        subplot(2,2,4)
        warning off
        imshow(undistortedImage)
        warning on
        title('Corrected image')
    end



function A = info_tiffs(root,local_path)
% Gathers information about the tiff files in a diretcory. Assumes that the
% capture time of each frame is store in the filename, with millisecond
% precision (e.g. "Dorsal Cal 1.12-07-28-673.tif")

% Get name of tiff files
a = dir([root filesep local_path filesep '*.tif']);

% Store local path
A.local_path = local_path;

% Check contents
if isempty(a)
    error('No tif images in requested directory');
end

% Gather time stamp from filenames
for i = 1:length(a) 
    
    %info = imfinfo([root filesep local_path filesep a(i).name]);
    
    % Find periods/dashes in filename
    iStart = regexp(a(i).name,'[.]');
    iDash = regexp(a(i).name,'[-]');
    
    % Check
    if length(iStart)<2
        error('File name should have a "." and end in ".tif"')
        
    %elseif (iStart(2)-iStart(1))~=13
    %    error('Timestamp in filename should be 13 characters')
        
    elseif length(iDash)~=3
        error('Filename should have 3 dashes')
        
    end
    
    % Time in filename
    time_str = a(i).name((end-15):(end-4));
    
    % Convert to time vector
    A.time(i,:) = datevec(time_str,'HH-MM-SS-FFF');
    
    % Store filename
    A.name{i} = a(i).name;    
end

% Check for non-consecutive time code
if min(diff(etime(A.time(2:end,:),repmat(A.time(1,:),length(A.time(:,1))-1,1))))<0
    error('Files not in consecutive time order')
end


    
    
    
