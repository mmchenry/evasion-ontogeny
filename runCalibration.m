function cal = runCalibration(calVid,calDir,sqr_size,num_digit)
% Conducts calibration for 3D kinematic experiments that use GoPro cameras
% Requires Matlab 2014b (or later), Computer Vision Toolbox and Image
% Processing Toolbox.


%% Parameter values

% Number of calibration images
numCalIm = 20;

% Difference in number of images between calibration runs  
num_increment = 10;

% Write results of calibration to disk
write_results = 1;

% Initial number of frames to be used for lens calibration
num_fr0 = 30;

% Get video information
M = videoInfo([calVid filesep 'Checkerboard video'],num_digit);
    

%% Creating calibration images

% If frames already exist . . .
if ~isempty(dir([calDir filesep 'frames']))
    
    % Prompt for what to do
    button = questdlg('Recreate frames to be analyzed?','Calibration','Yes',...
                      'No','Cancel','No');
    
    % Parse response
    if strcmp(button,'Yes')
        
        % Delete existing
        delete([calDir filesep 'frames' filesep '*.tif'])
        delete([calDir filesep 'frames_roi' filesep '*.tif'])
        
        % Set logical
        create_im = 1;
        
    elseif strcmp(button,'No')
        
        % Set logical
        create_im = 0;
        
    else
        return
    end
    
    clear button
  
% If no frames . . .
else
    create_im = 1;
end

% If creating images . . .
if create_im
    
    % Make directories
    [success,message,id] = mkdir(calDir,'frames');
    
    % Load frame numbers('ana_num'), if defined
    if ~isempty(dir([calDir filesep 'frame numbers.mat']))
        load([calDir filesep 'frame numbers.mat'])
        disp('Using frames from frames.mat file . . .')
        fr_num = ana_num+1;
        clear ana_num
    end
    
    fr_interval = round(length(M.frNums)/50);
    
    iFrame = fr_interval;
    
    % Initialize index
    k = 1;
    
    % Create figure window
    f = figure;
    
    % Loop thru calibration images
    %for i = 1:numCalIm
     while true
         
        % Read frame
        im = imread(M.path{iFrame});
             
        % Adjust grayscale values and convert to double
        im     = (imadjust(im));
        
        warning off
        imshow(im)
        title(['(y)es or (n)o? ' num2str(k-1) ' selected'])
        warning on
        
        ch = getkey;
        
        % If yes
        if ch==121
            
            % Frame string
            fr_str = ['000000' num2str(M.frNums(iFrame))];
            
            % image path
            im_path = [calDir filesep 'frames' filesep 'frame ' ...
                fr_str(end-5:end) '.tif'];
            
            % Write image files
            %imwrite(im,im_path,'TIFF');
            copyfile(M.path{iFrame},im_path)
            
            % Save paths
            %im_file.path{k} = im_path;
            
            % Advance to next frame
            iFrame = iFrame + fr_interval;
            
            if k==numCalIm
                close(f)
                break
            else
                k = k + 1;
            end
            
        % If no
        elseif ch==110
            
            iFrame = iFrame + fr_interval;
            
        else
            warning('keystroke not recognized')
            
        end
        
        % Start over, if exceeding end of video
        if iFrame > length(M.frNums)            
            iFrame = round(1.5*fr_interval);
        end
            
     end
    
end
         
% Survey existing frames (allows manual deleting)
a = dir([calDir filesep 'frames' filesep 'frame*']);

% Check contents
if isempty(a)
    error(['You need to create calibration frames in '...
        calDir filesep 'frames'])
end

% Loop thru calibration images
for i = 1:length(a)
    
    % Save paths
    im_file.name{i} = a(i).name;
    im_file.path{i} = [calDir filesep 'frames' filesep a(i).name];
end


%% Select ROI

if isempty(dir([calDir filesep 'frames_roi']))
    
    % Make directory
    mkdir([calDir filesep 'frames_roi']);

    % Create figure
    f = figure;
    
    % Loop thru images
    for i = 1:length(im_file.path)
        
        % Read frame
        im = imread(im_file.path{i});
        title(['Choose ROI points. ' num2str(i) ' of ' ...
               num2str(length(im_file.path))])
        
        % Show frame
        warning off
        imshow(imadjust(im))
        warning on
        
        % Interactively find ROI
        h = impoly;
        roi_poly = wait(h);
        
        % Return image within roi
        im2 = give_roi(im,roi_poly,'w');
        
        % Enhance contrast of cropped image
        im2 = imadjust(im2);
        
        % Write image file
        imwrite(im2,[calDir filesep 'frames_roi' filesep im_file.name{i}],'TIFF')
    end
    
    close
end

% Survey existing frames in roi folder (allows manual deleting)
a = dir([calDir filesep 'frames_roi' filesep 'frame*']);

% Check contents
if isempty(a)
    error(['You need to create calibration frames in '...
        calDir filesep 'frames_roi'])
end

% Loop thru calibration images
for i = 1:length(a)  
    % Save paths
    im_file.name_roi{i} = a(i).name;
    im_file.path_roi{i} = [calDir filesep 'frames' filesep a(i).name];
end


%% Run calibration

% Update status
disp(' Running calibration . . .')

% Detect the checkerboard corners in the images.
[imPts,brdSize,imUsed] = detectCheckerboardPoints(im_file.path_roi);

% Remove unused images
[im_file.path_roi,imUsed] = remove_unused(im_file.path_roi,imUsed);

% Create checkerboard points
worldPoints = generateCheckerboardPoints(brdSize, sqr_size);

% Find camera parameters
[cameraParam, imUsed, estErrors] = estimateCameraParameters(...
    imPts, worldPoints, ...
    'EstimateSkew', true, ...
    'EstimateTangentialDistortion', true, ...
    'NumRadialDistortionCoefficients', 3, ...
    'WorldUnits', 'm');

% Remove unused images
[im_file.path_roi,imUsed] = remove_unused(im_file.path_roi,imUsed);

% Translate detected points back into the original image coordinates
%refPts = bsxfun(@plus, imPts, mean(newOrigin{i}(imUsed),2));

% Calcualte mean extrinsic propterties of cameras
[R, t] = calc_extrinsics(imPts, worldPoints, imUsed, cameraParam);

% Calculate camera matrices from mean extrinsics
camMatrix = cameraMatrix(cameraParam, R, t);

% Write results of indivudal cameras to disk
write_cal_results(calDir, im_file.path_roi, imUsed, imPts, cameraParam);

% Save calibration data for camera A
cal.cameraParams  = cameraParam;
cal.R             = R;
cal.t             = t;
cal.camMatrix     = camMatrix;
cal.imPts         = imPts;
%
% Test calibration
im1 = imread(im_file.path{1});
[im, newOrigin] = undistortImage(im1, cameraParam,'OutputView','full');

if sum(im(:))==0
    error('Calibration failed.  Try deleting some of the calibration frames');
end

figure;
subplot(1,2,1)
imshow(im1)
title('Raw video frame')

subplot(1,2,2)
imshow(im)
title('Undistorted version')

% Save calibration data
save([calDir filesep 'calibration data.mat'],'cal')



function  im = give_roi(im,roi,field_clr)
% Return image in roi, with adjusted contrast

if nargin < 3
    field_clr = 'k';
end

% Define mask
BW = poly2mask(roi(:,1),roi(:,2),size(im,1),size(im,2)); 

% Update image, for verification
if strcmp(field_clr,'k')
    im(~BW) = 0;
else
    im(~BW) = 255;
end

% Adjust contrast
im = imadjust(im);


function write_cal_results(results_path, files, imUsed, imagePoints, cameraParams)
% Writes calibrtaion results to disk

% Clear or create directory
clear_dir(results_path,'Calibration images')

% Initiate index 
j = 1;

% Write images used for calibration
for i = 1:length(imUsed)
    if imUsed(i)        
        % Read image
        originalImage = imread(files{i});
        
        % Convert to gray, if necessary
        if size(originalImage,3)>1
            I = rgb2gray(originalImage);
        else
            I = originalImage;
        end
        
        % Add markers
        I = insertMarker(I, imagePoints(:,:,i), 'o', 'Color', 'green', 'Size', 5);
        I = insertMarker(I, imagePoints(1,:,i), 's', 'Color', 'yellow', 'Size', 8);
        
        % Frame string
        fr_str = ['0000' num2str(j)];
        
        % Write frame
        imwrite(uint8(I),[results_path filesep 'Calibration images' filesep 'image ' ...
            fr_str(end-3:end) '.jpeg'],'JPEG');
        
        % Add to index
        j = j + 1;
    end
end

% Visualize pattern locations
figure;
subplot(1,2,1)
showExtrinsics(cameraParams, 'CameraCentric');
view([0 -45])

% View reprojection errors
subplot(1,2,2)
showReprojectionErrors(cameraParams, 'BarGraph');
title('Single camera calibration')

% Capture graphs
I = getframe(gcf);
%close

% Write frame
imwrite(I.cdata,[results_path filesep 'calibration graph.jpeg'],'JPEG');


function  [R, t] = calc_extrinsics(imPts, worldPoints, imUsed, cameraParam)
% Calculate camera extrinsics for each frame, find mean
j = 1;

for i = 1:sum(imUsed)
    if imUsed(i)
        [R, t] = extrinsics(imPts(:,:,i), worldPoints, cameraParam);
        Rs(:,:,j)  = R;
        ts(j,:)    = t;
        j = j + 1;
    end
end

% Assuming fixed camera, take mean position
R = mean(Rs, 3);
t = mean(ts,1);


function [files,imUsed] = remove_unused(files,imagesUsed)
% Removes unused images from cell array

k = 1;
for j = 1:length(imagesUsed)
    if imagesUsed(j)
        imUsed(j) = imagesUsed(j);
        tmp{k} = files{j};
        k = k + 1;
    end
end
files = tmp;


function meanErrorsPerImage = computeMeanError(this)
% Code from 'showReproductionErrors' to calculate errors in calibration
errors = hypot(this.ReprojectionErrors(:, 1, :), ...
                this.ReprojectionErrors(:, 2, :));
meanErrorsPerImage = squeeze(mean(errors, 1));


function clear_dir(im_path,dir_name)
% Either removes existing images (if present), or creates empty diretory

% Delete image files within existing directory
if ~isempty(dir([im_path filesep dir_name]))
    % Delete image files within
    delete([im_path filesep dir_name filesep '*.tif'])
    delete([im_path filesep dir_name filesep '*.jpeg'])
else
    % Make directory
    [success,message,id] = mkdir(im_path,dir_name);
end


function pathh = exp_path(pred_sp,prey_sp,age,expt_type)
% Returns diretcory structure for certain experiments

pathh = [pred_sp ' pred' filesep ...
         prey_sp ' prey' filesep ...
         num2str(age) ' dpf' filesep ...
         expt_type];
  

function copy_file(file_source,file_dest)
% Copies files, reports results, deletes source 
[success,message,messageid] = copyfile(file_source,file_dest);
    
% Report results, check for failure
if success==1
    disp(' ');
    disp('Copy file completed:')
    disp(['From: ' file_source])
    disp(['To: ' file_dest])
else
    disp(['Failed file copy from ' file_source ' to ' file_dest]);
    disp(' ')
    disp(' ')
    error(message)
end

