function acqMaster
% Manages the workflow for the acquisition of kinematics


%% Parameter values

% Extension for image files
p.nameSuffix = 'tif';

% Number of digits in batch name
p.num_digit_batch = 4;

% Number of digits for frame number in file name
p.num_digit_frame = 6;

% Number of frames to visualize at start of acquisition
p.numVis = 50;

% Max number of frames for creating the mean image
p.max_frames = 100;

% Frame rate (fps)
p.framerate = 1000;

% Size of calibration squares
p.sqr_size = 10e-2


%% Path definitions

% Matt's computer
if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')])) 
    % Directory root
    root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    
elseif ispc && ~isempty(dir(['C:\Users\Amberle_2\Documents']))
    root = 'D:\Amberle\Looming stimulus experiments';
    
else
    error('This computer is not recognized')
end

% To raw video files 
paths.rawvid = [root filesep 'Raw video'];

% To data files
paths.data = [root filesep 'Data'];

% To video thumbnails (for post-processing)
paths.thumb  = [root filesep 'Thumbnail video'];

% To calibration files
paths.cal = [root filesep 'Calibrations'];


%% Set up a new batch of experiments

% Get listing of batches of raw video
batchList = dir([paths.rawvid filesep 'B*']);

% Check for batches
if ~isempty(batchList)   
    
    % Loop thru batches
    for i = 1:length(batchList)
        % Get batch number from directory name
        batchNums(i,1) = str2num(batchList(i).name(end-p.num_digit_batch+1:end));
        
        % Store in batchList
        batchList(i).batchNum = batchNums(i,1);
    end
    
    % Set up inputdlg
    prompt={'New batch number:'};
    name='Batch setup';
    defaultanswer={num2str(max(batchNums)+1)};
    
    % Ask about batch
    answer=inputdlg(prompt,name,1,defaultanswer);
    
    % Get current batch
    batchNum = str2num(answer{1});
    
    % Check requested batch number against existing
    if max(batchNum==batchNums)
        errordlg(['Batch ' num2str(batchNum) ' already exists. You need to' ...
            ' delete existing batch folders if you want to use this' ...
            ' batch number.','']);
    end
    
else
    warning(['No batch directories in: ' paths.rawvid]);
    batchNum = 1;
end

% Define batch directory name
batchName = ['0000000000' num2str(batchNum)];
batchName = ['B' batchName((end-p.num_digit_batch+1):end)];

% Make data directory
mkdir([paths.data filesep batchName])

% Make thumbnail directory
mkdir([paths.thumb filesep batchName])

% Make video directory
mkdir([paths.rawvid filesep batchName])

clear answer batchList i prompt name numLines defaultanswer batchNums


%% Run (or load) calibration

% Promp for what to do about calibrtaions 
calChoice = questdlg('What calibration do you want to use?', ...
                         'Calibration', ...
                         'New', 'Previous', 'Cancel', 'New');

% If canceling
if strcmp(calChoice,'Cancel')
    return
    
% If using previous 
elseif strcmp(calChoice,'Previous')
    
    % Get batch listing
    a = dir([paths.cal filesep 'B*']);
    
    % Check input
    if isempty(a)
        error('No calibrations completed!')
    end
    
    % Step thru batches
    j = 1;
    for i = 1:length(a)
        if ~isempty([paths.cal filesep batchName filesep ' calibration data.mat'])
            bNum(j,1) = str2num(a(j).name((end-p.num_digit_batch+1):end));   
            bPath{j} = a(j).name;
            j = j + 1;
        end
    end
    
    % Check 
    if j==1
        error('No calibration yet completed')
    end
    
    % Last batch with calibration completed
    iBatch = bNum==max(bNum);
    
    % Name of previous batch
    calBatchName = bPath{iBatch};
    
    % Report which used
    disp(['Using batch ' num2str(bNum(iBatch)) ' calibration'])
    
    % Load 'cal' data
    load([paths.cal filesep calBatchName filesep 'calibration data.mat'])
    
% If creating a new calibration
elseif strcmp(calChoice,'New')
    
    % Name of calibration same as experiment batch name
    calBatchName = batchName;
    
    % Path to checkerboard video
    check_path = [paths.cal filesep calBatchName filesep 'Checkerboard video'];
    
    % Check for batch directory
    if isempty(dir([paths.cal filesep calBatchName]))
        mkdir([paths.cal filesep calBatchName])
        %error(['No calibration found for ' batchName])
    end
    
    % Check for checkboard video
    if isempty(dir([check_path filesep '*.' p.nameSuffix]))
        error(['Video tifs need to be saved in a folder called '...
               ' "Checkerboard video" in:' paths.cal filesep batchName])
    end   
    
    % Run calibration
    cal = runCalibration([paths.cal filesep calBatchName],p.sqr_size,...
                          p.num_digit_frame);   
end

% Save calibration data in current batch data directory
save([paths.data filesep batchName filesep ...
     [calBatchName 'calibration data.mat']],'cal');


%% Choose ROI

if isempty(dir([paths.cal filesep batchName filesep 'roi_data.mat']))

    % Get video information
    M = videoInfo([paths.cal filesep batchName filesep 'Checkerboard video'] ...
        ,p.num_digit_frame);
    
    % Create figure
    f = figure;
    
    % Read frame
    im = imread(M.path{1});
    
    % Undistort
    im = undistortImage(imadjust(im), cal.cameraParams,'OutputView','full');
    
    % Show frame
    warning off
    imshow(im)
    warning on
    title(['Choose ROI points.'])
        
    % Interactively find ROI
    h = impoly;
    roi_poly = wait(h);
    
    % Store results
    tmp = getPosition(h);
    roi.x = tmp(:,1);
    roi.y = tmp(:,2);

    % Show results
    delete(h)
    hold on
    plot(roi.x,roi.y,'r-')
    pause(1)
    hold off
    
    % Save ROI data
    save([paths.cal filesep batchName filesep 'roi_data.mat'],'roi')
    
    close
end




