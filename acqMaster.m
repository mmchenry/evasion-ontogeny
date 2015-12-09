function acqMaster(batchName,seqName)
% Manages the workflow for the acquisition of kinematics

if nargin < 2
    batchName   = '2015-11-16';
    seqName     = 'S01';
end


%% Parameter values

% Include calibration
includeCalibration = 0;

% Extension for image files
p.nameSuffix = 'jpg';

% Number of digits in batch name
p.num_digit_batch = 4;

% Number of digits for frame number in file name
p.num_digit_frame = 4;

% Number of frames to visualize at start of acquisition
p.numVis = 50;

% Max number of frames for creating the mean image
p.max_frames = 100;

% Frame rate (fps)
p.framerate = 250;

% Frames to skip initially in analysis
%p.skipFrame_start = 5;

% Size of calibration squares
p.sqr_size = 0.8e-2;

% The rate of rotation that qualifies as a significnat turn (rad/s)
p.alpha_thresh = 20;

% Period after stim where the fast start needs to be executed to count (s)
p.FS_period = [0.25 1.5];

% Tolerance of smoothing spline for angular position data
p.tolAng = 0.1;

% Determines motion from minimum area of blob for moving pixels (pix)
p.minBlobArea = 5e3;

% Duration (in min) to wait before looking for another sequence to analyze
waitDur = 0.5;

% Duration (in s) to pause between each loop when waiting
loopDur = 2;   


%% Path definitions

% % Matt's computer
% if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')])) 
%     % Directory root
%     %root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
%     root     = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';
%     vid_root = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';
%     
% else
%     error('This computer is not recognized')
% end

% % % Alberto's MacMini
path = fullfile('Users','alberto','Documents','GitHub-SourceTree');

if ~isempty(dir([filesep path]))
    % Directory root
    root     = '/Volumes/VisualPred/ZF_visuomotor';
    vid_root = '/Volumes/VisualPred/ZF_visuomotor';
else
    error('This computer is not recognized')
end

% % % Alberto's MacBook
% path = fullfile('Users','A_Soto','Documents','Home');
% 
% if ~isempty(dir([filesep path]))
%     % Directory root
%     root     = '/Volumes/Backup/ZF_visuomotor';
%     vid_root = '/Volumes/Backup/ZF_visuomotor';
% else
%     error('This computer is not recognized')
% end



% To raw video files 
paths.rawvid = [vid_root filesep 'Raw video'];

% To calibration video files 
paths.calvid = [vid_root filesep 'Calibration video'];

% To data files
paths.data = [root filesep 'Data'];

% To video thumbnails (for post-processing)
paths.thumb  = [root filesep 'Thumbnail video'];

% To calibration files
paths.cal = [root filesep 'Calibrations'];

% Log files
paths.log = [root filesep 'logFiles']; 


%% Run (or load) calibration

if includeCalibration && ...
   isempty(dir([paths.cal filesep batchName filesep 'calibration data.mat']))
    
    % Update status
    update_status(paths.log,batchName,'running calibration','calibrating')
    
    % Prompt for what to do about calibrtaions
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
        
        % Report which used
        disp(['Using batch ' num2str(bNum(iBatch)) ' calibration'])
        
        % Load 'cal' data
        load([paths.cal filesep bPath{iBatch} filesep 'calibration data.mat'])
        
        % Make data directory
        mkdir([paths.cal filesep batchName])
    
        % Copy over roi data 
        if ~isempty(dir([paths.cal filesep bPath{iBatch} filesep ...
                         'roi_data.mat']))
            copyfile([paths.cal filesep bPath{iBatch} filesep ...
                         'roi_data.mat'],[paths.cal filesep batchName filesep ...
                         'roi_data.mat']);
        end
        
    % If creating a new calibration
    elseif strcmp(calChoice,'New')
        
        % Name of calibration same as experiment batch name
        %calBatchName = batchName;
        
        % Path to checkerboard video
        check_path = [paths.calvid filesep batchName filesep 'Checkerboard video'];
        
        % Check for batch directory
        if isempty(dir([paths.cal filesep batchName]))
            mkdir([paths.cal filesep batchName])
            %error(['No calibration found for ' batchName])
        end
        
        % Check for checkboard video
        if isempty(dir([check_path filesep '*.' p.nameSuffix]))
            error(['Video tifs need to be saved in a folder called '...
                ' "Checkerboard video" in:' paths.cal filesep batchName])
        end
        
        % Run calibration
        cal = runCalibration([paths.calvid filesep batchName],...
                             [paths.cal filesep batchName],p.sqr_size,...
                             p.num_digit_frame);
        
        % Make data directory
        mkdir([paths.cal filesep batchName])
    end

    % Save calibration data in current batch data directory
    save([paths.cal filesep batchName filesep 'calibration data.mat'],'cal');
    
    % Update status
    update_status(paths.log,batchName,'On standby','standby')

elseif includeCalibration
    
    % Load 'cal' structure of calibration data
    load([paths.cal filesep batchName filesep 'calibration data.mat'])
    
else
    disp(' ');disp('    Running analysis without a calibration . . .')
    
end


%% Choose (or load) ROI

if isempty(dir([paths.cal filesep batchName filesep seqName filesep 'roi_data.mat']))
    
    %TODO: Modify this code for Alberto's proejct, once we start running
    %      calibrations
    
    if includeCalibration
        % Get video information        
        M = videoInfo([paths.calvid filesep batchName filesep 'Checkerboard video'] ...
            ,p.num_digit_frame);
        
        % Create figure
        f = figure;
        
        % Read frame
        im = imread(M.path{1});
          
        % Undistort
        im = undistortImage(imadjust(im), cal.cameraParams,'OutputView','full');
        
    else
        % Get listing of video frames
        aVid = dir([paths.rawvid filesep batchName filesep seqName ...
                    filesep 'exp*' p.nameSuffix]);
        
        % Read first frame
        im = imread([paths.rawvid filesep batchName filesep seqName ...
                    filesep aVid(1).name]);
        
        % Adjust contrast
        %im = imadjust(im);
    end
    
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
    
    % Make data directory
    if isempty(dir([paths.cal filesep batchName filesep seqName]))
        mkdir([paths.cal filesep batchName filesep seqName])
    end
    
    % Save ROI data
    save([paths.cal filesep batchName filesep seqName filesep 'roi_data.mat'],'roi')
    
    close
    
else
    
    % Load 'roi'
    load([paths.cal filesep batchName filesep seqName filesep 'roi_data.mat'])
    
end


%% Run analysis

% Modified to treat each sequence independently, to create new 'mean image'
% for every image sequence

% Set logical
continue_analysis = true;

% Initialize motion index
motion = [];

% Initialize previous
expName_prev = [];

% Name of current experiment
% NOTE: expName same as seqName, so got rid of for loop
expName = seqName;

% Directory for current video
vPath = [paths.rawvid filesep batchName filesep expName];

% Directory for current thumbnails
tPath = [paths.thumb filesep batchName filesep expName];

% Directory for current data
dPath = [paths.data filesep batchName filesep expName];

% Directory path for calibration
cPath = [paths.cal filesep batchName];

% Make thumbnail directory, if not present
if isempty(dir(tPath))
    mkdir(tPath)
end

% Make data directory, if not present
if isempty(dir(dPath))
    mkdir(dPath)
end

% NOTE: got rid of if statement to make a new mean image for every sequence

% Update status at command line
disp(' '); disp(' ')
disp(['---------- Analyzing ' expName ' in ' batchName ' ----------'])

% Analyze for midlines
motion = anaFrames('sequential',dPath,vPath,tPath,cPath,p,roi,includeCalibration);

% Analyze sequence data
bStats = anaSeq('prelim',dPath,tPath,cPath,p);

% Save parameters ('p')
save([dPath filesep 'parameters'],'p')


