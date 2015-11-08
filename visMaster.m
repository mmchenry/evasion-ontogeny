function visMaster(batchNum)
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

% Frames to skip initially in analysis
p.skipFrame_start = 5;

% Size of calibration squares
p.sqr_size = 0.8e-2;

% The rate of rotation that qualifies as a significnat turn (rad/s)
p.alpha_thresh = 20;

% Period after stim where the fast start needs to be executed to count (s)
p.FS_period = [0.25 1.5];

% Tolerance of smoothing spline for angular position data
p.tolAng = 0.1;

% Determines motion from minimum area of blob for moving pixels (pix)
p.minBlobArea = 50e3;

% Duration (in min) to wait before looking for another sequence to analyze
waitDur = 20;

% Duration (in s) to pause between each loop when waiting
loopDur = 2;   


%% Path definitions

% Matt's computer
if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')])) 
    % Directory root
    root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    vid_root = root;
    
% Photron (i.e. Maingear) computer
elseif ~isempty(dir('C:\Users\experimentalist'))
    root = 'C:\Users\experimentalist\Dropbox\Looming Project';
    vid_root = 'C:\Users\experimentalist\Video\';
    
% Amberle's computer
elseif ispc && ~isempty(dir(['C:\Users\Amberle_2\Documents']))
    root = 'D:\Amberle\Looming stimulus experiments';
    vid_root = root;
    
else
    error('This computer is not recognized')
end

% To raw video files 
paths.rawvid = [vid_root filesep 'Raw video'];

% To data files
paths.data = [root filesep 'Data'];

% To video thumbnails (for post-processing)
paths.thumb  = [root filesep 'Thumbnail video'];

% To calibration files
paths.cal = [root filesep 'Calibrations'];

%% Query about mode

if nargin<1
    
    % Get listing of batches of data
    batchList = dir([paths.rawvid filesep 'B*']);
    
    % If there are batches . . .
    if ~isempty(batchList)
        % Loop thru batches
        for i = 1:length(batchList)
            % Get batch number from directory name
            batchNums(i,1) = str2num(batchList(i).name(end-p.num_digit_batch+1:end));
            
            % Store in batchList
            %batchList(i).batchNum = batchNums(i,1);
        end
        
        % Current batch is the max number
        curr_batch = max(batchNums);
        
        % If no batches . . .
    else 
        error(['No batches to analyze in: ' paths.rawvid]);
    end
       
    % Set up inputdlg
    prompt={'Batch to analyze:'};
    name='Batch analysis';
    defaultanswer={num2str(curr_batch)};
    
    % Ask about batch
    answer=inputdlg(prompt,name,1,defaultanswer);
    
    % Check answer
    if isempty(answer)
        return
    end
    
    % Get current batch
    batchNum = str2num(answer{1});
    
    % Check requested batch number against existing
    if runBatch && max(batchNum==batchNums)
        warning(['Batch ' num2str(batchNum) ' already exists.','']);
    end
    
    clear curr_batch answer batchList i prompt name numLines defaultanswer batchNums
     
end %nargin<4

% Define batch directory name
batchName = ['0000000000' num2str(batchNum)];
batchName = ['B' batchName((end-p.num_digit_batch+1):end)];

clear batchNum


%% Load calibration

    % Load 'cal' structure of calibration data
    load([paths.cal filesep batchName filesep ...
        ['calibration data.mat']])
    



%% Load ROI

    
    % Load 'roi'
    load([paths.cal filesep batchName filesep 'roi_data.mat'])
    


%% Interactively visualize results for one sequence
 
    % Play the data
    playData(paths,p,batchName,cal,roi) 


