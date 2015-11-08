function acqMaster(runBatch,batchNum)
% Manages the workflow for the acquisition of kinematics


%% Parameter values

% Wait for the stimulusPC to be initiated
waitForStim = 0;

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
p.minBlobArea = 20e3;

% Duration (in min) to wait before looking for another sequence to analyze
waitDur = 0.5;

% Duration (in s) to pause between each loop when waiting
loopDur = 2;   


%% Path definitions

% Matt's computer
if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')])) 
    % Directory root
    %root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    root = '/Users/mmchenry/Dropbox/Projects/Looming Project/';
    vid_root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    
% Photron (i.e. Maingear) computer
elseif ~isempty(dir('C:\Users\experimentalist'))
    root = 'C:\Users\experimentalist\Dropbox\Looming Project';
    vid_root = 'C:\Users\experimentalist\Desktop\';
    
% Amberle's computer
elseif ispc && ~isempty(dir(['C:\Users\Amberle_2\Documents']))
    root = 'D:\Amberle\Looming stimulus experiments';
    vid_root = root;
    
else
    error('This computer is not recognized')
end

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


%% Query about mode

if nargin<2
    % Query
    ButtonName1 = questdlg('Run experiments or analyze previous run?', ...
        'Task?', 'Experiments', 'Analysis', 'Experiments');
    
    % Check for just analysis
    if strcmp(ButtonName1,'Analysis')      
        % Set to false
        runBatch = 0;     
        
    % If running experiments
    elseif strcmp(ButtonName1,'Experiments')
        % Set to true
        runBatch = 1;
        
    else
        return
    end
    
    % Clear variables
    clear ButtonName2 ButtonName1
    
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
        if runBatch
            warning(['Creating new batch b/c no batch directories in: ' paths.data]);
            batchNum = 1;
        else
            error(['No batches to analyze in: ' paths.rawvid]);
        end
    end
    
    % Batch to analyze   
    if runBatch
        % Set up inputdlg
        prompt={'New batch number:'};
        name='Batch setup';
        defaultanswer={num2str(curr_batch+1)};
        
    else       
        % Set up inputdlg
        prompt={'Batch to analyze:'};
        name='Batch analysis';
        defaultanswer={num2str(curr_batch)};      
    end
    
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

    if runBatch
        % Make data directory
        mkdir([paths.data filesep batchName])
        
        % Make thumbnail directory
        mkdir([paths.thumb filesep batchName])
        
        % Make video directory
        mkdir([paths.rawvid filesep batchName])
    end
     
end %nargin<4

% Define batch directory name
batchName = ['0000000000' num2str(batchNum)];
batchName = ['B' batchName((end-p.num_digit_batch+1):end)];

clear batchNum

% Make log directory, if necessary
if isempty(dir([paths.log filesep batchName]))
    mkdir([paths.log filesep batchName]);
end


%% Run (or load) calibration

if isempty(dir([paths.cal filesep batchName filesep 'calibration data.mat']))
    
    % Update status
    update_status(paths.log,batchName,'running calibration','calibrating')
    
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
else
    
    % Load 'cal' structure of calibration data
    load([paths.cal filesep batchName filesep 'calibration data.mat'])
    
end


%% Choose (or load) ROI

if isempty(dir([paths.cal filesep batchName filesep 'roi_data.mat']))

    % Get video information

    M = videoInfo([paths.calvid filesep batchName filesep 'Checkerboard video'] ...
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
    
else
    
    % Load 'roi'
    load([paths.cal filesep batchName filesep 'roi_data.mat'])
    
end


%% Run analysis

if runBatch && waitForStim
    
    % Update status
    update_status(paths.log,batchName,'on standby','standby')
    
    % Read log
    lg = read_log(paths.log,batchName);
    
    % If StimPC is not ready 
    if ~strcmp(lg.status,'standby')
        % Loop until on standby
        while ~strcmp(lg.status,'standby')
            lg = read_log(paths.log,batchName);
            disp('    Waiting for StimPC . . .')
            pause(10);
        end
    end
end

% Set current sequence
if runBatch   
    currSeq = 0;
else
    currSeq = 1;
end
    
% Initial list of sequences
seqList = dir([paths.rawvid filesep batchName filesep 'C1S*']);
    
% Check list
if ~runBatch && isempty(seqList)
    warning('No sequences in the batch directory')
end
    
% Set logical
continue_analysis = true;

% Initialize motion index
motion = [];

% Initialize previous
expName_prev = [];
    
while continue_analysis
    
    % Wait for new sequence, if running a batch
    if runBatch        
        % Wait for new video
        currSeq = waitGame(paths,batchName,currSeq,p,waitDur,loopDur,expName_prev,motion);
        
        % Update status
        update_status(paths.log,batchName,'running analysis','analyzing',...
            expName_prev,motion)
        
        % Update list of sequences
        seqList = dir([paths.rawvid filesep batchName filesep 'C1S*']);
    end
    
    % Name of current experiment
    expName = seqList(currSeq).name;
    
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
    
    % If this is not the first recording and there is motion in prior movie
    if ~isempty(motion) && (motion==1)
        
        % Copy meanImage from prior
        copyfile([paths.data filesep batchName filesep ...
                  seqList(currSeq-1).name filesep 'meanImage.tif'],...
                 [dPath filesep 'meanImage.tif']);
        
    end
    
    % Update status at command line
    disp(' '); disp(' ')
    disp(['---------- Analyzing ' expName ', ' num2str(currSeq) ' of ' ...
        num2str(length(seqList)) ' in ' batchName ' ----------'])
    
    % Course pass on analysis
    motion = anaFrames('parallel',dPath,vPath,tPath,cPath,p,roi,...
                       p.skipFrame_start);
    
   %TODO: Figure out why motion is misdiagnosed
                   
    % Update status
    update_status(paths.log,batchName,...
            ['finished seq ' num2str(currSeq)],'analyzing',expName,motion);               
    
    % If fish moved . . .
    if motion
        % Update status
        %update_status(paths.log,batchName,'fish responded, starting FS analysis');
        
        % Analyze sequence data
        bStats = anaSeq('prelim',dPath,tPath,cPath,p);
        
        % Analyze every frame of the fast start
        if ~isnan(bStats.tFS(1))
            anaFrames('sequential',dPath,vPath,tPath,cPath,p,roi,0,...
                bStats.frFS(1),bStats.frFS(2));
        end
        
    % If no motion . . .
    else
        % Update status
        %update_status(paths.log,batchName,'fish motionless, rerunning expt');
        
        %update_status(paths.log,batchName,'fish motionless',...
        %    'running analysis',0)
    end
    
    % Update status
    update_status(paths.log,batchName,'on standby',...
            'standby',expName,motion);
        
    expName_prev = expName;
    
    % If all sequences are analyzed
    if ~runBatch && (currSeq==length(seqList))
        
        % Update status
        disp(['Completed all ' num2str(currSeq) ' sequences in ' batchName])
        
        % End execution
        continue_analysis = false;
        
    % Otherwise, continue analysis
    elseif ~runBatch
        currSeq = currSeq + 1;
        
    end
end


function currSeq = waitGame(paths,batchName,currSeq,p,waitDur,loopDur,expName,motion)

% Update status
update_status(paths.log,batchName,'waiting for new video','standby',expName,motion)

% Set logical for waiting loop
continue_waiting = true;

% Waiting loop
while continue_waiting
    
    % Get latest log file
    lg = read_log(paths.log,batchName);
    
    %if strcmp(lg.PC,'stimPC') && ~isempty(lg.last_seq)
        
    % Survey list of sequences recorded
    seqList = dir([paths.rawvid filesep batchName filesep 'C1S*']);
    
    % If there is a new sequence . . .
    if length(seqList)>currSeq
        
        % Set new current sequence
        currSeq = currSeq + 1;
        
        % Loop to make sure files aren't being written
        while true
            % List files
            fileList1 = dir([paths.rawvid filesep batchName filesep ...
                seqList(currSeq).name filesep '*.' p.nameSuffix]);
            
            % Wait before checking again
            pause(10)
            
            % List files
            fileList2 = dir([paths.rawvid filesep batchName filesep ...
                seqList(currSeq).name filesep '*.' p.nameSuffix]);
            
            % Break, if no new files
            if length(fileList1)==length(fileList2)
                break
            end
        end
        
        % Exit wait loop
        continue_waiting = false;
        
    % Otherwise (no new sequence) . . .
    else
        
%         % Initialize wait bar
%         hW = waitbar(0,'1','Name','Analysis paused',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
%         setappdata(hW,'canceling',0)
        
        disp('    Waiting for new video . . .')
        
        % Wait loop
        for i = 1:ceil(waitDur*60/loopDur)
            
%             % Check for Cancel button press
%             if getappdata(hW,'canceling')
%                 
%                 % Kill the waitbar
%                 close(hW,'force')
%                 
%                 % Stop all code
%                 return
%                 
%                 % Otherwise, update waitbar status
%             else
%                 waitbar(i/ceil(waitDur*60/loopDur),hW,...
%                     ['Waiting: ' sprintf('%6.1f',i*loopDur/60) ' of ' ...
%                     num2str(waitDur) ' min for seq ' num2str(currSeq)])
%             end
            
            % Pause to update status
            pause(loopDur)
        end
        
        % Kill waitbar for next pass thru loop
        %close(hW,'force')
    end
end


function update_status(log_path,batchName,event_txt,status_txt,...
                       lastseq,motion)
% Save a log file to update status 

% Formulate date and time strings
date_string = strcat(datestr(clock,'dd-mmm-yyyy-HH'),'_',...
                     datestr(clock, 'MM'),'_',datestr(clock,'ss'),'_',...
                     datestr(clock,'FFF'));
time_str = strcat([datestr(clock,'HH'),':',datestr(clock, 'MM'),'.',...
                   datestr(clock,'ss'),'.',datestr(clock,'FFF') ]);

% Store in log structure
lg.PC      = 'PhotronPC';
lg.status  = status_txt;
lg.msg     = strcat([time_str '  PhotronPC: ' event_txt]);

% Log motion
if nargin<5
    lg.lastseq = [];
else
    % Last sequence analyzed
    
    lg.lastseq = lastseq;
end

if nargin<6
    lg.motion = [];
else
    lg.motion = motion;
end

% Save
save([log_path filesep batchName filesep 'L' date_string '_PhotronPC'],'lg')


function log_out = read_log(log_path,batchName)

% Get log files
a = dir([log_path filesep batchName filesep 'L*']);

% Set default as 'standby'
log_out.status = [];

% Loop trhu log files
if ~isempty(a)
    for i = 1:length(a)
        % Load log structure
        load([log_path filesep batchName filesep a(i).name])
        
        if strcmp(lg.PC,'StimPC')
            log_out = lg;
        end       
    end
end



