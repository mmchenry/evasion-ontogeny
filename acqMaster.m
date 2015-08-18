function acqMaster(runAnalysis,setBatch,visData,batchNum,overwrite)
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
p.sqr_size = 0.8e-2;

% Duration (in min) to wait before looking for another sequence to analyze
waitDur = 20;

% Default for overwriting data
if nargin < 5
    overwrite = 0;
end


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


%% Query about mode

if nargin<4
    % Query
    ButtonName1 = questdlg('Run experiments, analyze, or visualize?', ...
        'Task?', 'Experiments', 'Analysis', 'Visualize', 'Experiments');
    
    % Check for just analysis
    if strcmp(ButtonName1,'Analysis')
        % Set to true
        runAnalysis = 1;
        
        % Set to false
        setBatch = 0;
        
        % Set to true
        visData = 1;
        
        % If running experiments
    elseif strcmp(ButtonName1,'Experiments')
        % Set to false
        runAnalysis = 0;
        
        % Set to true
        setBatch = 1;
        
        % Set to true
        visData = 0;
        
        % If visualizing results
    elseif strcmp(ButtonName1,'Visualize')
        % Set to false
        runAnalysis = 0;
        
        % Set to false
        setBatch = 0;
        
        % Set to true
        visData = 1;
        
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
        if runAnalysis
            error(['No batches to analyze in: ' paths.data]);
        else
            warning(['Creating new batch b/c no batch directories in: ' paths.data]);
            batchNum = 1;
        end
        %curr_batch = 1;
    end
    
    % Batch to anlayze
    if runAnalysis || visData
        % Set up inputdlg
        prompt={'Batch to analyze:'};
        name='Batch analysis';
        defaultanswer={num2str(curr_batch)};
        
    elseif setBatch
        % Set up inputdlg
        prompt={'New batch number:'};
        name='Batch setup';
        defaultanswer={num2str(curr_batch+1)};
        
    end
    
    % Ask about batch
    answer=inputdlg(prompt,name,1,defaultanswer);
    
    % Check answer
    if isempty(answer)
        return
    end
    
    % Get current batch
    batchNum = str2num(answer{1});
    
    % Get number of experiments
    %exptNum = str2num(answer{2});
    
    % Check requested batch number against existing
    if setBatch && max(batchNum==batchNums)
        warning(['Batch ' num2str(batchNum) ' already exists.','']);
    end
    
    clear curr_batch answer batchList i prompt name numLines defaultanswer batchNums

    if setBatch
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


%% Run (or load) calibration

if isempty(dir([paths.cal filesep batchName filesep ['calibration data.mat']]))
    
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
        
        % Make data directory
        mkdir([paths.cal filesep batchName])
    end
    
    
    
    % Save calibration data in current batch data directory
    save([paths.cal filesep batchName filesep ...
        ['calibration data.mat']],'cal');
else
    
    % Load 'cal' structure of calibration data
    load([paths.cal filesep batchName filesep ...
        ['calibration data.mat']])
    
end


%% Choose (or load) ROI

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
    
else
    
    % Load 'roi'
    load([paths.cal filesep batchName filesep 'roi_data.mat'])
    
end


%% Run experiments

% If running experiments 
if setBatch
    
    % Duration (in s) to pause between each loop when waiting
    loopDur = 2;    
    
    % Set logical
    continue_analysis = true;
    
    % Initialize current sequence
    currSeq = 0;
    
    while continue_analysis
        
        % Set logical for waiting loop
        continue_waiting = true;
        
        % Waiting loop
        while continue_waiting
            
            % Survey list of sequences
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
                    pause(15)
                    
                    % List files
                    fileList2 = dir([paths.rawvid filesep batchName filesep ...
                        seqList(currSeq).name filesep '*.' p.nameSuffix]);
                    
                    % Break, if no new files
                    if length(fileList1)==length(fileList2)
                        clear fileList1 fileList2
                        break
                    end
                end
                
                % Exit wait loop
                continue_waiting = false;
                
            % Otherwise (no new sequence) . . .
            else
                
                % Initialize wait bar
                hW = waitbar(0,'1','Name','Analysis paused',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(hW,'canceling',0)
                
                % Wait loop
                for i = 1:ceil(waitDur*60/loopDur)
                    
                    % Check for Cancel button press
                    if getappdata(hW,'canceling')
                        
                        % Kill the waitbar
                        close(hW,'force')
                        
                        % Stop all code
                        return
                        
                        % Otherwise, update waitbar status
                    else
                        waitbar(i/ceil(waitDur*60/loopDur),hW,...
                            ['Waiting: ' sprintf('%6.1f',i*loopDur/60) ' of ' ...
                            num2str(waitDur) ' min'])
                    end
                    
                    % Pause to update status
                    pause(loopDur)
                end
                
                % Kill waitbar for next pass thru loop
                close(hW,'force')
            end
        end
        clear continue_waiting        
         
        % Name of current experiment
        expName = seqList(currSeq).name;
        
        % Directory for current video
        currVid = [paths.rawvid filesep batchName filesep expName];
        
        % Directory for current calibration
        currCal = [paths.cal filesep batchName];
        
        % Directory for current thumbnails
        currThumb = [paths.thumb filesep batchName filesep expName];
        
        % Make directory, if not present
        if isempty(dir(currThumb))
            mkdir(currThumb)
        end
        
        % Directory for current data
        currData = [paths.data filesep batchName filesep expName];
         
        % Make directory, if not present
        if isempty(dir(currData))
            mkdir(currData)
        end
        
        % Update status
        disp(' '); disp(' ')
        disp(['---------- Analyzing ' expName ', ' num2str(currSeq) ' of ' ...
            num2str(length(seqList)) ' in ' batchName ' ----------'])
        
        % Track fish for the current sequence
        trackFish(currVid,currThumb,currData,currCal,roi,p,overwrite);
        
        % Analyze blobs for the current sequence
        anaBlobs(currVid,currThumb,currData,p,overwrite);

        clear seqList currData currCal currThumb expName 
    end
    
    clear loopDur 
end


%% Run analysis

if runAnalysis
    
    % List of sequences
    seqList = dir([paths.rawvid filesep batchName filesep 'C1S*']);
    
    % Check list
    if isempty(seqList)
        warning('No sequences in the batch directory')
    end
    
    % Set current sequence
    currSeq = 1;
    
    % Set logical
    continue_analysis = true;
    
    while continue_analysis
        
        % Name of current experiment
        expName = seqList(currSeq).name;
        
        % Directory for current video
        currVid = [paths.rawvid filesep batchName filesep expName];
        
        % Directory for current thumbnails
        currThumb = [paths.thumb filesep batchName filesep expName];
        
        % Make directory, if not present
        if isempty(dir(currThumb))
            mkdir(currThumb)
        end
        
        % Directory for current data
        currData = [paths.data filesep batchName filesep expName];
        
        % Directory for current calibration
        currCal = [paths.cal filesep batchName];
        
        % Make directory, if not present
        if isempty(dir(currData))
            mkdir(currData)
        end
        
        % Update status
        disp(' '); disp(' ')
        disp(['---------- Analyzing ' expName ', ' num2str(currSeq) ' of ' ...
            num2str(length(seqList)) ' in ' batchName ' ----------'])
        
        % Track fish for the current sequence
        trackFish(currVid,currThumb,currData,currCal,roi,p,overwrite);
        
        % Analyze blobs for the current sequence
        anaBlobs(currVid,currThumb,currData,p,overwrite);
        
        % Analyze sequence data
        %anaSeq(currVid,currThumb,currData,currCal,p,overwrite);
        
        % If running post-hoc analysis
        if runAnalysis
            % If all sequences are analyzed
            if currSeq==length(seqList)
                
                % Update status
                disp(['Completed all ' num2str(currSeq) ' sequences in ' batchName])
                
                % End execution
                continue_analysis = false;
                
                % Otherwise, continue analysis
            else
                currSeq = currSeq + 1;
            end
            
        % If running experiments . . .
        elseif setBatch
            
            % Set loop logical
            continue_waiting = true;
            
            % Waiting loop
            while continue_waiting
                
                % Survey list of sequences
                seqList = dir([paths.rawvid filesep batchName filesep 'C1S*']);
                
                % Identify if there is a new sequence
                if length(seqList)>currSeq
                    
                    % Set new current sequence
                    currSeq = currSeq + 1;
                    
                    % Loop to make sure files aren't being written
                    while true
                        % List files
                        fileList1 = dir([paths.rawvid filesep batchName filesep ...
                            seqList(currSeq).name filesep '*.' p.nameSuffix]);
                        
                        % Wait before checking again
                        pause(15)
                        
                        % List files
                        fileList2 = dir([paths.rawvid filesep batchName filesep ...
                            seqList(currSeq).name filesep '*.' p.nameSuffix]);
                        
                        % Break, if no new files
                        if length(fileList1)==length(fileList2)
                            clear fileList1 fileList2
                            break
                        end
                    end
                    
                    % Exit wait loop
                    continue_waiting = false;
                    
                    % Otherwise, wait
                else
                    
                    % Initialize wait bar
                    hW = waitbar(0,'1','Name','Analysis paused',...
                        'CreateCancelBtn',...
                        'setappdata(gcbf,''canceling'',1)');
                    setappdata(hW,'canceling',0)
                    
                    % Wait loop
                    for i = 1:ceil(waitDur*60/loopDur)
                        
                        % Check for Cancel button press
                        if getappdata(hW,'canceling')
                            
                            % Kill the waitbar
                            close(hW,'force')
                            
                            % Stop all code
                            return
                            
                            % Otherwise, update waitbar status
                        else
                            waitbar(i/ceil(waitDur*60/loopDur),hW,...
                                ['Waiting: ' sprintf('%6.1f',i*loopDur/60) ' of ' ...
                                num2str(waitDur) ' min'])
                        end
                        
                        % Pause to update status
                        pause(loopDur)
                    end
                    
                    % Kill waitbar for next pass thru loop
                    close(hW,'force')
                end
            end
            clear continue_waiting
        end
        
    end
    
end


%% Interactively visualize results for one sequence

if visData   
    
    
    playData(paths,p,batchName,cal,roi) 
             
end


 











