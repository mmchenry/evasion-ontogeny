function playData(paths,p,batchName,cal,roi) 
% Animates the results of the acqMaster data analysis

%% Parameters

% Declare global variables
global playOn cFrame playRate playDir gH

% fps
playRate = 30;

% Play direction 
playDir = 1;
   

%% Load data

% List of sequences
seqList = dir([paths.data filesep batchName filesep 'C1S*']);

% Check directory
if isempty(seqList)
    error(['No sequences in: ' paths.data filesep batchName]);
end

% If more than one sequence, prompt for selection
if length(seqList)>1
    % Format sequence list
    for i = 1:length(seqList)
        list_txt{i} = seqList(i).name;
    end
    
    % Prompt for sequence selection
    [currSeq,ok] = listdlg('PromptString','Choose sequence:',...
        'ListString',list_txt,'SelectionMode','single');
    
    % Quit if canceled
    if ok==0
        return
    end
    
    % Otherwise, skip the interaction
else
    currSeq = 1;
end

% Get Paths
[vPath,tPath,dPath,cPath] = setPaths(paths,batchName,seqList(currSeq).name);

% Load blob data (B)
load([dPath filesep 'blob data.mat'])

% Load sequence data (bStats)
load([dPath filesep 'Sequence stats.mat'])

% Load coordinate data (pd)
%load([dPath filesep 'fish_coords'])

% Load filenames for full frames
a = dir([vPath filesep '*.' p.nameSuffix]);

% Load list of thumbnails
%aT = dir([tPath filesep '*.' p.nameSuffix]);

% Check frames against data
if isempty(a)
    error(['No video frames in ' vPath]);
    
elseif length(a)~=length(B)
    error('Number of images does not match data');
end


%% Give commands
disp('Commands:')
disp('  spacebar      - Play/stop video')
disp('  right arrow   - Advance frame')
disp('  left arrow    - Back one frame')
disp('  up arrow      - Increase framerate')
disp('  down arrow    - Decrease framerate')
disp('  r             - Reverse play direction')


%% Create figure

% Initialize current frame number
cFrame = 1;

% Whether to play movie
playOn = 0;

% Calculate time series data
%D = calcTimeseries(B,p);

% Make figure window
f = figure('DoubleBuffer','on',...
           'WindowKeyPressFcn', {@but_down,bStats,B,cal,roi,p,tPath,vPath,a},...
           'WindowButtonMotionFcn', {@mouseMove},...
           'WindowButtonDownFcn', {@mouseDown,bStats,B,cal,roi,p,tPath,vPath,a});

gH = [];

update_fig(bStats,B,cal,roi,p,tPath,vPath,a)

        
  
function D = calcTimeseries(B,p)
% Make time vector
D.t = [1:length(B)]' ./ p.framerate;

% Loop thru time
for i = 1:length(B)
    if ~isnan(B(i).xMid(1))
        D.angl(i,1) = atan2(B(i).yMid(1)-B(i).yMid(2),...
                            B(i).xMid(1)-B(i).xMid(2));
        D.kappa(i,1) = hypot(B(i).yMid(1)-B(i).yMid(end),...
                             B(i).xMid(1)-B(i).xMid(end));                
    else
        D.angl(i,1)     = nan;
        D.kappa(i,1)    = nan;
    end
    
end
    
% Unwrap data
D.angl = unwrap(D.angl);
  

function play_vid(bStats,B,cal,roi,p,tPath,vPath,a)

% Declare global variables
global playOn cFrame playRate playDir gH

% Start timer
tStart = tic;

% Update figure window
update_fig(bStats,B,cal,roi,p,tPath,vPath,a)

% Duration for rendering
tRender = toc(tStart);

while playOn

        % If at the end of frames . . .
        if cFrame==length(B)
            playOn = 0;
            
        % Otherwise, advance
        else
            % If play period is longer than time to render 
            if (1/playRate) > tRender
                
                % Pause for necessary duration
                pause(1/playRate-tRender);
                
                % Advance to next
                cFrame = cFrame + playDir;
                
            else
                % Advance by necessary number of frames
                cFrame = cFrame + playDir*round(playRate * tRender);
            end        
        end
        
        % Update figure window
        update_fig(bStats,B,cal,roi,p,tPath,vPath,a)
     
end


function but_down(fig, key, bStats, B, cal, roi, p, tPath, vPath, a)
% Actions when a key is pressed

% Declare global variables
global playOn cFrame playRate playDir gH

% Set figure to current
figure(fig);


% NUMBER ----------------------------------------------------------
if sum(key.Key(1)==num2str([1:9]))==1
    % Stop playing
    playOn = 0;
    
    % Set relative frame number
    cFrame = round((str2num(key.Key)-1)/9*length(B));
    cFrame = max([cFrame 1]);
    cFrame = min([cFrame length(B)]);
    
% SPACEBAR --------------------------------------------------------
elseif strcmp(key.Key,'space')
    
    % Toggle play
    if playOn==1
        playOn = 0;
    else
        playOn = 1;
         play_vid(bStats,B,cal,roi,p,tPath,vPath,a);
    end
    
% LEFT ARROW ------------------------------------------------------
elseif strcmp(key.Key,'leftarrow')
    % Subtract by one
    cFrame = max([1 cFrame-1]);
    
% RIGHT ARROW -----------------------------------------------------
elseif strcmp(key.Key,'rightarrow')
    %Add by one
    cFrame = min([length(B) cFrame+1]);
    
% UPARROW ---------------------------------------------------------
elseif strcmp(key.Key,'uparrow')
    if playRate>90
        playRate = playRate + 60;
    elseif (playRate>15) 
        playRate = playRate + 15;
    else
        playRate = playRate + 5;
    end
    
% DOWNARROW -------------------------------------------------------
elseif strcmp(key.Key,'downarrow')
    if playRate >150 
        playRate = playRate - 60;
    elseif playRate>30
        playRate = playRate - 15;
    else
        playRate = max([1 playRate-5]);
    end

elseif strcmp(key.Key,'r') || strcmp(key.Key,'R')
    playDir = -1*playDir;
    
end

update_fig(bStats,B,cal,roi,p,tPath,vPath,a);


function update_fig(bStats,B,cal,roi,p,tPath,vPath,a)

global cFrame playOn playRate playDir gH

% Current time (s)
cTime = cFrame/p.framerate;

% Read full frame
[im,cmap] = imread([vPath filesep a(cFrame).name]);

% Read blob image (imBlob)
%load([tPath filesep a(cFrame).name]);
fName = [a(cFrame).name(1:(end-4)) '.mat'];

% Load 'blob' data, if present
if ~isempty(dir([tPath filesep fName]))
    load([tPath filesep fName]);
else
    blob.xMid = nan;
end

% Undistort
im = undistortImage(imadjust(im), cal.cameraParams,'OutputView','full');

% Blob region of interest
if ~isnan(blob.xMid(1))
    roiB.x = [blob.roi_blob(1) ...
        blob.roi_blob(1)+blob.roi_blob(3) ...
        blob.roi_blob(1)+blob.roi_blob(3) ...
        blob.roi_blob(1) ...
        blob.roi_blob(1)];
    roiB.y = [blob.roi_blob(2) ...
        blob.roi_blob(2) ...
        blob.roi_blob(2)+blob.roi_blob(4) ...
        blob.roi_blob(2)+blob.roi_blob(4) ...
        blob.roi_blob(2)];
else
    roiB.x = nan; roiB.y = nan;
end


% Show frame with overlaid data
gH(1) = subplot(5,1,1:2);

warning off
hFull = imshow(imadjust(im),cmap);
warning on
hold on
plot([roi.x;roi.x(1)],[roi.y;roi.y(1)],'y-')
plot(roiB.x,roiB.y,'r-')

if ~isnan(blob.xMid(1))
    h = plot(blob.xMid+roiB.x(1),blob.yMid+roiB.y(1),'g-',...
             blob.xEye+roiB.x(1),blob.yEye+roiB.y(1),'yo');
    
%     h = plot([blob.origin(1) blob.rost(1)]+roiB.x(1),...
%         [blob.origin(2) blob.rost(2)]+roiB.y(1),'y-', ...
%         blob.rost(1)+roiB.x(1),blob.rost(2)+roiB.y(1),'yo');
    set(h(2),'MarkerFaceColor','y')
    set(h(2),'MarkerSize',2)
end

hold off
if playDir==1
    title(['     Frame ' num2str(cFrame) ' (' num2str(playRate) 'fps) ->'])
else
    title(['<- Frame ' num2str(cFrame) ' (' num2str(playRate) 'fps)'])
end

gH(2) = subplot(5,1,3);
if ~isnan(blob.xMid(1))
    warning off
    imshow(imadjust(blob.im))
    warning on
    hold on
    
    % Plot border
    plot([1 size(blob.im,2) size(blob.im,2) 1 1],...
        [1 1 size(blob.im,1) size(blob.im,1) 1],'r-')
    
else
   imshow(uint8([255]),[],'InitialMagnification','fit')
end

%plot(yPerim,xPerim,'r-')
if ~isnan(blob.xMid(1))
    h = plot(blob.xMid,blob.yMid,'g-',...
             blob.xEye,blob.yEye,'yo');
    
%    h = plot(blob.xMid,blob.yMid,'go')
%     h = plot([blob.origin(1) blob.rost(1)],...
%         [blob.origin(2) blob.rost(2)],'y-', ...
%         blob.rost(1),blob.rost(2),'yo');


    set(h(2),'MarkerFaceColor','y')
    set(h(1),'LineWidth',2)
    set(h(2),'MarkerSize',7)
%     h = plot(blob.eye.R(1),blob.eye.R(2),'yo',...
%         blob.eye.L(1),blob.eye.L(2),'yo');
%     set(h,'MarkerFaceColor','y')
end
hold off

gH(3) = subplot(5,1,4);
plot(bStats.t,bStats.bSpan,'k.')
hold on
plot([1 1].*cTime,ylim,'r-')
hold off
xlabel('Time (s)')
ylabel('Tail curvature (rad/pix)')

gH(4) = subplot(5,1,5);
plot(bStats.t,bStats.angl,'k.')
hold on
plot([1 1].*cTime,ylim,'r-')
hold off
xlabel('Time (s)')
ylabel('Head angle (rad)')

% Pause to render
pause(1e-2);


%% Mouse position Callback
function [X,Y] = mouseMove(fig,~)
% Action when the cursor is moved over the figure window

global cFrame playOn playRate playDir gH

% Get subplot data
xLim3 = get(gH(3),'XLim');
yLim3 = get(gH(3),'YLim');

% Collect current point
C3 = get(gH(3), 'CurrentPoint');

% Get subplot data
xLim4 = get(gH(4),'XLim');
yLim4 = get(gH(4),'YLim');

% Collect current point
C4 = get(gH(4), 'CurrentPoint');

% If cursor in first  or second timeseries . . .
if ((C3(1,1)>=xLim3(1)) && (C3(1,1)<=xLim3(2)) && ...
    (C3(1,2)>=yLim3(1)) && (C3(1,2)<=yLim3(2))) || ...
   ((C4(1,1)>=xLim4(1)) && (C4(1,1)<=xLim4(2)) && ...
    (C4(1,2)>=yLim4(1)) && (C4(1,2)<=yLim4(2)))
   
    % Set cursor
     set(gcf,'Pointer','crosshair')
     

% Otherwise . . .
else
    % Set cursor
     set(gcf,'Pointer','arrow')
end


function [X,Y] = mouseDown(fig,~,bStats,B,cal,roi,p,tPath,vPath,a)
% Action when the cursor is moved over the figure window

global cFrame playOn playRate playDir gH

% Get subplot data
xLim3 = get(gH(3),'XLim');
yLim3 = get(gH(3),'YLim');

% Collect current point
C3 = get(gH(3), 'CurrentPoint');

% Get subplot data
xLim4 = get(gH(4),'XLim');
yLim4 = get(gH(4),'YLim');

% Collect current point
C4 = get(gH(4), 'CurrentPoint');

% If cursor in first  or second timeseries . . .
if ((C3(1,1)>=xLim3(1)) && (C3(1,1)<=xLim3(2)) && ...
    (C3(1,2)>=yLim3(1)) && (C3(1,2)<=yLim3(2))) || ...
   ((C4(1,1)>=xLim4(1)) && (C4(1,1)<=xLim4(2)) && ...
    (C4(1,2)>=yLim4(1)) && (C4(1,2)<=yLim4(2)))
   
    % Normalize current time
    rel_time = C3(1,1)./bStats.t(end);
    
    % Update current frame number
    cFrame = round(rel_time.*B(end).fr_num);
end

% Update figure
update_fig(bStats,B,cal,roi,p,tPath,vPath,a);

