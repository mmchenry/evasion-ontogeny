function anaData(batchName,expName)

% anaData loads raw midline, prey, and eye data and runs a simple analysis
% that finds where turns occur and extracts data from these portions. 

if nargin < 2
    batchName   = '2016-04-14';
    expName     = 'S04';
end

% Indicator for visualizing spline fits
vis = 0;

close all;

%% Path definitions

% % Matt's computer
if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')]))
    
    % Directory root
    %root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    root     = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';
    vid_root = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';

    % % Alberto's MacMini (office)
elseif ~isempty(dir([filesep fullfile('Users','alberto','Documents','GitHub-SourceTree')]))
    % Directory root
    root     = '/Volumes/VisualPred/ZF_visuomotor';
    vid_root = '/Volumes/VisualPred/ZF_visuomotor';
        
else
    % % Alberto's MacBook (laptop)
    path = fullfile('Users','A_Soto','Documents','Home');
    if ~isempty(dir([filesep path]))
        
        % Directory root
        root     = '/Volumes/BackUp/ZF_visuomotor';
        vid_root = '/Volumes/BackUp/ZF_visuomotor';
    
    else
        error('This computer is not recognized')
    end
    
end

% To raw video files 
paths.rawvid = [vid_root filesep 'Raw video'];

% To data files
paths.data = [root filesep 'Data'];

% To video thumbnails (for post-processing)
paths.thumb  = [root filesep 'Thumbnail video'];

% Directory for current data
dPath = [paths.data filesep batchName filesep expName];

%% Load data

% Load predator midline data
load([dPath filesep 'midline data.mat'])

% Load eye & heading data
load([dPath filesep 'eye data.mat'])

% Load prey data
load([dPath filesep 'prey data.mat'])

%% Take care of NaNs

% indices of digitized frames
indxPrey = find(isfinite(prey.xPrey));
indxPred = find(isfinite(mid.xRost));

% limited by Prey data; choose pred data to coincide with Prey

% indices for first digitized frame in common
ind1 = max(indxPrey(1),indxPred(1));
% ind1 = 90;

% index for last Prey frame
ind2 = max(indxPrey);
% ind2 = 865;

%% Load/Create Pred turning data parameters

% Check for (pred) turn data parameters, if they exist, load
if ~isempty(dir([dPath filesep 'turn data.mat']))
    
    load([dPath filesep 'turn data.mat'])
    
    % Check if start & end frames have been set
    if isfield(sp,'frStart')
        
        % Set start & end frames
        frStart = max(sp.frStart,ind1);
        frEnd   = min(sp.frEnd,ind2);
        
        % if not...set them
    else
        frStart = max(1,ind1);
        frEnd   = min(length(eyes.t),ind2);
        
        sp.frStart = frStart;
        sp.frEnd = frEnd;
    end
    
    % Put all tolerances into subfield 'tol'
    if isfield(sp,'tolHead')
        
%         % save tolerance values into 'tol' field
%         sp.tol.Head = sp.tolHead;
%         sp.tol.Reye = sp.tolReye;
%         sp.tol.Leye = sp.tolLeye;
        
        sp = rmfield(sp,{'tolHead','tolReye','tolLeye'});
    else
    end
    
    % Get field names in structure 'eyes'
    fNames = fieldnames(eyes);
    
    % Extract eye/heading data range specified by start & end frames
    for j=1:numel(fNames)
        if ~strcmp(fNames{j},'tol')
            eyes.(fNames{j}) = eyes.(fNames{j})(frStart:frEnd,:);
        else
        end
    end
    
    % otherwise...create it
else
%     
    % Set threshold value for peak angular velocity
    sp.thresh = 5;              % (rad/s)

    % Set start & end frames
    frStart = max(1,ind1);
    frEnd   = min(length(eyes.t),ind2);
%     % frEnd = min(ind2,1208);
   
    % Save start & end frames in 'sp' structure
    sp.frStart = frStart;
    sp.frEnd = frEnd;
    
    % Get field names in structure 'eyes'
    fNames = fieldnames(eyes);
    
    % Extract eye/heading data range specified by start & end frames
    for j=1:numel(fNames)
        if ~strcmp(fNames{j},'tol')
            eyes.(fNames{j}) = eyes.(fNames{j})(frStart:frEnd,:);
        else
        end
    end
    
end

% Unwrap heading angle
eyes.hdAngle = unwrap(eyes.hdAngle);

% Compute gaze angle (world coordinates)
gazeR = ((eyes.hdAngle + eyes.rAngle) - pi/2);
gazeL = ((eyes.hdAngle + eyes.lAngle) + pi/2);

%% Extract relevant Prey and Pred position data

% Prey field names
fNames = fieldnames(prey);

% Extract prey data range specified by start & end frames
for j=1:numel(fNames)
    prey.(fNames{j}) = prey.(fNames{j})(frStart:frEnd,:);
end

% Mid field names
fNames = fieldnames(mid);

% Extract pred data range specified by start & end frames
for j=1:numel(fNames)
    mid.(fNames{j}) = mid.(fNames{j})(frStart:frEnd,:);
end

%% Set smoothing tolerance parameters with smoothGUI 

% Time vector
sp.time    = eyes.t;

% Call interactive smoothing GUI to obtain smoothing tolerances
if ~isfield(sp,'tolSet') || 0
    
    % Clear out old 'tol' subfield
    sp = rmfield(sp,'tol');
    
    % Set up initial spline parameters
    sp.tol.Head = 0.0005;
    sp.tol.Reye = 0.0001;
    sp.tol.Leye = 0.0001;
    sp.tol.Prey = 10;
    sp.tol.Pred = 1.5;
    sp.tol.Eyes = 5;
    
    % Get tolerance parameter names
    varNames = fieldnames(sp.tol);
    
    % Data structure for splines
    data.Head  = eyes.hdAngle;
    data.Reye  = eyes.rAngle;
    data.Leye  = eyes.lAngle;
    data.Prey  = prey.xPrey;
    data.Pred  = mid.xRost;
    data.xReye = eyes.xReye;
    
    sp = smoothGUI(sp,varNames,data);
else
    % They've already been set...continue
end


%% Pred Spline fits and derivatives

% Spline fit heading angle
sp.hdAngle  = spaps(eyes.t,eyes.hdAngle,sp.tol.Head);

% Spline fit right eye angle
sp.rAngle   = spaps(eyes.t,eyes.rAngle,sp.tol.Reye);

% Spline fit left eye angle
sp.lAngle   = spaps(eyes.t,eyes.lAngle,sp.tol.Leye);

% Spline fit gaze angles
sp.gazeR    = spaps(eyes.t,gazeR,sp.tol.Head);
sp.gazeL    = spaps(eyes.t,gazeL,sp.tol.Head);

% Spline fit eye positions
sp.xReye    = spaps(eyes.t,eyes.xReye,sp.tol.Eyes);
sp.yReye    = spaps(eyes.t,eyes.yReye,sp.tol.Eyes);
sp.xLeye    = spaps(eyes.t,eyes.xLeye,sp.tol.Eyes);
sp.yLeye    = spaps(eyes.t,eyes.yLeye,sp.tol.Eyes);

% Check accurary of spline fits
if vis
    
    close all;
    
    % Generate smooth heading, and eye angle values
    hdAngle     = fnval(sp.hdAngle,eyes.t);
    rAngle      = fnval(sp.rAngle,eyes.t);
    lAngle      = fnval(sp.lAngle,eyes.t);
    
    figure,
    ax1 = subplot(2,1,1);
    % Plot heading angle
    plot(ax1,eyes.t,eyes.hdAngle * 180/pi,'.')
    hold on
    % Plot spline fit
    plot(ax1,eyes.t,hdAngle * 180/pi,'-k','LineWidth',2)
    hold off;
    
    ax2 = subplot(2,1,2);
    % Plot right eye angle
    plot(ax2,eyes.t,eyes.rAngle * 180/pi,'.')
    hold on
    % Plot left angle
    plot(ax2,eyes.t,eyes.lAngle * 180/pi,'.')
    % Plot spline fits
    plot(ax2,eyes.t,rAngle * 180/pi,'LineWidth',2)
    plot(ax2,eyes.t,lAngle * 180/pi,'LineWidth',2)
    hold off
    
    beep
    pause
end

% Compute angular velocity (rad/s) & accleration of heading
sp.hdAngle_D1    = fnder(sp.hdAngle,1);
% hdAngle_D       = fnval(sp.hdAngle_D,eyes.t);

sp.hdAngle_D2   = fnder(sp.hdAngle_D1);

% Compute angular velocity (rad/s) & accleration of right eye angle
sp.rAngle_D1     = fnder(sp.rAngle,1);
% rAngle_D        = fnval(sp.rAngle_D,eyes.t);

sp.rAngle_D2    = fnder(sp.rAngle_D1);

% Compute angular velocity (rad/s) & accleration of left eye angle
sp.lAngle_D1     = fnder(sp.lAngle,1);
% lAngle_D        = fnval(sp.lAngle_D,eyes.t);

sp.lAngle_D2    = fnder(sp.lAngle_D1);

if vis
    fnplt(sp.hdAngle_D1)
    pause
end


%% Prey & Pred Position Data Splines 

% Spline fit prey position
sp.xPrey  = spaps(prey.t,prey.xPrey,sp.tol.Prey);  
sp.yPrey  = spaps(prey.t,prey.yPrey,sp.tol.Prey);  

% Spline fit pred rostrum position
sp.xRost  = spaps(mid.t,mid.xRost,sp.tol.Pred);  
sp.yRost  = spaps(mid.t,mid.yRost,sp.tol.Pred);

% xRost       = fnval(sp.xRost,mid.t);

% Save spline and turning data up to this point
save([dPath filesep 'turn data.mat'],'sp', 'D')

%% Bearing Angle

% Values for x-coordinate
xPrey_val = fnval(sp.xPrey,prey.t);
xPred_val = fnval(sp.xRost,mid.t);

% Values for y-coordinate (transform origin to bottom-left)
yPrey_val = 1024 - fnval(sp.yPrey,prey.t);
yPred_val = 1024 - fnval(sp.yRost,mid.t);

% x-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
x_R = xPrey_val - xPred_val; 

% y-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
y_R = yPrey_val - yPred_val;

% Direction of range/baseline vector
alpha_R = atan2(y_R, x_R);

% Magnitude of range/baseline vector (distance between prey & pred)
rangeMag = sqrt(sum([x_R, y_R].^2,2));

% Bearing Angle
phi = alpha_R - fnval(sp.hdAngle,eyes.t);

% Check angles
if 0
    subplot(3,1,1), plot(eyes.t,phi*180/pi)
    ylabel('Bearing Angle (deg)')
    
    subplot(3,1,2), fnplt(fncmb(sp.hdAngle,180/pi))
    ylabel('Predator Heading Angle (deg)')
    
    subplot(3,1,3), plot(eyes.t,alpha_R*180/pi)
    ylabel('Direction of Range vector (deg)')
    xlabel('time (s)')
end

% Spline fit bearing angles (with heading tolerance)
sp.bearAngl = spaps(prey.t,phi,sp.tol.Head);

% Spline fit range vector direction (alpha)
sp.alpha    = spaps(prey.t,alpha_R,sp.tol.Head);

%% Prelim analysis (Find peaks in angular velocity & time intervals)

% If intervals have not been set...
if ~isfield(D,'intSet')
    
    % Clear old data structure D
    clear D
    
    % Roots of first derivative (angular velocity)
    hd_D1roots = fnzeros(sp.hdAngle_D1);
    hd_D1roots = hd_D1roots(:);
    
    % Time points for local min and max of heading angle
    tHd = unique(hd_D1roots);
    
    % Roots of second derivative (angular acceleration)
    hd_D2roots = fnzeros(sp.hdAngle_D2);
    hd_D2roots = hd_D2roots(:);
    
    % Get peaks and troughs of angular velocity
    peaks = abs(fnval(sp.hdAngle_D1,hd_D2roots));
    
    % Indices of peaks (and troughs) in angular velocity
    idx1 = peaks > sp.thresh & peaks < sp.thresh*15;
    
    % Time points of peak angular velocity
    tAV = unique(hd_D2roots(idx1));
    
    % Values of peak angular velocity
    hdAngle_max = fnval(sp.hdAngle_D1,tAV);
    
    % Store max angular velocity
    D.hdAngle_maxAV = hdAngle_max;
    
    % Total number of turns
    numTurns    = length(tAV);
    D.numTurns  = numTurns;
    
    % For each peak in angular velocity, find the time interval for the turn
    for k=1:numTurns
        
        % Time points of zeros before peak angular velocity
        tAV_before = hd_D2roots(hd_D2roots < tAV(k));
        
        % Time points of zeros after peak angular velocity
        tAV_after = hd_D2roots(hd_D2roots > tAV(k));
        
        % Time points of local max/min heading angle before peak AV
        tHD_before = hd_D1roots(hd_D1roots < tAV(k));
        
        % Time points of local max/min heading angle after peak AV
        tHD_after = hd_D1roots(hd_D1roots > tAV(k));
        
        % Left endpoint of time interval of turn
        D.tInt(k,1) = max([tAV_before; tHD_before; prey.t(1)]);
        
        % Right endpoint of time interval of turn
        if isempty(tAV_after)
            D.tInt(k,2) = eyes.t(end);
        else
            D.tInt(k,2) = min([tAV_after; tHD_after]);
        end
        
        clear tAV_before tAV_after tHD_before tHD_after
    end
end

%% Interactive GUI for time intervals corrections

% Call interactive time intervals GUI
if ~isfield(D,'intSet') 
    D = intervalsGUI(D,numTurns,sp.hdAngle);
    
else
    % intervals have already been adjusted...continue
end

%% Prior time intervals

% Preallocate array for time intervals before turns (prior interval)
D.priorInt = zeros(size(D.tInt));

% Left endpoint of first prior interval (initial time)
D.priorInt(1,1) = mid.t(frStart);

% Left endpoints of prior interval 
D.priorInt(2:numTurns,1) = D.tInt(1:numTurns-1,2);

% Right endpoints of prior interval
D.priorInt(:,2) = D.tInt(:,1);

%% Compute changes in angles & absolute values before turn

% Change in orientation during turn
hdDelta         = diff(fnval(sp.hdAngle,D.tInt),1,2);
D.hdDelta       = hdDelta;

% Change in gaze angle in prior interval
gazeDelta       = diff(fnval(sp.gazeR,D.priorInt),1,2);
D.gazeDelta     = gazeDelta;

% Change in bearing angle in prior interval
bearDelta       = diff(fnval(sp.bearAngl,D.priorInt),1,2);
D.bearDelta     = bearDelta;

% Change in range vector direction (alpha) in prior interval
alphaDelta      = diff(fnval(sp.alpha,D.priorInt),1,2);
D.alphaDelta    = alphaDelta;

% Compute bearing angle at start of turn
bearingPre      = fnval(sp.bearAngl,D.tInt(:,1));
D.bearingPre    = bearingPre;

% Compute angular velocity (derivative of bearing) at start of turn
bearD1Pre       = fnval((fnder(sp.bearAngl)),D.tInt(:,1));
D.bearD1Pre     = bearD1Pre;

% matrix of all change in angle data
allData = [hdDelta gazeDelta bearDelta alphaDelta bearingPre bearD1Pre];

%% Save Data

% fileID = fopen([paths.data filesep 'DeltaAngle-Data.txt'],'a+');
% 
% % Number format
% fmt = '%.4f\t %.4f\t %.4f\t %.4f\n';
% 
% % print headers
% fprintf(fileID,'%s\t %s\t %s\t %s\n %s\n',...
%     'hdDelta', 'gazeDelta', 'bearDelta', 'alphaDelta',...
%     [batchName, '-', expName]);
% 
% % Go to end of file
% fseek(fileID,0,'eof')
% 
% % print data
% fprintf(fileID,fmt, allData);
% fclose(fileID);

% Save spline and turning data
save([dPath filesep 'turn data.mat'],'sp', 'D')


%% Visualize the results

if vis
    
    totFrames   = length(gazeR);
    frames      = eyes.t(1:totFrames)*250;
    
    figure
    subplot(2,1,1);
    plot(frames,unwrap(eyes.hdAngle)./pi*180,'-k');
    hold on;
    plot(frames,gazeR./pi*180,'-',frames,gazeL./pi*180,'-');
    grid on;
    legend('Heading','R Gaze','L Gaze');
    % xlabel('Time (s)')
    xlabel('Frame number')
    ylabel('Head/Gaze angle (deg)')
    
    subplot(2,1,2);
    plot(frames,eyes.rAngle./pi*180,'-')
    hold on
    plot(frames,eyes.lAngle./pi*180,'-')
    grid on
    legend('R','L');
    % xlabel('Time (s)')
    xlabel('Frame number')
    ylabel('Eye angle (deg)')

end

close all;

disp([batchName, '; exp ', expName])

