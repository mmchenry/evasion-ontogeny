function anaData(batchName,expName)

% anaData loads raw midline, prey, and eye data and runs a simple analysis
% that finds where turns occur and extracts data from these portions. 

if nargin < 2
    batchName   = '2016-03-14';
    expName     = 'S01';
end

% Indicator for visualizing spline fits
vis = 1;

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

%% Load/Create Pred turning data

% check for (pred) turn data, if it exists, load it
if ~isempty(dir([dPath filesep 'turn data.mat']))
    
    load([dPath filesep 'turn data.mat'])
    
    % Set start & end frames
    frStart = sp.frStart;
    frEnd = sp.frEnd;
    
    % otherwise...create it
else
    %% Pred Spline fits and derivatives
    
    % Threshold value for peak angular velocity
    sp.thresh = 5;              % (rad/s)
    
    % Set up spline parameters
    sp.tol.Head   = 0.0005;
    sp.tol.Reye   = 0.0009;
    sp.tol.Leye   = 0.0009;
    
    % Set start & end frames
    frStart = 1;
    frEnd = length(eyes.t);
    % frEnd = 1208;
    
    % Save start & end frames in 'sp' structure
    sp.frStart = frStart;
    sp.frEnd = frEnd;
    
    % Get field names in structure 'eyes'
    fNames = fieldnames(eyes);
    
    % Extract data range specified by start & end frames
    for j=1:numel(fNames)
        eyes.(fNames{j}) = eyes.(fNames{j})(frStart:frEnd,:);
    end
    
    % Unwrap heading angle
    eyes.hdAngle = unwrap(eyes.hdAngle);
    
    % Compute gaze angle (world coordinates)
    gazeR = ((eyes.hdAngle + eyes.rAngle) - pi/2);
    gazeL = ((eyes.hdAngle + eyes.lAngle) + pi/2);
    
    % Spline fit heading angle
    sp.hdAngle  = spaps(eyes.t,eyes.hdAngle,sp.tol.Head);
    
    % Spline fit right eye angle
    sp.rAngle   = spaps(eyes.t,eyes.rAngle,sp.tol.Reye);
    
    % Spline fit left eye angle
    sp.lAngle   = spaps(eyes.t,eyes.lAngle,sp.tol.Leye);
    
    % Spline fit gaze angles
    sp.gazeR    = spaps(eyes.t,gazeR,sp.tol.Head);
    sp.gazeL    = spaps(eyes.t,gazeL,sp.tol.Head);
    
    
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
    
    fnplt(sp.hdAngle_D1)
    pause
    
end

%% Prey & Pred Position Data Splines

sp.tol.Prey = 1;
sp.tol.Pred = 1.5;

% Prey field names
fNames = fieldnames(prey);

% Extract data range specified by start & end frames
for j=1:numel(fNames)
    prey.(fNames{j}) = prey.(fNames{j})(frStart:frEnd,:);
end

% Mid field names
fNames = fieldnames(mid);

% Extract data range specified by start & end frames
for j=1:numel(fNames)
    mid.(fNames{j}) = mid.(fNames{j})(frStart:frEnd,:);
end

% Spline fit prey position
sp.xPrey  = spaps(prey.t,prey.xPrey,sp.tol.Prey);  
sp.yPrey  = spaps(prey.t,prey.yPrey,sp.tol.Prey);  

% Spline fit pred rostrum position
sp.xRost  = spaps(mid.t,mid.xRost,sp.tol.Pred);  
sp.yRost  = spaps(mid.t,mid.yRost,sp.tol.Pred);

% xRost       = fnval(sp.xRost,mid.t);

%% Bearing Angle

% x-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
x_R = fnval(sp.xPrey,prey.t) - fnval(sp.xRost,mid.t); 

% y-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
y_R = 1024 - (fnval(sp.yPrey,prey.t) - fnval(sp.yRost,mid.t));

% Direction of range/baseline vector
theta_R = atan2(y_R, x_R);

% Magnitude of range/baseline vector (distance between prey & pred)
rangeMag = sqrt(sum([x_R, y_R].^2,2));

% Bearing Angle
phi = theta_R - fnval(sp.hdAngle,eyes.t);

%% Prelim analysis

% Roots of first derivative (angular velocity)
hd_D1roots = fnzeros(sp.hdAngle_D1);
hd_D1roots = hd_D1roots(:);

% Time points of local min and max heading angle
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

% Value of peak angular velocity
hdAngle_max = fnval(sp.hdAngle_D1,tAV);

% Store max angular velocity
D.hdAngle_max = hdAngle_max;

% For each peak, find the time interval for the turn

for k=1:length(tAV)
    
    % Time points of zeros before peak angular velocity
    tAV_before = hd_D2roots(hd_D2roots < tAV(k));
    
    % Time points of zeros after peak angular velocity
    tAV_after = hd_D2roots(hd_D2roots > tAV(k));
    
    % Time points of local max/min heading angle before peak AV
    tHD_before = hd_D1roots(hd_D1roots < tAV(k));
    
    % Time points of local max/min heading angle after peak AV
    tHD_after = hd_D1roots(hd_D1roots > tAV(k));
    
    % Time interval of turn
    D.tInt(k,1) = max([tAV_before; tHD_before]);
    
    if isempty(tAV_after)
        D.tInt(k,2) = eyes.t(end);
    else
        D.tInt(k,2) = min([tAV_after; tHD_after]);
    end
    
    % Change in orientation during turn
    D.hdDelta(k,1) = diff(fnval(sp.hdAngle,D.tInt(k,:)));
    
    clear tAV_before tAV_after tHD_before tHD_after
end

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


disp([batchName, '; exp ', expName])

