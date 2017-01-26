function D = anaData(mid,eyes,prey,batchName)
% anaData loads raw midline, prey, and eye data and runs a simple analysis
% that finds where turns occur and extracts data from these portions. 

% Default experiment
if nargin < 2
    batchName   = '2016-04-14';
    expName     = 'S01';
end
%Test


%% Parameter values

% Theshold deviation for spline fit (smaller values give tighter fit) 
D.p.dev_thresh = 0.1;

% Temporal density of roots that spline fit will not exceed
D.p.rootThresh = 40;

% High threshold rate of rotation for a turn (rad/s)
D.p.yawRate_high = 5;

% Low threshold rate of rotation for a turn (rad/s)
D.p.yawRate_low = 2;

% Height of video frame (pixels)
D.p.vid_height = 1024;

% Minimum duration for a single beat
D.p.min_dur = 0.1;

% Minimum duration for a scoot in prey
D.p.min_durPy = 0.15;

% Minium change in direction for a large tailbeat in predators (rad)
D.p.yawDelta = 10*pi/180;

% High threshold change in speed for a prey (?/s)
D.p.highSpdPy = 0.75;

% Low threshold change in speed for a prey (?/s)
D.p.lowSpdPy = 0.5;

% Vergence angle (rad)
D.p.vergAng = 33*pi/180;

% Field of view (rad)
D.p.fov = 186*pi/180;

% Conversion factor to convert from pixels to cm.
if str2num(batchName(1:4))==2016
    switch batchName(6:7)
        case {'02','03'}
            D.p.cF = 0.0064;        % for experiments before 03-23-16
            
            if strcmp(batchName,'2016-03-23')
                D.p.cF = 0.005;        % for experiments on 03-23-16
            end
            
        case {'04','06'}
            D.p.cF = 0.0050;        % for experiments in April 2016
            
        case {'08'}
            D.p.cF = 0.0065;        % for experiments in August 2016
            
        case {'09','10'}
            D.p.cF = 0.00763;       % for experiments in Sept. 2016
            
            if strcmp(batchName,'2016-09-08')
                D.p.cF = 0.00625;       % for experiments on 09-08-16
            end
            
%     if str2num(batchName(6:7))>3
%         D.p.cF = 0.0050;        % for experiments after 03-23-16
%     elseif str2num(batchName(6:7))>7
%         D.p.cF = 0.0065;        % for experiments after 06-17-16
%     elseif str2num(batchName(6:7))>9
%         D.p.cF = 0.00763;       % for experiments after 09-12-16
%     elseif strcmp(batchName,'2016-09-08')
%         D.p.cF = 0.00625;       % for experiments on 09-08-16
%     else
%         D.p.cF = 0.0064;        % for experiments before 03-23-16
    end
else
    error('This calibration code assumes all experiments run in 2016');
end

% Indicator for visualizing spline fits
vis = 0;


%% Take care of NaNs & extract relevant Prey and Pred position data

% indices of digitized frames
indxPrey = find(isfinite(prey.xPrey));
indxPred = find(isfinite(mid.xRost));
indxEyes = find(isfinite(eyes.xReye));

% limited by Prey data; choose pred data to coincide with Prey

% indices for first digitized frame in common
ind1 = max([indxPrey(1) indxPred(1)]);% indxEyes(1)]);

% index for last frame in common
ind2 = min([indxPrey(end),indxPred(end)]);%,indxEyes(end)]);

% Store frame numbers
sp.frStart = max(1,ind1);
sp.frEnd   = min(length(eyes.t),ind2);
      
% Prey field names
fNames = fieldnames(prey);

% Extract prey data range specified by start & end frames
% Begin with 3rd variable, first two are scalars
for j=3:8
    prey.(fNames{j}) = prey.(fNames{j})(sp.frStart:sp.frEnd,:);
end

% Mid field names
fNames = fieldnames(mid);

% Extract pred data range specified by start & end frames
for j=1:numel(fNames)
    mid.(fNames{j}) = mid.(fNames{j})(sp.frStart:sp.frEnd,:);
end

% Get field names in structure 'eyes'
fNames = fieldnames(eyes);

% Extract eye/heading data range specified by start & end frames
for j=1:numel(fNames)
    if ~(strcmp(fNames{j},'tol') || strcmp(fNames{j},'eyeData'))
        eyes.(fNames{j}) = eyes.(fNames{j})(sp.frStart:sp.frEnd,:);
    end
end

% Check length of fields
if (length(mid.xRost) ~= length(eyes.xReye)) || ...
    (length(mid.xRost) ~= length(prey.xPrey))
    error('The "mid", "eyes" & "prey" vectors need to all be the same length');
end

% Check start time 
if (mid.t(1) ~= eyes.t(1)) || (mid.t(1) ~= prey.t(1))
    error('The "mid", "eyes" & "prey" vectors need to start at the same time');
end

% Check end time 
if (mid.t(end) ~= eyes.t(end)) || (mid.t(end) ~= prey.t(end))
    error('The "mid", "eyes" & "prey" vectors need to end at the same time');
end

% Time vector
D.t    = eyes.t;

% Calculate tail bending index
mid.bend = findBending(D,mid);

% Calculate midline arclength and store in 'D' (cm)
D.bodyL = findLength(mid) * D.p.cF;

clear indexPrey indxPred ind1 ind2 indxPrey indxPred indxEyes fNames j 


%% Set smoothing tolerance parameters

% Get tolerance parameter names
varNames = {'Head','Reye','Leye','prey_ang','Prey','Pred','Eyes','Bend'};

% Data structure for splines
data.Head  = unwrap(eyes.hdAngle);
data.Reye  = eyes.rAngle;
data.Leye  = eyes.lAngle;
data.prey_ang = unwrap(prey.thetaPrey/180*pi);

% Data structure continued (all length units: cm)
data.Prey  = prey.xPrey.*D.p.cF;
data.Pred  = mid.xRost.*D.p.cF;
data.xReye = eyes.xReye.*D.p.cF;
data.bend  = mid.bend;

% Non-interactive method for finding tolerance
try
    D.sp = findSplineTols(D,sp,varNames,data);
    
    % Turn OFF inidicator for ignoring eye data
    spEye = 0;
catch
    % Eye data not complete for some sequences, but may not be relevant 
    % This can happen for sequences where eye or prey codes don't work 
    
    % Display warning
    disp('Eye and/or Prey data may be incomplete...')
    
    % Turn ON inidicator for ignoring eye data
    spEye = 1;
    
    % Tolerance parameter names needed (ignoring eye data)
    varNames    = {'Head','Prey','Pred','Bend'};
    
    % Remove nonessential fields from data structure
    data = rmfield(data,{'Reye','Leye','prey_ang','xReye'});
    
    % Find tolerance parameters on subset of all variables
    D.sp        = findSplineTols(D,sp,varNames,data);
    
    % Set tolerance parameters for other variables (arbitrary)  
    D.sp.tol.prey_ang   = D.sp.tol.Prey;
    D.sp.tol.Eyes       = D.sp.tol.Pred; 
    D.sp.tol.Reye       = D.sp.tol.Pred; 
end

clear sp data varNames


%% Smooth data with splines 

% Turn off warnings due to spline fitting
warning off

% Spline fit prey position
D.sp.xPy    = spaps(D.t,prey.xPrey.*D.p.cF,D.sp.tol.Prey);  
D.sp.yPy    = spaps(D.t,(D.p.vid_height-prey.yPrey).*D.p.cF,D.sp.tol.Prey);

% Split prey heading angle if data is complete
if ~spEye
    D.sp.angPy  = spaps(D.t,unwrap(prey.thetaPrey/180*pi),D.sp.tol.prey_ang);
D.posPy     = [fnval(D.sp.xPy,D.t) fnval(D.sp.yPy,D.t) fnval(D.sp.angPy,D.t)];
else
    D.posPy = [fnval(D.sp.xPy,D.t) fnval(D.sp.yPy,D.t)];
end


% Spline fit predator rostrum position & orientation
D.sp.xPd    = spaps(D.t,mid.xRost.*D.p.cF,D.sp.tol.Pred);  
D.sp.yPd    = spaps(D.t,(D.p.vid_height-mid.yRost).*D.p.cF,D.sp.tol.Pred);
D.sp.angPd  = spaps(D.t,unwrap(eyes.hdAngle),D.sp.tol.Head);
D.posPd     = [fnval(D.sp.xPd,D.t) fnval(D.sp.yPd,D.t) fnval(D.sp.angPd,D.t)];

% Get tail position
[D.sp.bend,D.bend]  = spaps(D.t,mid.bend,D.sp.tol.Bend);
D.bend = D.bend';

% Spline fit eye positions
try
D.sp.xReye    = spaps(eyes.t,eyes.xReye.*D.p.cF,D.sp.tol.Eyes);
D.sp.yReye    = spaps(eyes.t,(D.p.vid_height-eyes.yReye).*D.p.cF,D.sp.tol.Eyes);
D.sp.xLeye    = spaps(eyes.t,eyes.xLeye.*D.p.cF,D.sp.tol.Eyes);
D.sp.yLeye    = spaps(eyes.t,(D.p.vid_height-eyes.yLeye).*D.p.cF,D.sp.tol.Eyes);

% Spline fit eye angles
D.sp.angR   = spaps(D.t,eyes.rAngle,D.sp.tol.Reye);
D.sp.angL   = spaps(D.t,eyes.lAngle,D.sp.tol.Reye);

% Store discrete eye values
D.posR = [fnval(D.sp.xReye,D.t) fnval(D.sp.yReye,D.t) fnval(D.sp.angR,D.t)];
D.posL = [fnval(D.sp.xLeye,D.t) fnval(D.sp.yLeye,D.t) fnval(D.sp.angL,D.t)];

catch 
% NOTE: Figure out what to do with eye/prey data that is incomplete
% Quick and dirty fix: make all subfields empty
    if spEye
        D.sp.xReye = [];
        D.sp.yReye = [];
        D.sp.xLeye = [];
        D.sp.yLeye = [];
        D.sp.angR  = [];
        D.sp.angL  = [];
        D.posR     = [];
        D.posL     = [];
    end
end

warning on

% Prey diameter if assumed a sphere: weighted average of long and short axes
D.lenPy = (0.30*mean(prey.MajorAxis)*D.p.cF + 0.70*mean(prey.MinorAxis))*D.p.cF;

% Check accurary of spline fits
if vis
    % Plot x-position
    subplot(5,1,1);
    plot(mid.t,mid.xRost*D.p.cF,'r.',D.t,D.posPd(:,1),'k-')
    grid on
    ylabel('xRost')
    title('Spline fits')
    
    % Plot heading angle (unwrapped heading)
    subplot(5,1,2);
    plot(mid.t,unwrap(eyes.hdAngle),'r.',D.t,D.posPd(:,3),'k-')
    grid on 
    ylabel('Head angle (rad)')

    % Plot eye angle
    subplot(5,1,3);
    plot(mid.t,eyes.rAngle,'r.',D.t,D.posR(:,3),'k-')
    grid on 
    ylabel('Right eye angle (rad)')
    
    % Plot bending
    subplot(5,1,4);
    plot(mid.t,mid.bend,'r.',D.t,D.bend,'k-')
    grid on 
    ylabel('Bending index')
    
    % Plot prey position
    subplot(5,1,5);
    plot(prey.t,prey.xPrey*D.p.cF,'r.',D.t,D.posPy(:,1),'k-')
    grid on 
    ylabel('Prey x-position')
   
    %beep
%     pause
end

% Clear rawdata variables
clear eyes mid prey posTail


%% Define the prey in the predator's FOR

if ~spEye
    
    for i = 1:length(D.t)
        
        % Local coord system, Both eyes
        D.tformR{i} = local_system_R(D.posR(i,1:2),D.posL(i,1:2));
        D.tformL{i} = local_system_L(D.posR(i,1:2),D.posL(i,1:2));
        
        % Position of prey in the eye FORs
        [xPyR,yPyR] = global_to_local(D.tformR{i},D.posPy(i,1),D.posPy(i,2));
        [xPyL,yPyL] = global_to_local(D.tformL{i},D.posPy(i,1),D.posPy(i,2));
        
        angR = atan2(yPyR,xPyR);
        angL = atan2(yPyL,xPyL);
        
        D.posPyR(i,:) =  [xPyR yPyR angR];
        D.posPyL(i,:) =  [xPyL yPyL angL];
        
        % Visual check on transformation
        if 0
            subplot(2,2,[1 3])
            plot([D.posR(i,1),D.posL(i,1)],[D.posR(i,2),D.posL(i,2)],'k-o',...
                D.posPy(i,1),D.posPy(i,2),'ro',D.posR(i,1),D.posR(i,2),'k*')
            axis equal;grid on
            title('Intertial FOR')
            
            subplot(2,2,2)
            plot(0,0,'k*',D.posPyR(i,1),D.posPyR(i,2),'ro')
            title('Right eye')
            axis equal;grid on
            
            subplot(2,2,4)
            plot(0,0,'ko',D.posPyL(i,1),D.posPyL(i,2),'ro')
            title('Left eye')
            axis equal; grid on
        end
    end
end



%% Prelim analysis (Find peaks in angular velocity & time intervals)

% Find periods of tail beating (predator)
D = find_pred_intervals(D);

% Find periods of tail beating (prey)
D = find_prey_intervals(D);

% Visualize results
if 0
    vis_beats(D,'prey')
    
    vis_beats(D,'pred')
end



function sp = findSplineTols(D,sp,varNames,data)
% Algrotihm to determine the smoothing tolerance for a spline.  It starts
% at a high value and stops when either the number of rootThresh is
% exceeded or the dev_thresh drops below some level.

% Get data fieldnames
dataNames = fieldnames(data);


% Loop thru each variable
for j = 1:length(varNames)
    
    % Indicator for setting tolerance (0 means tolerance is not yet set)
    tolset = 0;
    
    % Name of current variable
    currVar = varNames{j};
     
    % Initial smoothing tolerance for current variable
    %spTol = sp.tol.(currVar);
    spTol = 1e2;
    
    % Time vector
    t = D.t;
    
    % Current data
    d = data.(dataNames{j});
    
    % Start timer
    tstart = tic;
    
    % Check for vector of zeros
    if sum(d)==0
        spTol = 0;
        
    else
        while true
            
            % Reduce tolerance value
            spTol = spTol*0.9;
            
            % Current smoothing spline fit
            warning off
            spCurr  = spaps(t,d,spTol);
            warning on
            Dsp     = fnder(spCurr,1);
            D2sp    = fnder(spCurr,2);
            
            % Roots of second derivative (time of peak/troughs)
            D2roots = fnzeros(D2sp);
            D2roots = D2roots(1,:)';
            
            % Normalized deviation
            dev = sqrt(sum((d - fnval(spCurr,t)).^2))/range(d);
            
            % Root density
            root_den = length(D2roots)/range(t) ;
            
            % Elapsed time since start of loop
            telapsed = toc(tstart);
            
            % Check to quit
            if (root_den>D.p.rootThresh) || (dev<D.p.dev_thresh)
                break
                
            elseif telapsed > 120
                error('Time out: Spline optimization taking too long')
            end
        end
    end
    
    % Store results
    sp.tol.(currVar)    = spTol;
    
    % If plotting data . . .
    if 0
        figure
        % Plot raw data and overlay smoothed data
        subplot(2,1,1)
        
        % Plot raw data
        plot(t,d,'.');hold on; legend(num2str(spTol))
        
        % Plot smooth data
        fnplt(spCurr);hold off;grid on
        xlabel('Time (s)');ylabel(currVar)
        
        subplot(2,1,2)
        fnplt(Dsp);set(gca,'ColorOrderIndex',1);hold on
        
        % Plot peaks and troughs
        plot(D2roots,fnval(Dsp,D2roots),'k+');hold off;grid on        
        drawnow
    
        subplot(2,1,1)
        title(['tol = ' num2str(spTol) ', dev = ' num2str(dev) ...
            ', root_den = ' num2str(root_den)])
    end
    
    clear Dsp D2sp
end

%pause(1)

% Create the field 'tolSet', indicating that the tolerances have been set.
sp.tolSet = 1;


function bend = findBending(D,mid)
% Defines metric of bending

% If no midline data . . .
if ~isfield(mid,'sMid')
    bend = zeros(length(mid.t),1);
    
else
    % Initialize bend parameter vector
    bend = zeros(length(mid.t),1);
    
    % Loop thru midline points
    for i = 1:length(mid.sMid)
        
        % Check for empty midline data
        if isempty(mid.sMid{i})
            try
                bend(i,1) = mean(bend(i-4:i-1,1));
            catch
                bend(i,1) = 0;
            end
        else
            
            % Extract
            s = mid.sMid{i}./mid.sMid{i}(end);
            x = mid.xMid{i}./mid.sMid{i}(end);
            y = mid.yMid{i}./mid.sMid{i}(end);
            
            cX = polyfit(s,x,1);
            cY = polyfit(s,y,1);
            
            bend(i,1) = sum(sqrt((x-polyval(cX,s)).^2 + (y-polyval(cY,s)).^2));
        end
            if 0
                plot(x,y,'k-',polyval(cX,s),polyval(cY,s),'-r')
                title(['t = ' num2str(D.t(i)) ': bend = ' num2str(bend(i))])
                pause(0.3)
            end
    end
end

function bodyL = findLength(mid)
% Computes length of midline segments, to be used for computing fish length

% If no midline data . . .
if ~isfield(mid,'sMid')
    bodyL = zeros(length(mid.t),1);
    disp('   No midline data, cannot compute body length!')
    
else
    % Initialize body length vector
    bodyL = zeros(length(mid.t),1);
    
    % Loop thru midline points
    for i = 1:length(mid.sMid)
        
        % Check for empty midline data
        if isempty(mid.sMid{i})
            try
                bodyL(i,1) = mean(bodyL(i-4:i-1,1));
            catch
                bodyL(i,1) = 0;
            end
        else
            
            % Extract arc-length data
            bodyL(i,1) = mid.sMid{i}(end);
            
            % Sum up the arc length values for current frame
%             bodyL(i,1) = sum(s);
        end
%             if 0
%                 plot(x,y,'k-',polyval(cX,s),polyval(cY,s),'-r')
%                 title(['t = ' num2str(D.t(i)) ': bend = ' num2str(bend(i))])
%                 pause(0.3)
%             end
    end
end
    


function D = find_prey_intervals(D)
% Algorithm for determining the time intervals for tail beats and tail
% flicks.  These are saved in the 'tBeat' field of D in a 3xn matrix, with
% the first column denoting tail beats (1) and flicks (0).

% Current smoothing spline fit
spd = give(D,'prey spd');
sp      = spline(D.t,spd);
Dsp     = fnder(sp,1);

% Time of changes in sign of rate of change
tHd = fnzeros(Dsp);
tHd = unique(tHd(1,:)');

% Time of scoots above threshold spd (candidate)
tBig = tHd(fnval(sp,tHd)>D.p.highSpdPy);

% Check that tBig is nonempty
if isempty(tBig)
    disp('tBig is empty...')
end

% Time of change in spd below threshold
tSmall = tHd(fnval(sp,tHd)<D.p.lowSpdPy);

% Start of changes in angular direction before peak
tPrior = tSmall(tSmall<tBig(1));

% Ending of changes in angular direction after peak
tFollow = tSmall(tSmall>(tBig(1)+D.p.min_durPy));

% Set starting time for first tail beat
if isempty(tPrior)
    tBeat = D.t(1);
else
   tBeat = tPrior(end); 
end

% Set time for end of first tail beat
if isempty(tFollow)
   tBeat(1,2) = D.t(end);
else
   tBeat(1,2) = tFollow(1);
end

% Loop thru candidate times
for i = 2:length(tBig)
    
    % Starts of changes in prior to current iBig
    tPrior = tSmall((tSmall<tBig(i)) & (tSmall>=(tBig(i-1))));
    
    % Ending of changes in angular direction after current peak
    tFollow = tSmall(tSmall>(tBig(i)+D.p.min_durPy));
    
    % If next candidate is beyond the min duration . . .
    if ~isempty(tPrior) && ~isempty(tFollow)
        tBeat = [tBeat; tPrior(end) tFollow(1)];
    elseif ~isempty(tPrior) && isempty(tFollow)
        tBeat = [tBeat; tPrior(end) D.t(end)];
    end
end

% Code all scoots as '1'
D.tBeatPy   = [ones(size(tBeat,1),1) tBeat];

% Visualize
if 0
    subplot(3,1,[1:2]); 
    plot(D.t,spd,'-',tHd,fnval(sp,tHd),'+k',tBig,fnval(sp,tBig),'ro',...
         tSmall,fnval(sp,tSmall),'go')
    addbeats(D,'prey') 
    subplot(3,1,3); 
    plot(D.t,fnval(Dsp,D.t),'-',tHd,fnval(Dsp,tHd),'+k',...
         tBig,fnval(sp,tBig),'ro',tSmall,fnval(sp,tSmall),'go')
    addbeats(D,'prey') 
end


function D = find_pred_intervals(D)
% Algorithm for determining the time intervals for tail beats and tail
% flicks.  These are saved in the 'tBeat' field of D in a 3xn matrix, with
% the first column denoting tail beats (1) and flicks (0).

% Current smoothing spline fit
sp      = D.sp.angPd;
Dsp     = fnder(sp,1);
D2sp    = fnder(sp,2);

% Time of changes in sign of rate of change
tHd = fnzeros(D2sp);
tHd = unique(tHd(1,:)');

% Time of peaks in angular rate of change
tPeak = fnzeros(D2sp);
tPeak = unique(tPeak(1,:)');

% Time of big tail beats (candidate)
tBig = tPeak(abs(fnval(Dsp,tPeak))>D.p.yawRate_high);

% Check that tBig is nonempty
if isempty(tBig)
    % tBig is empty, use a smaller threshold value
    disp('tBig is empty...')
    tBig = tPeak(abs(fnval(Dsp,tPeak))>1.5);
end

% Start of changes in angular direction before peak
tPrior = tHd(tHd<tBig(1));

% Ending of changes in angular direction after peak
tFollow = tHd(tHd>tBig(1));

% Set starting time for first tail beat
if isempty(tPrior)
    tBeat = tBig(1);
else
    tBeat = tPrior(end);
end

% Set time for end of first tail beat
if isempty(tFollow)
    error('No times following!')
else
    tBeat(1,2) = tFollow(1);
end

% Loop thru canidate times
for i = 2:length(tBig)
    
    % Starts of changes in prior to current iBig
    tPrior = tHd((tHd<tBig(i)) & (tHd>=(tBig(i-1)+D.p.min_dur)));
    
    % Ending of changes in angular direction after current peak
    tFollow = tHd(tHd>tBig(i));
    
    % If next candidate is beyond the min duration . . .
    if ~isempty(tPrior) && ~isempty(tFollow)
        tBeat = [tBeat; tPrior(end) tFollow(1)];
    end
end

tBeat2 = [];

% Loop thru candidate beat times, eliminate small turns
for i = 1:size(tBeat,1)
    
    if size(tBeat,1)==1
        iP = D.t<tBeat(1,1);
        iF = D.t>tBeat(1,2);
        
    elseif i==1
        iP = D.t<tBeat(1,1);
        iF = (D.t>tBeat(1,2)) & (D.t<tBeat(2,1));
        
    else
        % Set prior to previous follow
        iP = iF;
        
        % If no other beat, follow is remaining time
        if i==size(tBeat,1)
            iF = D.t>tBeat(i,2);
            
            % Otherwise, follow stops at next beat
        else
            iF = (D.t>tBeat(i,2)) & (D.t<tBeat(i+1,1));
        end
    end
    
    % Change in heading
    dAng = mean(fnval(sp,D.t(iF))) - mean(fnval(sp,D.t(iP)));
    
    % Store if beat created above-threshold change in heading
    if abs(dAng) > D.p.yawDelta
        tBeat2 = [tBeat2; tBeat(i,:)];
    end
end

% Redefine tBeat
tBeat = tBeat2; clear tBeat2

% Time of small tail beats (candidate)
tSmall = tPeak((abs(fnval(Dsp,tPeak))>D.p.yawRate_low));
tFlick = [];

% Loop thru canidate times
for i = 1:length(tSmall)
    
    if i>1
        % Starts of changes in angular direction prior to current iBig
        tPrior = tHd((tHd<tSmall(i)) & (tHd>=(tSmall(i-1)+D.p.min_dur/4)));
    else
        tPrior = tHd(tHd<tSmall(1));
    end
    
    % Ending of changes in angular direction after current peak
    tFollow = tHd(tHd>tSmall(i));
    
    if ~isempty(tPrior) && ~isempty(tFollow)
        
        % Starting and end times for current tailbeat
        tStart = tPrior(end);
        if length(tFollow)>1
            tEnd = tFollow(2);
        else
            tEnd = tFollow(1);
        end
        
        % Says if start is within any tailbeat
        inBeat1 = max((tStart<=(tBeat(:,2)+D.p.min_dur)) & ...
            (tStart>=(tBeat(:,1))));
        
        % Says if end (+min_dur) is within tailbeat
        inBeat2 = max(((tEnd+D.p.min_dur)<=tBeat(:,2)) & ...
            ((tEnd+D.p.min_dur)>=(tBeat(:,1))));
        
        % Says if end is within tailbeat
        inBeat3 = max((tEnd<=tBeat(:,2)) & (tEnd>=(tBeat(:,1))));
        
        % Skip logging, if within a tail beat
        if ~inBeat1 && ~inBeat2 && ~inBeat3
            
            tFlick = [tFlick; tStart tEnd];
        end
    end
end

% Remove overlapping flicks
if size(tFlick,1)>1
    tFlick2 = tFlick(1,:);
    for i = 2:length(tFlick)
        if ~(tFlick(i,1)<tFlick(i-1,2))
            tFlick2 = [tFlick2; tFlick(i,:)];
        end
    end
    
    tFlick = tFlick2; clear tFlick2
end

tmp   = [ones(size(tBeat,1),1) tBeat];
tmp  = [tmp;[zeros(size(tFlick,1),1) tFlick]];

[tmp2,idx] = sort(tmp(:,2));

% Store in 'D'

D.tBeat = tmp(idx,:);


function vis_beats(D,fish)

% Current smoothing spline fit
if strcmp(fish,'pred')
    sp      = D.sp.angPd;
    
    % Predator speed
    spd = give(D,'pred spd');
    
    % Get beats and flicks
    tBeat  = D.tBeat(D.tBeat(:,1)==1,2:3);
    tFlick = D.tBeat(D.tBeat(:,1)==0,2:3);
else
    sp      = D.sp.angPy;
    
    % Prey speed
    spd = give(D,'prey spd');
    
    % Get beats and flicks
    tBeat  = D.tBeatPy(D.tBeatPy(:,1)==1,2:3);
    tFlick = D.tBeatPy(D.tBeatPy(:,1)==0,2:3);
end

Dsp     = fnder(sp,1);
D2sp    = fnder(sp,2);



figure
subplot(5,1,1:2) %--------------------------
fnplt(sp);grid on;yL=ylim;hold on;
ylabel('Heading (rad)')

addbeats(D,'pred')

subplot(5,1,3) %--------------------------
fnplt(Dsp);hold on; 
ylabel('Dheading (rad/s)')
plot(tBeat,fnval(Dsp,tBeat),'k+');yL=ylim;grid on

addbeats(D,'pred')

subplot(5,1,4) %--------------------------
plot(D.t,spd,'-');
ylabel('spd (cm/s)')
grid on;yL=ylim;hold on;

addbeats(D,'prey')

subplot(5,1,5) %--------------------------
plot(D.t,D.bend,'-');
ylabel('Bending (a.u.)')
grid on;yL=ylim;hold on;


function addbeats(D,fish)
% Overlays periods of tail beats and tail flicks

if strcmp(fish,'pred')
    % Get beats and flicks
    tBeat  = D.tBeat(D.tBeat(:,1)==1,2:3);
    tFlick = D.tBeat(D.tBeat(:,1)==0,2:3);
    clr1 = 'r';
elseif strcmp(fish,'prey')
    % Get beats and flicks
    tBeat  = D.tBeatPy(D.tBeatPy(:,1)==1,2:3);
    tFlick = D.tBeatPy(D.tBeatPy(:,1)==0,2:3);
    clr1 = 'g';
end

grid on;
yL=ylim;
hold on;

for i=1:size(tBeat,1),
    h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],clr1);
    set(h,'EdgeColor','none'); alpha(h,0.2)
end
for i=1:size(tFlick,1),
    h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
    set(h,'EdgeColor','none'); alpha(h,0.1)
end
hold off


function R = local_system(origin,rost)

% Check dimensions
if size(origin,1)~=1 || size(origin,2)~=2 || size(rost,1)~=1 || size(rost,2)~=2 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = rost(1) - origin(1);
xaxis(1,2) = rost(2) - origin(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% Format for matlab
tform = affine2d(R);


function tform = local_system_R(eyeR,eyeL)
% Defines a local coordinate system for the right eye where the y-axis
% points anteriorly and the x-axis points to the right of the body and the
% origin is at the center of the right eye. NOTE: Does not take into
% account the orientation of the eye.

% Right eye serves as origin
origin = eyeR;

% Define R eye relative to L eye
xaxis(1,1) = eyeR(1) - eyeL(1);
xaxis(1,2) = eyeR(2) - eyeL(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% Format for matlab
tform = affine2d(R);


function tform = local_system_L(eyeR,eyeL)
% Defines a coordinate system for the left eye where the y-axis points
% posteriorly and is aligned with the heading. The x-axis points toward the
% left side of the body.  NOTE: Does not take into account the orientation
% of the eye.

% Left eye serves as origin
origin = eyeL;

% Define L eye relative to R eye
xaxis(1,1) = eyeL(1) - eyeR(1);
xaxis(1,2) = eyeL(2) - eyeR(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% Format for matlab
tform = affine2d(R);


function [xT,yT] = global_to_local(tform,x,y)
% Assumes column vectors for coordinates

% Check dimensions
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end

pts = [x y];

% Translate
pts(:,1) = pts(:,1) - tform.T(3,1);
pts(:,2) = pts(:,2) - tform.T(3,2);

% Rotate points
ptsT = [tform.T(1:2,1:2) * pts']';

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);

function [xT,yT] = local_to_global(tform,x,y)
% Assumes columns vectors for coordinates

% Check dimensions
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end

% Loop thru columns of coordinates
for i = 1:size(x,2)
    
    pts = [x(:,i) y(:,i)];
    
    % Rotate points
    ptsT = (tform.T(1:2,1:2) \ pts')';
    
    % Translate global coordinates wrt origin
    ptsT(:,1) = ptsT(:,1) + tform.T(3,1);
    ptsT(:,2) = ptsT(:,2) + tform.T(3,2);
    
    % Extract columns of points
    xT(:,i) = ptsT(:,1);
    yT(:,i) = ptsT(:,2);
    
    clear ptsT pts
end

