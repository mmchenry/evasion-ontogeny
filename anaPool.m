function anaPool
% Collects data that was created by anaData


%TODO: Check bearing calculation.


%% Code execution

% Run anaData on all eligible sequences
run_anaData = 0;

% Force anaData to be re-run on all eligible sequences
rerun_anaData = 0;

% Force acqMaster to run anaFrames to get midline data
run_getmidline = 0; 

% Make comparisons between flicks and tail beats
run_compare = 1;

% Compare targeting vs. spontaneous modes
run_modes = 1;

% Visualize timeseries data for all sequences
vis_timeseries = 0;

% Visualize eye data for all sequences
vis_eyetimeseries = 0;

% Visualize eye data for all sequences
vis_eyepool = 0;

% Visualize timeseries data for the predator and prey for all sequences
vis_timepredprey = 0;

% Visualize correlation btwn bearing before a turn and the change in
% heading
vis_bearing_vs_Dheading = 1;

% Visualize glide stats
vis_glidestats = 1;


%% Path definitions

% Give general path definitions
paths = givePaths;

% List of batches
batches = dir([paths.data filesep '20*']);


%% Run anaData

if run_anaData || rerun_anaData

k = 1; 
runner = 1;
runnerPred = 0;

% Loop thru batches
for i = 1:length(batches)
    
    % List of sequences
    seqs = dir([paths.data filesep batches(i).name filesep 'S*']);
    
    % Loop thru experiments
    for j = 1:length(seqs)
        
        % Directory for current data
        dPath = [paths.data filesep batches(i).name filesep seqs(j).name];
        
        % If "merge data" not present (i.e. anaData has not been run) . . .
        if isempty(dir([dPath filesep 'merged data.mat'])) || rerun_anaData
            
            % Title text for plots
            title_txt = [batches(i).name ': ' seqs(j).name];
            
            disp(' '); disp(['Starting ' title_txt])
            
            % Load midline data, if present
            if ~isempty(dir([dPath filesep 'midline data.mat']))
                % Load predator midline data ('mid')
                load([dPath filesep 'midline data.mat'])
            else
                disp('    Skippping anaData: no "midline data.mat"');
                runner = 0;
            end
            
            % Load eye data, if present
            if ~isempty(dir([dPath filesep 'eye data.mat']))
                % Load eye & heading data ('eyes')
                load([dPath filesep 'eye data.mat'])
            else
                disp('    Skippping anaData: no "eye data.mat"');
                runner = 0;
            end
            
            if ~isempty(dir([dPath filesep 'prey data.mat']))
                % Load prey data ('prey')
                load([dPath filesep 'prey data.mat'])
            else
                disp('    Skippping anaData: no "prey data.mat"');
                runner = 0;
            end
            
            % Check fields of eye data
            if ~isfield(eyes,'xReye')
                disp('    Skippping anaData: eye data incomplete');
                runner = 0;
            % Check if eyeData contains only heading data
            elseif isfield(eyes,'eyeData')
                if ~eyes.eyeData
                    runnerPred = 1;
                end
            end
            
            % Check fields of midline data
            if ~isfield(mid,'sMid')
                disp('    Note: No midline data');
                if run_getmidline
                    acqMaster(batches(i).name,seqs(j).name)
                    disp('    Note: Getting midline data');
                end
            end
            
            if runner
                % Run anaData
                D = anaData(mid,eyes,prey,batches(i).name);
                
                % Save D structure
                save([dPath filesep 'merged data.mat'],'D')
                
                % Report success
                disp('    Completed anaData!');
            elseif runnerPred
                if isempty(dir([dPath filesep 'merged data(pred).mat']))
                    
                    disp('   Running anaDataPred');
                    % Run anaDataPred
                    eyeData = eyes.eyeData;
                    D = anaDataPred(mid,eyes,batches(i).name,eyeData);
                
                    % Save D structure
                    save([dPath filesep 'merged data(pred).mat'],'D')
                
                    % Report success
                    disp('    Completed anaDataPred!')
                else
                    % load data ('D')
                    load([dPath filesep 'merged data(pred).mat'])
                end
                
            end
            
            % Reset runner (prey) for next sequence
            runner = 1;
            % Reset runner (pred) for next sequence
            runnerPred = 0;
        end
        

    end
end
end


%% Compare timing btwn turns for moving and stationary prey

if vis_glidestats
    
    nbins = 7;
    
    % Run ana_priorbeat for all sequences
    dOut = execute_action(paths,batches,'dOut = ana_durations(D,dOut);');
    
    % Unpack results
    dur_still      = dOut(dOut(:,1)==1,2);
    dur_moving     = dOut(dOut(:,1)==0,2);
    
    histogram(dur_still,nbins,'Normalization','pdf');
    hold on
    histogram(dur_moving,nbins,'Normalization','pdf');
    hold off
    xlabel('Duration of glide (s)')
    ylabel('PDF')
    
    legend('Still prey','Moving prey')
end


%% Analyze events before turn

if vis_bearing_vs_Dheading
    
    % Run ana_priorbeat for all sequences
    dOut = execute_action(paths,batches,'dOut = ana_priorbeat(D,dOut);');
    
    % Unpack results
    bear_still      = dOut(dOut(:,1)==1,2)*180/pi;
    Dhead_still     = dOut(dOut(:,1)==1,3)*180/pi;
    bear_move       = dOut(dOut(:,1)==0,2)*180/pi;
    Dhead_move      = dOut(dOut(:,1)==0,3)*180/pi;
    bear_all        = dOut(:,2)*180/pi;
    Dhead_all       = dOut(:,3)*180/pi;
    
    
    figure
    subplot(3,2,[1:4])
    [stats,slope,intercept] = reducedMajorAxis(bear_all,Dhead_all,1,0.05,1,...
        'All prey: ');
    xlabel('Bearing (deg)'); ylabel('Change in heading (deg)')
    
    subplot(3,2,5)
    [stats,slope,intercept] = reducedMajorAxis(bear_still,Dhead_still,1,0.05,1,...
        'Still prey: ');
    xlabel('Bearing (deg)'); ylabel('Change in heading (deg)')
    
    subplot(3,2,6)
    [stats,slope,intercept] = reducedMajorAxis(bear_move,Dhead_move,1,0.05,1,...
        'Moving prey: ');
    xlabel('Bearing (deg)'); ylabel('Change in heading (deg)')
end


%% Visualize timeseries data

if vis_timeseries
    execute_action(paths,batches,'vis_beats(D,title_txt)');
    if run_modes
        execute_action(paths,batches,'vis_beats(D,title_txt)',1);
    end
end


%% Visualize predator-prey data

if vis_timepredprey
    execute_action(paths,batches,'vis_predprey(D,title_txt,k)');
end


%% Visualize eye data

if vis_eyetimeseries
    dOut = execute_action(paths,batches,'vis_eyes(D,title_txt);');
end


%% Pool eye data

if vis_eyepool
    dOut = execute_action(paths,batches,'dOut = pool_eyes(D,dOut);');
    
    % Index of glides where prey moved
    iMove = dOut(:,1)==1;
    
    % Change in angular position of prey
    angPy = dOut(:,2)*180/pi;
    
    % Change in angular position of eyes
    angEye = dOut(:,3)*180/pi;
    
    subplot(1,2,1)
    plot(angPy(~iMove),angEye(~iMove),'o',...
        [min(angPy(~iMove)) max(angPy(~iMove))],...
        [min(angPy(~iMove)) max(angPy(~iMove))],'-');
        axis equal
        title('Motionless Prey')
        xlabel('Change in pos. of prey (deg)')
        ylabel('Change in pos. of eye (deg)')
        
    subplot(1,2,2)
    plot(angPy(iMove),angEye(iMove),'o',...
        [min(angPy(iMove)) max(angPy(iMove))],...
        [min(angPy(iMove)) max(angPy(iMove))],'-');
    axis equal
    title('Moving Prey')
    xlabel('Change in pos. of prey (deg)')
    ylabel('Change in pos. of eye (deg)')
end


%% Compare behaviors

if run_compare
    
    % Gather tailbeat stats
    t_stats = execute_action(paths,batches,...
        'dOut = [dOut; give(D,''turn stats'')];');
    
    % Gather tail flick stats
    f_stats = execute_action(paths,batches,...
        'dOut = [dOut; give(D,''flick stats'')];');
    
    % Plot results
    figure
    subplot(1,2,1)
    h1 = histogram(abs(t_stats(:,1))*180/pi,'Normalization','pdf');
    hold on
    h2 = histogram(abs(f_stats(:,1))*180/pi,'Normalization','pdf');
    hold off
    xlabel('Change in heading');
    ylabel('PDF')
    axis square
    
    subplot(1,2,2)
    h1 = histogram(t_stats(:,2),'Normalization','pdf');
    hold on
    h2 = histogram(f_stats(:,2),'Normalization','pdf');
    hold off
    xlabel('Max speed (?/s)');
    axis square
    
    legend('Turns','Flicks')  
end

%% Compare targeted vs. spontaneous modes

if run_compare && run_modes
    
    % Gather tailbeat stats for spontaneous swimming
    t_stats2 = execute_action(paths,batches,...
        'dOut = [dOut; give(D,''turn stats'')];',1);
    
    % Gather tail flick stats for spontaneous swimming
    f_stats2 = execute_action(paths,batches,...
        'dOut = [dOut; give(D,''flick stats'')];',1);
    
    % Plot results
    figure
    
    % Plot change in heading
    subplot(2,2,1)
    h1 = histogram(abs(t_stats(:,1))*180/pi,'Normalization','pdf');
    hold on
    h2 = histogram(abs(f_stats(:,1))*180/pi,'Normalization','pdf');
    hold off
    xlabel('Change in heading');
    ylabel('PDF')
    axis square
    title('Targeted Swimming')
    
    subplot(2,2,3)
    h1 = histogram(abs(t_stats2(:,1))*180/pi,'Normalization','pdf');
    hold on
    h2 = histogram(abs(f_stats2(:,1))*180/pi,'Normalization','pdf');
    hold off
    xlabel('Change in heading');
    ylabel('PDF')
    axis square
    title('Spontaneous Swimming')
    
    % Plot max speed
    subplot(2,2,2)
    h1 = histogram(t_stats(:,2),'Normalization','pdf');
    hold on
    h2 = histogram(f_stats(:,2),'Normalization','pdf');
    hold off
    xlabel('Max speed (?/s)');
    axis square
    
    legend('Turns','Flicks')  
    
    subplot(2,2,4)
    h1 = histogram(t_stats2(:,2),'Normalization','pdf');
    hold on
    h2 = histogram(f_stats2(:,2),'Normalization','pdf');
    hold off
    xlabel('Max speed (?/s)');
    axis square
    
    legend('Turns (spont.)','Flicks(spont.)')  
    
    % Plot glide durations
    figure
    histogram(t_stats(:,3),'Normalization','pdf');
    hold on
    histogram(t_stats2(:,3),'Normalization','pdf');
    hold off
    xlabel('Glide duration (s)');
    axis square
    
    legend('Targeted', 'Spontaneous')
end



function vis_eyes(D,title_txt)

% Get adjusted values for eye angle
tmp = give(D,'eye angs');
angsR = unwrap(tmp(:,1));
angsL = unwrap(tmp(:,2));

% Visual range for right & left eyes
rangeR = [(pi/2+D.p.vergAng/2)-D.p.fov (pi/2+D.p.vergAng/2)];
rangeL = [rangeR(2)-D.p.vergAng rangeR(2)-D.p.vergAng+D.p.fov];

rangeR = rangeR*180/pi;
rangeL = rangeL*180/pi;

figure
subplot(4,1,1)
plot(D.t,angsR*180/pi,'-')
addbeats(D,'pred');
xlabel('time');ylabel('R ang');grid on

subplot(4,1,2)
plot(D.t,angsL*180/pi,'-')
addbeats(D,'pred');
xlabel('time');ylabel('L ang');grid on

subplot(4,1,[3:4])
plot(D.t,unwrap(D.posPyR(:,3))*180/pi,'k-',...
     [D.t(1) D.t(end)],rangeR(1).*[1 1],'r--',...
     [D.t(1) D.t(end)],rangeR(2).*[1 1],'r--',...
     [D.t(1) D.t(end)],rangeL(1).*[1 1],'b--',...
     [D.t(1) D.t(end)],rangeL(2).*[1 1],'b--');
addbeats(D,'pred');
addbeats(D,'prey');
legend('R','L')
xlabel('time');ylabel('prey pos (R ang)');grid on

%TODO: Calculate the angular positon of prey realtive to predator &
%incorportate fov for both eyes

ttt=3;

function vis_predprey(D,title_txt,k)
% Visualizes kinematic data for the predator and prey

% Current smoothing spline fit for predator heading
%spPd      = D.sp.angPd;

% Predator speed
%spd_pd = give(D,'pred spd');

% Prey speed
spd_py = give(D,'prey spd');

figure
subplot(4,1,1) %----- PREDATOR HEADING
plot(D.t,D.posPd(:,3)*180/pi,'-')
addbeats(D,'pred')
ylabel('Pred heading (deg)')

title(title_txt)

subplot(4,1,2) %----- PREY speed
plot(D.t,spd_py,'-')
addbeats(D,'prey')
ylabel('Prey Speed (?/s)')

subplot(4,1,3) %----- BEARING 
plot(D.t,give(D,'bearing')*180/pi,'-');
addbeats(D,'pred')
ylabel('Bearing (deg)')

subplot(4,1,4) %----- DISTANCE 
plot(D.t,give(D,'distPredPrey')*100,'-');
addbeats(D,'pred')
addbeats(D,'prey')
ylabel('Distance (?)')

%preBearStill = give(D,'bearingPre still');

% Write images to downloads folder
if 0
    fname = ['000' num2str(k)];
    fname = fname(end-3:end);
    drawnow
    im = getframe(1);
    imwrite(im.cdata,['~/Downloads' filesep fname '.jpg'],'jpg');  
    close
end

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

function vis_beats(D,title_txt)
% Visualizes kinematic data for the predator

% Current smoothing spline fit
sp      = D.sp.angPd;
Dsp     = fnder(sp,1);
D2sp    = fnder(sp,2);

% Predator speed
spd_pd = give(D,'pred spd');

% Get beats and flicks
tBeat  = D.tBeat(D.tBeat(:,1)==1,2:3);
tFlick = D.tBeat(D.tBeat(:,1)==0,2:3);

figure
subplot(5,1,1:2) %--------------------------
fnplt(sp);grid on;yL=ylim;hold on;
ylabel('Heading (rad)')
for i=1:size(tBeat,1),
    h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
    set(h,'EdgeColor','none'); alpha(h,0.2)
end
for i=1:size(tFlick,1),
    h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
    set(h,'EdgeColor','none'); alpha(h,0.1)
end

title(title_txt)

subplot(5,1,3) %--------------------------
fnplt(Dsp);hold on; 
ylabel('Dheading (rad/s)')
plot(tBeat,fnval(Dsp,tBeat),'k+');yL=ylim;grid on
for i=1:size(tBeat,1),
    h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
    set(h,'EdgeColor','none'); alpha(h,0.2)
end
for i=1:size(tFlick,1),
    h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
    set(h,'EdgeColor','none'); alpha(h,0.1)
end

subplot(5,1,4) %--------------------------
plot(D.t,spd_pd,'-');
ylabel('Pred spd (?/s)')
grid on;yL=ylim;hold on;

for i=1:size(tBeat,1),
    h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
    set(h,'EdgeColor','none'); alpha(h,0.2)
end
for i=1:size(tFlick,1),
    h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
    set(h,'EdgeColor','none'); alpha(h,0.1)
end

subplot(5,1,5) %--------------------------
plot(D.t,D.bend,'-');
ylabel('Bending (a.u.)')
grid on;yL=ylim;hold on;

function dOut = pool_eyes(D,dOut)
% Visualizes eye kinematics for the predator

% Spline fit to predator heading
sp      = D.sp.angPd;

% Predator speed
spd_pd = give(D,'pred spd');

% Get beats and flicks
tBeat  = D.tBeat(D.tBeat(:,1)==1,2:3);
tFlick = D.tBeat(D.tBeat(:,1)==0,2:3);

tmp = give(D,'angular comparison');
dOut = [dOut; tmp];

% % Gaze of each eye
% gazeR = give(D,'gazeR');
% gazeL = give(D,'gazeL');
% 
% % Angular position of eyes
% angR = give(D,'angR');
% angL = give(D,'angL');

if 0
    figure
    subplot(5,1,1) %--------------------------
    fnplt(sp);grid on;yL=ylim;hold on;
    ylabel('Heading (rad)')
    for i=1:size(tBeat,1),
        h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
        set(h,'EdgeColor','none'); alpha(h,0.2)
    end
    for i=1:size(tFlick,1),
        h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
        set(h,'EdgeColor','none'); alpha(h,0.1)
    end
    
    title(title_txt)
    
    subplot(5,1,2) %--------------------------
    plot(D.t,spd_pd,'-');
    ylabel('Pred spd (?/s)')
    grid on;yL=ylim;hold on;
    
    for i=1:size(tBeat,1),
        h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
        set(h,'EdgeColor','none'); alpha(h,0.2)
    end
    for i=1:size(tFlick,1),
        h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
        set(h,'EdgeColor','none'); alpha(h,0.1)
    end
    
    subplot(5,1,3:4) %--------------------------
    h = plotyy(D.t,unwrap(angR)*180/pi,D.t,unwrap(angL)*180/pi);
    ylabel('Eye angle (deg)')
    legend('R','L')
    grid on;yL=ylim;hold on;
    
    for i=1:size(tBeat,1),
        h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
        set(h,'EdgeColor','none'); alpha(h,0.2)
    end
    for i=1:size(tFlick,1),
        h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
        set(h,'EdgeColor','none'); alpha(h,0.1)
    end
    
    subplot(5,1,5) %--------------------------
    plot(D.t,unwrap(gazeR)*180/pi,'-',D.t,unwrap(gazeL)*180/pi,'-');
    ylabel('Gaze (deg)')
    legend('R','L')
    grid on;yL=ylim;hold on;
    
    for i=1:size(tBeat,1),
        h = fill([tBeat(i,:) tBeat(i,2) tBeat(i,1)],[yL(1) yL(1) yL(2) yL(2)],'r');
        set(h,'EdgeColor','none'); alpha(h,0.2)
    end
    for i=1:size(tFlick,1),
        h = fill([tFlick(i,:) tFlick(i,2) tFlick(i,1)],[yL(1) yL(1) yL(2) yL(2)],'b');
        set(h,'EdgeColor','none'); alpha(h,0.1)
    end
end

function dOut = execute_action(paths,batches,action,mode)
% Performs 'action' on all sequences found in the given batches

% Set default mode as targeted swimming (mode=0)
if nargin < 4
    mode = 0;
end

% Set turn data name based on 'mode' indicator
if mode
    dataName = 'merged data(pred).mat';
else
    dataName = 'merged data.mat';
end

% Set empty default
dOut = [];

k = 1;

% Loop thru batches
for i = 1:length(batches)
    
    % List of sequences
    seqs = dir([paths.data filesep batches(i).name filesep 'S*']);
    
    % Loop thru experiments
    for j = 1:length(seqs)
        
        % Directory for current data
        dPath = [paths.data filesep batches(i).name filesep seqs(j).name];
        
        % If turn data exists (i.e. anaData has been run) . . .
        if ~isempty(dir([dPath filesep dataName]))
            
            % Title text for plots
            title_txt = [batches(i).name ': ' seqs(j).name];
            
            % Load 'D' for sequence
            load([dPath filesep dataName])
            
            % Execute the action
            eval(action)         
            
            k = k + 1;
        end
    end
end


function dOut = ana_durations(D,dOut)

% Data for prey staying still
d_still = give(D,'beat duration still');

% Data for moving prey
d_move = give(D,'beat duration moving');

% Log into dOut
dOut = [dOut; ones(size(d_still,1),1) d_still];
dOut = [dOut; zeros(size(d_move,1),1) d_move];


function dOut = ana_priorbeat(D,dOut)

% Data for prey staying still
d_still = give(D,'bearingPre still');

% Data for moving prey
d_move = give(D,'bearingPre moving');

% Log into dOut
dOut = [dOut; ones(size(d_still,1),1) d_still];
dOut = [dOut; zeros(size(d_move,1),1) d_move];


