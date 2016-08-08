function anaPool
% Collects data that was created by anaData


%TODO: Check bearing calculation.
%TODO: Start to analyze the eye data.


%% Code execution

% Run anaData on all eligible sequences
run_anaData = 0;

% Force anaData to be re-run on all eligible sequences
rerun_anaData = 1;

% Force acqMaster to run anaFrames to get midline data
run_getmidline = 0; 

% Make comparisons between flicks and tail beats
run_compare = 0;

% Visualize timeseries data for all sequences
vis_timeseries = 0;

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
            
            % Load midine data, if present
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
            end
            
            % Reset runner for next sequence
            runner = 1;
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
end


%% Visualize predator-prey data

if vis_timepredprey
    execute_action(paths,batches,'vis_predprey(D,title_txt,k)');
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

function dOut = execute_action(paths,batches,action)
% Performs 'action' on all sequences found in the given batches

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
        if ~isempty(dir([dPath filesep 'merged data.mat']))
            
            % Title text for plots
            title_txt = [batches(i).name ': ' seqs(j).name];
            
            % Load 'D' for sequence
            load([dPath filesep 'merged data.mat'])
            
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


