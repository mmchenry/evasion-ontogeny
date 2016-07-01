function anaPool
% Collects data that was created by anaData


%% Code execution

% Run anaData on all eligible sequences
run_anaData = 1;

% Force anaData to be re-run on all eligible sequences
rerun_anaData = 1;


%% Path definitions

% Give general path definitions
paths = givePaths;

% List of batches
batches = dir([paths.data filesep '20*']);



%% Run anaData to create merged dataset

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

return

%% Plot data

% Loop thru batches
for i = 1:length(batches)
    
    % List of sequences
    seqs = dir([paths.data filesep batches(i).name filesep 'S*']);
    
    % Loop trhu experiments
    for j = 1:length(seqs)
        
        % Directory for current data
        dPath = [paths.data filesep batches(i).name filesep seqs(j).name];
        
        % If turn data exists (i.e. anaData has been run) . . .
        if ~isempty(dir([dPath filesep 'turn data.mat'])) 
            
            % Title text for plots
            title_txt = [batches(i).name ': ' seqs(j).name];
            
            % Load 'sp' and 'D' for sequence
            load([dPath filesep 'turn data.mat'])
            
            % Load predator midline data ('mid')
            load([dPath filesep 'midline data.mat'])
            
            % Load eye & heading data ('eyes')
            load([dPath filesep 'eye data.mat'])
            
            % Load prey data ('prey')
            load([dPath filesep 'prey data.mat'])
            
            % Check time vector presence
            if ~isfield(sp,'time')
                warning(['structure sp does not have "time" field in ' ...
                    batches(i).name ': ' seqs(j).name])
                
            elseif ~isfield(D,'tAV')
                warning(['structure D does not have "tAV" field in ' ...
                    batches(i).name ': ' seqs(j).name])
            else
                [batches(i).name filesep seqs(j).name]
                % Plot basic data from sequence
                if 0
                    % Plot trajectories and angle/speed data
                    plot_basics(sp,eyes,mid,prey,D,title_txt)
                    %close all
                end
                
            end
            
        end
    end
end






function  plot_basics(sp,eyes,mid,prey,D,title_txt)
% Plots trajectories of predator and prey

figure

% Time vector to evaluate spline
t = linspace(min(sp.time),max(sp.time),500);

subplot(6,2,[1:2:12])
% Prey trajectory
plot(fnval(sp.xPrey,t),fnval(sp.yPrey,t),'-','Color',[.73 .83 .96]);
hold on
plot(prey.xPrey,prey.yPrey,'b.')

% Predator trajectory
plot(fnval(sp.xRost,t),fnval(sp.yRost,t),'-','Color',.8.*[.8 .6 .6]);
plot(mid.xRost,mid.yRost,'.r')

% Plot first points
plot(fnval(sp.xRost,t(1)),fnval(sp.yRost,t(1)),'ro',...
    fnval(sp.xPrey,t(1)),fnval(sp.yPrey,t(1)),'bo');

% Overlay connecting lines for beat times
for k = 1:length(D.tAV)
    cPdX = fnval(sp.xRost,D.tAV(k));
    cPdY = fnval(sp.yRost,D.tAV(k));
    cPyX = fnval(sp.xPrey,D.tAV(k));
    cPyY = fnval(sp.yPrey,D.tAV(k));
    
    plot([cPdX cPyX],[cPdY cPyY],'-','Color',.8.*[1 1 1])
    
    clear cPdX cPdY cPyX cPyY
end

xlabel('x-position (pix)')
ylabel('y-position (pix)')
axis equal
hold off
title(title_txt)

% Generate values from splines
hdAngle     = fnval(sp.hdAngle,eyes.t);
hdX         = fnval(fnder(sp.xRost),sp.time);
hdY         = fnval(fnder(sp.yRost),sp.time);
hdSpd       = hypot(hdX,hdY);
pyX         = fnval(fnder(sp.xPrey),sp.time);
pyY         = fnval(fnder(sp.yPrey),sp.time);
pySpd       = hypot(pyX,pyY);


subplot(6,2,[2 4]) %------------------------

% Plot spline fit
plot(eyes.t,hdAngle * 180/pi,'-','LineWidth',1,...
    'Color',.6.*[1 1 1])
hold on

% Plot heading angle
plot(eyes.t,unwrap(eyes.hdAngle) * 180/pi,'r.')
yL = ylim;

% Overlay turn markers
for k = 1:length(D.tAV)
    plot(D.tAV(k).*[1 1],yL,'-','Color',.8.*[1 1 1]);
end

xlabel('t (s)')
ylabel('Head orientation (inertial coords, deg)')
title([title_txt ', Predator'])
hold off;

subplot(6,2,6) %------------------------

% Plot spline fit
plot(sp.time,hdSpd,'r-')
hold on
yL = ylim;

% Overlay turn markers
for k = 1:length(D.tAV)
    plot(D.tAV(k).*[1 1],yL,'-','Color',.8.*[1 1 1]);
end

xlabel('t (s)')
ylabel('Rostrum speed (pix/s)')
hold off

subplot(6,2,[8 10]) %------------------------

plot(prey.t,unwrap(prey.thetaPrey).*180/pi,'b-')
hold on
yL = ylim;

% Overlay turn markers
for k = 1:length(D.tAV)
    plot(D.tAV(k).*[1 1],yL,'-','Color',.8.*[1 1 1]);
end

xlabel('t (s)')
ylabel('Prey orientation (inertial coords, deg)')
title('Prey')
hold off

subplot(6,2,12) %------------------------


% Plot spline fit
plot(sp.time,pySpd,'b-')
hold on
yL = ylim;

% Overlay turn markers
for k = 1:length(D.tAV)
    plot(D.tAV(k).*[1 1],yL,'-','Color',.8.*[1 1 1]);
end

xlabel('t (s)')
ylabel('Prey speed (pix/s)')
hold off
