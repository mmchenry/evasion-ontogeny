function bStats = anaSeq(action,dPath,tPath,cPath,p)
% Runs analysis of a sequence


%% Paths & data loading 

% Load data for blobs (B)
load([dPath filesep 'blob data.mat']);

% Load calibration data ('cal')
%load([cPath filesep 'calibration data.mat'])

% Load filenames of frames
aT = dir([tPath  filesep '*.mat']);

if isempty(aT)
    error('No video frames found');
end


%% Prelim analysis

if strcmp(action,'prelim')
    % Make time vector
    t = [1:length(B)]' ./ p.framerate;
    
    % Initialize index
    k = 1;
    
    % Loop thru frames
    for i = 1:length(B)
        
        % Filename
        fName = [B(i).filename(1:(end-4)) '.mat'];
        
        % If thumbnail exists
        if ~isempty(dir([tPath filesep fName]))
            % Load blob data
            load([tPath filesep fName])
        else
            blob.xMid = nan;
        end
   
 %%%%%%%%------------------ Do eye computations here ----------
        % If not-nan data
        if ~isnan(blob.xMid)
            
            % Store time
            bStats.t(k,1) = t(i);
            
            % Store frame number
            bStats.frame(k,1) = B(i).fr_num;
            
            % Store angle
            bStats.angl(k,1) = atan2(blob.yMid(1)-blob.yMid(2),...
                                     blob.xMid(1)-blob.xMid(2));
            
            % Store distance btwn rostrum and tail
            bStats.bSpan(k,1) = hypot(blob.yMid(1)-blob.yMid(end),...
                                      blob.xMid(1)-blob.xMid(end));
            
%             % Store coordinate data
%             B(i).sMid       = blob.sMid;
%             B(i).xMid       = blob.xMid;
%             B(i).yMid       = blob.yMid;
%             B(i).xEye       = blob.xEye;
%             B(i).yEye       = blob.yEye;
%             B(i).roi_blob   = blob.roi_blob;
            
            % Advance index
            k = k + 1;
        else
            %B(i).xMid = nan;
        end
    end
    
    % Unwrap data
    bStats.angl = unwrap(bStats.angl);
    
    % Smoothing spline fit
    bStats.spAng = spaps(bStats.t,bStats.angl,p.tolAng);
    
    % Roots of second derivative (time of peak/troughs)
    D2roots = fnzeros(fnder(bStats.spAng,2),p.FS_period); 
    D2roots = D2roots(:);
    
    % Index of peaks within the window for a FS
    idx = abs(fnval(fnder(bStats.spAng),D2roots)) > p.alpha_thresh;     
    
    % If any fast starts
    if max(idx)
        
        % Time of peak angular rotation
        tFS = D2roots(find(idx,1,'first'));
        
        % Rerun second derivative for early roots
        D2roots = fnzeros(fnder(bStats.spAng,2),[0 p.FS_period(2)]); 
        D2roots = D2roots(:);
        
        afterRoot   = D2roots(D2roots>tFS);
        beforeRoot  = D2roots(D2roots<tFS);
        
        % Start and end times for the fast start
        bStats.tFS = [beforeRoot(end) afterRoot(1)];
     
        % Frame number of 
        tmp = abs(bStats.t - bStats.tFS(1));        
        bStats.frFS         = bStats.frame(find(tmp==min(tmp),1,'first'));
        
        tmp = abs(bStats.t - bStats.tFS(2));        
        bStats.frFS(1,2)    = bStats.frame(find(tmp==min(tmp),1,'first'));
        
    % Otherwise, no fast start
    else
        bStats.tFS = nan;
    end
    
    % Roots of second derivative
    %Droots = fnzeros(fnder(bStats.spAng,2));
    
    % Translate into real-world coordinates
    %seq.originRW = pointsToWorld(cal.cameraParams,cal.R,cal.t,seq.origin);
    
    
    if 0
        figure;
        subplot(2,1,1)
        % PLot to get y-axis limit
        plot(bStats.t,bStats.angl,'.');
        yL = ylim;
     
        tRange = bStats.tFS;
        
        % Fill in a region
        h = fill([tRange(1) tRange(2) tRange(2) tRange(1) tRange(1)],...
                 [yL(1) yL(1) yL(2) yL(2) yL(1)],'k');
        alpha(h,0.2)
        set(h,'EdgeColor','none')
        hold on
        fnplt(bStats.spAng);
        plot(bStats.t,bStats.angl,'.');
        hold off
        grid on
        xlabel('Time (s)')
        ylabel('Head angle (rad)')
        
        subplot(2,1,2)
        fnplt(fnder(bStats.spAng));
        grid on
        xlabel('Time (s)')
        ylabel('Rate of head rotation (rad/s)')
        
    end
    
    % Save data for blobs (B)
    save([dPath filesep 'blob data.mat'],'B');
    
    % Save sequence stats
    save([dPath filesep 'Sequence stats.mat'],'bStats');
    
end
    


function D = calcBlobStats(blob,p)

    if ~isnan(blob.sMid(1))
        D.angl(i,1) = atan2(blob.yMid(1)-blob.yMid(2),...
                            blob.xMid(1)-blob.xMid(2));
        D.kappa(i,1) = hypot(blob.yMid(1)-blob.yMid(end),...
                             blob.xMid(1)-blob.xMid(end));                
    else
        D.angl(i,1)     = nan;
        D.kappa(i,1)    = nan;
    end

    
% Unwrap data
D.angl = unwrap(D.angl);



function [ti,kappa] = calcTailCurve(B)

num_pts = length(B.xMid)*4;

x = B.xMid;
y = B.yMid;
t = B.sMid;
ti = linspace(t(1),t(end),num_pts)';

% Smoothing spline, calc coords, curvature (kappa)   
tol = 1.e-2;
sp = spaps(t',[x y]',tol);
dsp = fnder(sp);
dspt = fnval(dsp,ti);
ddspt = fnval(fnder(dsp),ti);
kappa = (abs(dspt(1,:).*ddspt(2,:)-dspt(2,:).*ddspt(1,:))./...
    (sum(dspt.^2)).^(3/2))';

if 0
    subplot(2,1,1)
    plot(B.xMid,B.yMid,'ko')
    axis equal
    hold on
    fnplt(sp)
    hold off
    
    subplot(2,1,2)
    plot(ti,kappa,'k-')
    pause(.1)
end