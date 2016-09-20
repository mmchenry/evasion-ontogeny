function dOut = give(D,param)
% Returns variable values from raw data

% Default output is empty
dOut = [];


%% Indicies for interval before, during, and after a beat

if isfield(D,'tBeat')
    for i = 1:size(D.tBeat,1)
        
        %Log type
        beat_type(i,1) = D.tBeat(i,1);
        
        % Index for current beat
        iCurr{i} = (D.t>=D.tBeat(i,2)) & (D.t<=D.tBeat(i,3));
        
        % If first loop . . .
        if i==1
            iPre{i} = D.t<D.tBeat(i,2);
            
        % If after first . . .
        else
            iPre{i} = (D.t>D.tBeat(i-1,3)) & (D.t<D.tBeat(i,2));
        end
        
        % If last loop . . .
        if i==size(D.tBeat,1)
            iPost{i} = D.t>D.tBeat(i,3);
            
        % Else . . .
        else
            iPost{i} = (D.t>D.tBeat(i,3)) & (D.t<D.tBeat(i+1,2));
        end
    end
end


%% Return requested parameters
switch param

    case 'turn stats'
        
        spd = give(D,'pred spd');
        
        % Initial index
        j = 1;
        
        % Step thru all beats
        for i = 1:length(beat_type)

            if i==1 
                % Handle cases when there is no glide prior to a turn
                if isempty(D.t(iPre{i}))
                    % Change in heading
                    dOut(1,1) = mean(D.posPd(iPost{i},3));
                    % Max speed during turn
                    dOut(1,2) = max(spd(iCurr{i}));
                    % Duration of glide
                    dOut(1,3) = NaN;
                else
                    % Change in heading
                    dOut(1,1) = mean(D.posPd(iPost{i},3)) - mean(D.posPd(iPre{i},3));
                    
                    % Max speed during turn
                    dOut(1,2) = max(spd(iCurr{i}));
                    
                    % Duration of glide prior to beat
                    dOut(1,3) = range(D.t(iPre{i}));
                end
                % Advance index
                j = j + 1;
            else
                % Change in heading
                dOut(j,1) = mean(D.posPd(iPost{i},3)) - mean(D.posPd(iPre{i},3));
                
                % Max speed during turn
                dOut(j,2) = max(spd(iCurr{i}));
                
                try
                    % Duration of glide prior to beat
                    dOut(1,3) = range(D.t(iPre{i}));
                catch
                    % No glide in prior beat (fish not moving?)
                    dOut(1,3) = NaN;
                end
                
                % Advance index
                j = j + 1;
                
            end
        end
        
        
    case 'flick stats'
        
        spd = give(D,'pred spd');
        
        % Initial index
        j = 1;
        
        % Step thru all beats
        for i = 1:length(beat_type) 
            
            % If tail flick . . .
            if beat_type(i)==0
                
                % Log change in heading
                dOut(j,1) = mean(D.posPd(iPost{i},3)) - mean(D.posPd(iPre{i},3));
                
                % Log max speed
                dOut(j,2) = max(spd(iCurr{i}));
                
                % Advance index
                j = j + 1;
            end
        end
    
    case 'pred spd'
        % Rate fo change in x and y via splines
        dXsp = fnder(D.sp.xPd);
        dYsp = fnder(D.sp.yPd);
        
        % Resultant speed
        dOut = sqrt(fnval(dXsp,D.t).^2 + fnval(dXsp,D.t).^2);
        
        % Test against discrete calculation
        if 0, figure
            spd = sqrt(diff(D.posPd(:,1)).^2 + diff(D.posPd(:,2)).^2)./diff(D.t);
            plot(D.t,dOut,'-',D.t(2:end),spd,'.')
        end
        
    case 'prey spd'
        % Rate fo change in x and y via splines
        dXsp = fnder(D.sp.xPy);
        dYsp = fnder(D.sp.yPy);
        
        % Resultant speed
        dOut = sqrt(fnval(dXsp,D.t).^2 + fnval(dXsp,D.t).^2);
        
        % Test against discrete calculation
        if 0, figure
            spd = sqrt(diff(D.posPy(:,1)).^2 + diff(D.posPy(:,2)).^2)./diff(D.t);
            plot(D.t,dOut,'-',D.t(2:end),spd,'.')
        end    
        
    case 'ang dev'  % Angular deviation between gaze and prey
        
        % Indices for positive (or 0) bearing angle; i.e., prey on left side
        indPos = D.posPd(:,3) >= 0;
        
        % Indices for negative bearing angle; i.e., prey on right side
        indNeg = ~indPos;
        
        % Preallocate vector
        xtheta_E = zeros(size(D.t));
        ytheta_E = zeros(size(D.t));
        theta_E  = zeros(size(D.t));
        
        % Coordinates of vector from eye to prey (for left side/left eye)
        xtheta_E(indPos) = D.posPy(indPos,1) - D.posL(indPos,1);
        ytheta_E(indPos) = D.posPy(indPos,2) - D.posL(indPos,2);
        
        % Coordinates of vector from eye to prey (for right side/right eye)
        xtheta_E(indNeg) = D.posPy(indNeg,1) - D.posR(indNeg,1);
        ytheta_E(indNeg) = D.posPy(indNeg,2) - D.posR(indNeg,2);
        
        % Direction of angle from eye to prey (inertial coordinates)
        theta_EG = atan2(ytheta_E, xtheta_E);
        
        % Gaze angle values
        gazeR = give(D,'gazeR');
        gazeL = give(D,'gazeL');
        
        % Angle between gaze and vector from pred to prey: theta_E
        theta_E(indPos) = theta_EG(indPos) - gazeL(indPos);
        theta_E(indNeg) = theta_EG(indNeg) - gazeR(indNeg);
        
        % Spline fit theta_E
        [sp,dOut] = spaps(D.t,theta_E,D.sp.tol.Head);
        
        % Spline fit theta_EG
        %[sp,dOut] = spaps(prey.t,theta_EG,D.sp.tol.Head);
        
    case 'ang size' % Angular size of prey (angle subtended on retina)
 
        % Distance between prey & position of eye
        distPreyEye = sqrt(sum([D.posPy(:,1)-D.posPd(:,1) ...
                                D.posPy(:,2)-D.posPd(:,2)].^2,2));
        
        % Angular size of prey (angle subtended on retina)
        delta = 2*atan(D.lenPy./(2*distPreyEye));
        
        % Save angular size of prey
        [sp,dOut] = spaps(D.t,delta,D.sp.tol.Head*5);

    case 'eye angs'
 
        % Current angular eye positions
        angR = D.posR(:,3);
        angL = D.posL(:,3);
        
        angR2 = angR - mean([angR; angL]);
        angL2 = angL - mean([angR; angL]);
        
        if 0
            plot(D.t,angR*180/pi,'--',D.t,angL*180/pi,'--',...
                 D.t,angR2*180/pi,'-',D.t,angL2*180/pi,'-');
            grid on; 
            
            ttt=2;
        end
        
        % Output
        dOut = [angR2 angL2];      
        
    case 'bearing unwrapped' % Bearing angle
        
        % Calculate unwrapped data
        dOut = give(D,'bearing');
        
        % Loop until all values greater than pi pr less than -pi are gone
        while true
            
            % Index for outside values
            idxPos = dOut>pi;
            idxNeg = dOut<-pi;
            
            if max(idxPos)==0 && max(idxNeg)==0
                break
            end
            
            dOut(idxPos) = dOut(idxPos) - 2*pi;
            dOut(idxNeg) = 2*pi + dOut(idxNeg);
        end
        
    case 'bearing' % Bearing angle
        
        % x-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        x_R = D.posPy(:,1) - D.posPd(:,1);
        
        % y-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        y_R = D.posPy(:,2) - D.posPd(:,2);
        
        % Direction of range/baseline vector
        alpha_R = atan2(y_R, x_R);
        
        % Magnitude of range/baseline vector (distance between prey & pred)
        rangeMag = sqrt(sum([x_R, y_R].^2,2));
        
        % Bearing Angle
        %phi = unwrap(alpha_R) - unwrap(D.posPd(:,3));
        dOut = unwrap(alpha_R) - unwrap(D.posPd(:,3));
        
        % Spline fit bearing angles (with heading tolerance)
        %[sp,dOut] = spaps(D.t,phi,D.sp.tol.Head);
        
    case 'distPredPrey' % Distance btwn predator and prey
        
        % x-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        x_R = D.posPy(:,1) - D.posPd(:,1);
        
        % y-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        y_R = D.posPy(:,2) - D.posPd(:,2);
        
        % Direction of range/baseline vector
        alpha_R = atan2(y_R, x_R);
        
        % Magnitude of range/baseline vector (distance between prey & pred)
        %rangeMag = sqrt(sum([x_R, y_R].^2,2));
        dOut = sqrt(sum([x_R, y_R].^2,2));
        
        % Save distance between prey & pred
        %[sp,dOut] = spaps(D.t,rangeMag,D.sp.tol.Head);
    
    case 'angR' % ANgular position of right eye
        
        dOut = fnval(D.sp.angR,D.t);
        
        
    case 'angL' % Angular position of left eye
        
        dOut = fnval(D.sp.angL,D.t);    
        
    case 'angular comparison' % Angular position of prey in R eye
        
        % Initial index
        j = 1;
        
        % Indicies for turning tail beats that don't have motion during glide
        iMoving = glideOverlap(D);
        
        % Get adjusted values for eye angle
        tmp = give(D,'eye angs');
        angsR = tmp(:,1);
        angsL = tmp(:,2);
        
        % Step thru all beats
        for i = 1:length(beat_type) 

            % If tail beat . . .
            if beat_type(i)==1
                
                % Current time values
                t = D.t(iPost{i});
                              
                % Current angular eye positions
                angR = angsR(iPost{i});
                angL = angsL(iPost{i});
                
                % Current angular positions of prey
                angPyR = D.posPyR(iPost{i},3);
                angPyL = D.posPyL(iPost{i},3);
                
                % Visual range for right & left eyes
                rangeR = [(pi/2+D.p.vergAng/2)-D.p.fov (pi/2+D.p.vergAng/2)];
                rangeL = -[rangeR(2) rangeR(1)];
                
                % If always in the FOV of either eye . . .
                if min((angPyR>=rangeR(1) & angPyR<=rangeR(2)) | ...
                        (angPyL>=rangeL(1) & angPyL<=rangeL(2)))
                    
                    % If prey is closer to the diretcion of R eye
                    if mean(abs(angPyR)) < mean(abs(angPyR))
                        % Use angles from right
                        angC   = unwrap(angR);
                        angPyC = unwrap(angPyR);
                        
                    % Otherwise, use the left
                    else
                        angC   = unwrap(angL);
                        angPyC = unwrap(angPyL);
                    end
                        
                    % Fit to change prey angular position wrt both eyes
                    cPy = polyfit(t,angPyC,1);
                    
                    % Fit to change in eye angles
                    cAng = polyfit(t,angC,1);
                    
                    % Change in prey and eye
                    Dpy  = polyval(cPy,t(end)) - polyval(cPy,t(1));
                    Dang = polyval(cAng,t(end)) - polyval(cAng,t(1));    
 
                    % Visual check on angular values
                    if 0            
                        subplot(2,1,1)
                        plot(t,angPyC*180/pi,'ro',t,polyval(cPy,t)*180/pi,'r-');
                        title('Angular prey position (rad)')
                        %legend('R','R fit','L','L fit')
                        
                        subplot(2,1,2)
                        plot(t,angC*180/pi,'ro',t,polyval(cAng,t)*180/pi,'r-')
                        title('Eye angle (rad)')
                        ttt=1;
                        %legend('R','R fit','L','L fit')
                    end
                    
                    % Store whether prey are moving, change in angles
                    dOut(j,:) = [iMoving(i) Dpy Dang];

                    % Advance index
                    j = j + 1;
                    
                end
            end
        end

         
    case 'gazeR' % Gaze angle (right eye)
        % Eye angles (world coordinates: [-pi pi])
        ReyeG = (D.posPd(:,3) + D.posR(:,3));
        
        % Get negative angles
        negR = ReyeG < 0;
        
        % Rewrite eye angles to be in [0 2Pi]
        ReyeG(negR) = 2*pi + ReyeG(negR);
        
        % Compute gaze angle (world coordinates: [0 2Pi])
        gazeR = (ReyeG - pi/2);
        
        % Spline fit gaze angles
        [sp,dOut]   = spaps(D.t,gazeR,D.sp.tol.Head);
        
    case 'gazeL' % Gaze angle (left eye)
        % Eye angles (world coordinates: [-pi pi])
        LeyeG = (D.posPd(:,3) + D.posL(:,3));
        
        % Get negative angles
        negL = LeyeG < 0;
        
        % Rewrite eye angles to be in [0 2Pi]
        LeyeG(negL) = 2*pi + LeyeG(negL);
        
        % Compute gaze angle (world coordinates: [0 2Pi])
        gazeL = (LeyeG + pi/2);
        
        % Spline-fit gaze angles
        [sp,dOut] = spaps(D.t,gazeL,D.sp.tol.Head);
        
    case 'bearingPre still' % bearing angle at start of turn (motionless prey)
        
        % Bearing values
        bearing  = give(D,'bearing unwrapped');
        
        % Turn stats
        stats = give(D,'turn stats');
        
        % Indicies for turning tail beats that d
        idx = D.tBeat(:,1)==1 & ~glideOverlap(D);
        
        % Get start (of turn) times 
        tStart = D.tBeat(idx,2);
        
        % Output bearing and change in heading
        dOut = [interp1(D.t,bearing,tStart) stats(idx,1)];

    case 'bearingPre moving' % bearing angle at start of turn (motionless prey)
        
        % Bearing values
        bearing  = give(D,'bearing unwrapped');
        
        % Turn stats
        stats = give(D,'turn stats');
        
        % Indicies for turning tail beats that d
        idx = D.tBeat(:,1)==1 & glideOverlap(D);
        
        % Get start (of turn) times 
        tStart = D.tBeat(idx,2);
        
        % Output bearing and change in heading
        dOut = [interp1(D.t,bearing,tStart) stats(idx,1)];
        
    case 'beat duration still' % bearing angle at start of turn (motionless prey)
               
        % Turn stats
        stats = give(D,'turn stats');
        
        % Indicies for turning tail beats that don't have motion during glide
        idx = D.tBeat(:,1)==1 & ~glideOverlap(D);
        
        % Get duration of glides
        dOut = stats(idx,3);       
        
    case 'beat duration moving' % bearing angle at start of turn (motionless prey)
               
        % Turn stats
        stats = give(D,'turn stats');
        
        % Indicies for turning tail beats that don't have motion during glide
        idx = D.tBeat(:,1)==1 & glideOverlap(D);
        
        % Get duration of glides
        dOut = stats(idx,3);          
        
    case 'hdDelta' % Change in orientation during turn
        dOut  = diff(fnval(D.sp.angPd,D.tInt),1,2);
        
    case 'hdDelta_pre' % Change in heading in prior interval (proxy for coasting)
        dOut  = diff(fnval(D.sp.angPd,D.priorInt),1,2);
        
    case 'gazeDelta' % Change in gaze angle in prior interval
        dOut       = diff(fnval(D.sp.gazeR,D.priorInt),1,2);
        
    case 'bearDelta' % Change in bearing angle in prior interval
        dOut       = diff(fnval(D.sp.bearAngl,D.priorInt),1,2);
        
    case 'alphaDelta' % Change in range vector direction (alpha) in prior interval
        dOut      = diff(fnval(D.sp.alpha,D.priorInt),1,2);
        
    case 'thetaE_Delta' % Change in gaze/prey deviation (theta_E) vector in prior interval
        dOut      = diff(fnval(D.sp.thetaE,D.priorInt),1,2);
        
    case 'thetaEG_Delta' % Change in gaze/prey deviation (theta_EG) vector in prior interval
        dOut      = diff(fnval(D.sp.thetaEG,D.priorInt),1,2);
        
    case 'thetaE_D1Delta'% Change in derivative of theta_E in prior interval
        dOut      = diff(fnval(fnder(D.sp.thetaE),D.priorInt),1,2);
        
    case 'thetaE_Pre' % value of theta_E at start of turn
        dOut        = fnval(D.sp.thetaE,D.tInt(:,1));
        
    case 'thetaEG_Pre' % value of theta_EG at start of turn
        dOut        = fnval(D.sp.thetaEG,D.tInt(:,1));
        
    case 'thetaE_D1Pre' % derivative of theta_E at start of turn
        dOut        = fnval(fnder(D.sp.thetaE),D.tInt(:,1));
        
    case 'bearingPre' % bearing angle at start of turn
        dOut      = fnval(D.sp.bearAngl,D.tInt(:,1));
        
    case 'bearingPost' % bearing angle at end of turn
        dOut      = fnval(D.sp.bearAngl,D.tInt(:,2));
        
    case 'bearD1Pre' % angular velocity (derivative of bearing) at start of turn
        dOut       = fnval(fnder(D.sp.bearAngl),D.tInt(:,1));
        
    case 'deltaPre' % Angular size of prey at start of turn
        dOut        = fnval(D.sp.delta,D.tInt(:,1));
        
    case 'deltaPost' % Angular size of prey at end of turn
        dOut       = fnval(D.sp.delta,D.tInt(:,2));
        
    case 'delta_DeltaPre' % Change in angular size during previous turn
        dOut        = diff(fnval(D.sp.delta,D.priorInt),1,2);
        
    case 'distPre' % Distance between predator and prey at start of turn (cm)
        dOut         = fnval(D.sp.distPredPrey,D.tInt(:,1)) .* D.p.cF;
        
    case 'distPost' % Distance between predator and prey at end of turn (cm)
        dOut        = fnval(D.sp.distPredPrey,D.tInt(:,2)) .* D.p.cF;
        
    otherwise
        error('Parameter requested not recognized');
end


function idx = beatOverlap(D)
% Return indicies for tail beats (or flicks) that overlap with prey swimming

% Loop thru turning tailbeats
for i = 1:size(D.tBeat,1)
    
    % True for when the start of turn is within a scoot
    startOverlap = max(repmat(D.tBeat(i,2),size(D.tBeatPy,1),1)>...
        D.tBeatPy(:,2) & ...
        repmat(D.tBeat(i,2),size(D.tBeatPy,1),1)<...
        D.tBeatPy(:,3));
    
    % True for when the end of turn is within a scoot
    endOverlap = max(repmat(D.tBeat(i,3),size(D.tBeatPy,1),1)>...
        D.tBeatPy(:,2) & ...
        repmat(D.tBeat(i,3),size(D.tBeatPy,1),1)<...
        D.tBeatPy(:,3));
    
    % True for when the start of turn is before a scoot
    startBefore = repmat(D.tBeat(i,2),size(D.tBeatPy,1),1)<...
        D.tBeatPy(:,2);
    
    % True for when the end of turn is after a scoot
    endAfter = repmat(D.tBeat(i,3),size(D.tBeatPy,1),1)>...
        D.tBeatPy(:,3);
    
    idx(i,1) = startOverlap || endOverlap || ...
               max(startBefore & endAfter);
end


function idx = glideOverlap(D)
% Return indicies for period before a tail beat that overlaps with prey swimming

% Start time for first glide
tStart = D.t(1);

% Loop thru turning tailbeats
for i = 1:size(D.tBeat,1)
    
    % True for when a scoot starts within current glide
    startOverlap = max(...
        (D.tBeatPy(:,2) > repmat(tStart,size(D.tBeatPy,1),1)) & ...
        (D.tBeatPy(:,2) < repmat(D.tBeat(i,2),size(D.tBeatPy,1),1)));
    
    % True for when a scoot ends within current glide
    endOverlap = max(...
        (D.tBeatPy(:,3) > repmat(tStart,size(D.tBeatPy,1),1)) & ...
        (D.tBeatPy(:,3) < repmat(D.tBeat(i,2),size(D.tBeatPy,1),1)));
    
    % True of scoots starting before current glide
    startBefore = D.tBeatPy(:,2) < repmat(D.tBeat(i,2),size(D.tBeatPy,1),1);
    
    % True of scoots ending after current glide
    endAfter = D.tBeatPy(:,3) > repmat(D.tBeat(i,3),size(D.tBeatPy,1),1);
    
    % True if there is no glide before the beat
    if D.tBeat(i,2)==D.t(1)
        idx(i,1) = 1;
        
        % True if scoot starts or ends in glide or spans glide
    else
        idx(i,1) = startOverlap || endOverlap || max(startBefore & endAfter);
    end
    
    % Start time for next glide
    tStart = D.tBeat(i,3);
end
