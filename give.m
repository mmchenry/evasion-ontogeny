function dOut = give(D,param)
% Returns variable values from raw data


switch param
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
        phi = unwrap(alpha_R) - unwrap(D.posPd(:,3));
        
        % Spline fit bearing angles (with heading tolerance)
        [sp,dOut] = spaps(D.t,phi,D.sp.tol.Head);
        
    case 'distPredPrey' % Distance btwn predator and prey
        
        % x-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        x_R = D.posPy(:,1) - D.posPd(:,1);
        
        % y-coordinate of vector (Range/baseline vector) from pred Rost to prey COM
        y_R = D.posPy(:,2) - D.posPd(:,2);
        
        % Direction of range/baseline vector
        alpha_R = atan2(y_R, x_R);
        
        % Magnitude of range/baseline vector (distance between prey & pred)
        rangeMag = sqrt(sum([x_R, y_R].^2,2));
        
        % Save distance between prey & pred
        [sp,dOut] = spaps(D.t,rangeMag,D.sp.tol.Head);
        
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

% transpose
if size(dOut,2)>size(dOut,1)
    dOut = dOut';
end

% % Duration of turns (sec)
% tTurn   = diff(D.tInt,1,2);
%
% % Duration between turns (sec); exlude first interval
% interTurn(1,1) = 100;
% interTurn = [interTurn; diff(D.priorInt(2:end,:),1,2)];
%
% % matrix of all change in angle data
% allData = [hdDelta, hdDelta_pre, gazeDelta, bearDelta, alphaDelta, ...
%      thetaE_Delta, thetaEG_Delta, thetaE_D1Delta,...
%      thetaE_Pre, thetaEG_Pre, thetaE_D1Pre,...
%      bearingPre, bearingPost, bearD1Pre, ...
%      deltaPre, deltaPost, delta_DeltaPre, ...
%      distPre, distPost, tTurn, interTurn];