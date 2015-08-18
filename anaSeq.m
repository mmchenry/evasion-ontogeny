function anaSeq(vPath,tPath,dPath,cPath,p,overwrite)
% Runs analysis of a sequence





%% Paths & data loading 

% Load data for blobs (B)
load([dPath filesep 'blob data.mat']);

% Load calibration data ('cal')
load([cPath filesep 'calibration data.mat'])


for i = 1:length(B)
    
    % Store origin points
    seq.origin(i,:) = B(i).origin;
    
end

% Translate into real-world coordinates
seq.originRW = pointsToWorld(cal.cameraParams,cal.R,cal.t,seq.origin);

ttt=3