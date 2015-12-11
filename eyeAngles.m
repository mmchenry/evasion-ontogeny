function S = eyeAngles(turnFrames)

close all

% get size of input vector
numTurn = length(turnFrames);

% check that input vector size is even 
if ~mod(numTurn,2)
    % continue
else
    error('Input array must be even')
end


%% Load data and setup variables

% dialog box to select the data file
[filename,path] = uigetfile('*.mat','Pick a MATLAB data file');

% load the data as a structure
data = load([path filesep filename]);

% rename fields
t       = data.bStats.t;
frame   = data.bStats.frame;
heading = data.bStats.angl;
Reye    = data.bStats.ReyeAngl;
Leye    = data.bStats.LeyeAngl;

% get total number of frames
N = length(frame);

% initialize index arrays
idxStart = zeros(1,numTurn);

% get indices for intervals between turns 
for j = 1:numTurn
    idxStart(j) = find(frame==turnFrames(j));
end

%% run some calculations

% compute average angles for first interval
S.Angl(1) = mean(heading(1:idxStart(1)))*180/pi;
S.Reye(1) = mean(Reye(1:idxStart(1)))*180/pi;
S.Leye(1) = mean(Leye(1:idxStart(1)))*180/pi;

% compute average angle for last interval
S.Angl(numTurn/2+1) = mean(heading(idxStart(numTurn):N))*180/pi;
S.Reye(numTurn/2+1) = mean(Reye(idxStart(numTurn):N))*180/pi;
S.Leye(numTurn/2+1) = mean(Leye(idxStart(numTurn):N))*180/pi;

for k = 1:(numTurn/2-1)  
    
    % compute average HEADING ANGLE during interval k
    S.Angl(k+1) = mean(heading(idxStart(2*k):idxStart(2*k+1)))*180/pi;
    
    % compute average EYE ANGLES during interval k
    S.Reye(k+1) = mean(Reye(idxStart(2*k):idxStart(2*k+1)))*180/pi;
    S.Leye(k+1) = mean(Leye(idxStart(2*k):idxStart(2*k+1)))*180/pi;
    
end



