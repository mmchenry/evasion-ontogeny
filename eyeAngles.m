function S = eyeAngles(turnFrames)

close all

if nargin < 1
    numTurn = [];
    turns = 0;
else
    % get size of input vector
    numTurn = length(turnFrames);
    
    % check that input vector size is even
    if ~mod(numTurn,2)
        % continue
    else
        error('Input array must be even')
    end

end

%% Load data and setup variables

% dialog box to select the data file
[filename,path] = uigetfile('*.mat','Pick a MATLAB data file');

% load the data as a structure
data1 = load([path filesep filename]);

% rename fields
t       = data1.bStats.t;
frame   = data1.bStats.frame;
angl    = data1.bStats.angl;
Reye    = data1.bStats.ReyeAngl;
Leye    = data1.bStats.LeyeAngl;
m1      = data1.bStats.m1;

% get total number of frames
N = length(frame);

% dialog box to select the data file
[filename,path] = uigetfile('*.mat','Load eye points data file');

% load the data as a structure
data2 = load([path filesep filename]);

% rename variables
frame_num   = data2.frame_num;
x1_left     = data2.x1_left;
x1_right    = data2.x1_right;
x2_left     = data2.x2_left;
x2_right    = data2.x2_right;
y1_left     = data2.y1_left;
y1_right    = data2.y1_right;    
y2_left     = data2.y2_left; 
y2_right    = data2.y2_right; 


%% run some calculations

if turns
    
    % initialize index arrays
    idxStart = zeros(1,numTurn);
    
    % get indices for intervals between turns
    for j = 1:numTurn
        idxStart(j) = find(frame==turnFrames(j));
    end
    
    % compute average angles for first interval
    S.Angl(1) = mean(angl(1:idxStart(1)))*180/pi;
    S.Reye(1) = mean(Reye(1:idxStart(1)))*180/pi;
    S.Leye(1) = mean(Leye(1:idxStart(1)))*180/pi;
    
    % compute average angle for last interval
    S.Angl(numTurn/2+1) = mean(angl(idxStart(numTurn):N))*180/pi;
    S.Reye(numTurn/2+1) = mean(Reye(idxStart(numTurn):N))*180/pi;
    S.Leye(numTurn/2+1) = mean(Leye(idxStart(numTurn):N))*180/pi;
    
    for k = 1:(numTurn/2-1)
        
        % compute average HEADING ANGLE during interval k
        S.Angl(k+1) = mean(angl(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        
        % compute average EYE ANGLES during interval k
        S.Reye(k+1) = mean(Reye(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        S.Leye(k+1) = mean(Leye(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        
    end
else
    
end

%% compute and plot angles

% compute slope defining eye angles
m_right = (y2_right-y1_right)./(x2_right-x1_right);
m_left = (y2_left-y1_left)./(x2_left-x1_left);

% find common frame number values
[~,ia,ib] = intersect(frame,frame_num,'stable');

% extract values corresponding to indices 'ia'
frameNew    = frame(ia);
anglNew     = angl(ia);
m1New       = m1(ia);

% extract values corresponding to indices 'ib'
m_right = m_right(ib);
m_left  = m_left(ib);







