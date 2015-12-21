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
Reye_loc  = data1.bStats.ReyeAngl;
Leye_loc  = data1.bStats.LeyeAngl;
m1      = data1.bStats.m1;
anglBod = data1.bStats.anglBody;

% get total number of frames
N = length(frame);

% dialog box to select the data file
[filename,path] = uigetfile('*.mat','Load eye points data file');

% load the data as a structure
data2 = load([path filesep filename]);

% rename variables
frame_num   = data2.frame;
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
    S.Reye(1) = mean(Reye_loc(1:idxStart(1)))*180/pi;
    S.Leye(1) = mean(Leye_loc(1:idxStart(1)))*180/pi;
    
    % compute average angle for last interval
    S.Angl(numTurn/2+1) = mean(angl(idxStart(numTurn):N))*180/pi;
    S.Reye(numTurn/2+1) = mean(Reye_loc(idxStart(numTurn):N))*180/pi;
    S.Leye(numTurn/2+1) = mean(Leye_loc(idxStart(numTurn):N))*180/pi;
    
    for k = 1:(numTurn/2-1)
        
        % compute average HEADING ANGLE during interval k
        S.Angl(k+1) = mean(angl(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        
        % compute average EYE ANGLES during interval k
        S.Reye(k+1) = mean(Reye_loc(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        S.Leye(k+1) = mean(Leye_loc(idxStart(2*k):idxStart(2*k+1)))*180/pi;
        
    end
else
    
end

%% compute and plot angles

% change in y (to compute arctan)
y_right = (y2_right-y1_right);
y_left  = (y2_left-y1_left);

% change in x (to compute arctan)
x_right = (x2_right-x1_right);
x_left  = (x2_left-x1_left);

% right eye angles (world coordinates)
Reye_w = unwrap(atan2(y_right,x_right));

% left eye angles (world coordinates)
Leye_w = unwrap(atan2(y_left,x_left));

% find common frame number values between digitized & autotracking
[~,ia,ib] = intersect(frame,frame_num,'stable');

% extract autotracking values corresponding to indices 'ia'
frameNew    = frame(ia);
anglNew     = angl(ia);
m1New       = m1(ia);
anglBodNew  = anglBod(ia);

% extract digitized values corresponding to indices 'ib'
Reye_w  = Reye_w(ib);
Leye_w  = Leye_w(ib);

% smooth heading angles
angl_smooth = smooth(anglNew);

% plot smooth heading angles
plot(frameNew,angl_smooth*180/pi,'.','MarkerSize',10)
hold on;

% plot eye orienation in world coordinates
plot(frameNew,(Reye_w)*180/pi,'.','MarkerSize',10)
plot(frameNew,(Leye_w)*180/pi,'.','MarkerSize',10)

% plot eye orientation in local coordinates
plot(frameNew,(Reye_w - angl_smooth)*180/pi,'.','MarkerSize',10)
plot(frameNew,(Leye_w - angl_smooth)*180/pi,'.','MarkerSize',10)

Pi = pi;







