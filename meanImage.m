function imMean = meanImage(vDir,num_digit,maxFrames)
% Creates a mean image from a sequence of video frames
% vDir      - directory containing tif images for video
% num_digit - number of digits at the end of the filename for frame numbers
% maxFrames - maximum number of frames to use for mean image



%% Cue up the processing

% Get video information
M = videoInfo(vDir,num_digit);

% Define list of files to be used, depending on max number of frames
if length(M.frNums) > maxFrames
    dframe = floor(length(M.frNums)/maxFrames);
    frIdx = 1:dframe:length(M.frNums);
    clear dframe
else
    frIdx = 1:length(M.frNums);
end

% Create waitbar
h = waitbar(0,...
    ['Mean image: ' num2str(i)],...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% Create sum image based on first frame
[imCurr,tmp] = imread([vDir filesep M.filename{1}]);
imSum = double(imCurr);

clear imCurr tmp


%% Execute processing

% Loop through frames
for j = 1:length(frIdx)
    
    % Add current frame to sum image
    [imCurr,tmp] = imread([vDir filesep filesep M.filename{frIdx(j)}]);
    imSum        = imSum + double(imCurr);
    
    clear tmp imCurr
    
    % Update status bar
    h = waitbar(j/length(frIdx),h,...
        ['Mean image: ' num2str(j) ' of ' ...
        num2str(length(frIdx)) ' frames']);
    
    % Quit m-file, if cancel button pushed
    if getappdata(h,'canceling')
        close force
        return
    end
    
end

% Determine image bit depth
[imCurr,tmp] = imread([vDir filesep filesep M.filename{frIdx(j)}]);

% 16 bit
if strcmp(class(imCurr),'uint16')
    imMean = uint16(round(imSum./length(frIdx)));

% 8 bit
elseif strcmp(class(imCurr),'uint8')
    imMean = uint8(round(imSum./length(frIdx)));
    
else
    error('Cannot recognize bit depth of images')
end

% Calculate mean from sum image
imMean = imMean(:,:,1);

% Get rid of status bar
close force

