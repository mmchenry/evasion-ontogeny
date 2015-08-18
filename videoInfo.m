function M = videoInfo(vDir,num_digit)
% Returns information about video, stored as a series of tifs
% vDir        - directory containing tif files
% num_digit   - number of digits that specific the frame number at end of
%               filename

% Ending to filename for images
nameSuffix = 'tif';

% Load filenames for frames
a = dir([vDir filesep '*.' nameSuffix]);

% Check that not empty
if isempty(a)
    error(['No video frames found in: ' vDir]);
    
elseif length(a)==1
    error(['Only a single tif file in: ' vDir])
end

% Loop thru tif filenames
for j = 1:length(a)
    
    % Read filename for frame number
    frNum = str2num(a(j).name(end-num_digit-length(nameSuffix):...
                    end-length(nameSuffix)-1));
    
    % Check ordering
    if (j>1) && (frNum ~= M.frNums(j-1)+1)
        error('Frame numbers not consecutive')
    end
    
    % Store info
    M.frNums(j,1) = frNum;
    M.filename{j} = a(j).name;
    M.path{j} = [vDir filesep a(j).name];
end
