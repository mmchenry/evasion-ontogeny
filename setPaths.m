function [vPath,tPath,dPath,cPath] = setPaths(paths,batchName,expName)

% Directory for current video
vPath = [paths.rawvid filesep batchName filesep expName];

% Directory for current thumbnails
tPath = [paths.thumb filesep batchName filesep expName];

% Directory for current data
dPath = [paths.data filesep batchName filesep expName];

% Directory path for calibration
cPath = [paths.cal filesep batchName];

% Make thumbnail directory, if not present
if isempty(dir(tPath))
    mkdir(tPath)
end

% Make data directory, if not present
if isempty(dir(dPath))
    mkdir(dPath)
end