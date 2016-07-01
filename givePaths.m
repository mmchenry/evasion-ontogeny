function paths = givePaths

% % Matt's computer
if ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')]))
    
    % Directory root
    %root = '/Users/mmchenry/Documents/Projects/Ontogeny of evasion/Batch experiments';
    root     = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';
    vid_root = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';

    % % Alberto's MacMini (office)
elseif ~isempty(dir([filesep fullfile('Users','alberto','Documents','GitHub-SourceTree')]))
    % Directory root
    root     = '/Volumes/VisualPred/ZF_visuomotor';
    vid_root = '/Volumes/VisualPred/ZF_visuomotor';
        
else
    % % Alberto's MacBook (laptop)
    path = fullfile('Users','A_Soto','Documents','Home');
    if ~isempty(dir([filesep path]))
        
        % Directory root
        root     = '/Volumes/BackUp/ZF_visuomotor';
        vid_root = '/Volumes/BackUp/ZF_visuomotor';
    
    else
        error('This computer is not recognized')
    end
    
end

% Save root
paths.root = root;

% To raw video files 
paths.rawvid = [vid_root filesep 'Raw video'];

% To data files
paths.data = [root filesep 'Data'];
