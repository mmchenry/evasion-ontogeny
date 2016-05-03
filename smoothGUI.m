function sp = smoothGUI(sp,varNames,data)
%
% smoothGUI is an interactive GUI for selecting a good smoothing tolerance
% value for raw data traces.
%
% INPUTS:
%   1) 'sp' is a structure with spline parameters 
%   2) 'varNames' is a cell array of strings indicating variables
%   that need a tolerance parameter
%
% OUTPUTS:
%   'sp' is the data structure with field 'tol' which contains the
%   smoothing tolernace parameters set by the user. 

%% Interactively choose smoothing tolerance

% Instructions
disp('Selecting period before and after response')
disp('     up: more smoothing');
disp('     down: less smoothing');
disp('     right: increase threshold');
disp('     left: decrease threshold');
disp('     return: set value');
disp('     esc: quit');

% Get data fieldnames
dataNames = fieldnames(data);

% Make figure window
figure

% Loop thru each variable

for j = 1:length(varNames)
    
    % Indicator for setting tolerance (0 means tolerance is not yet set)
    tolset = 0;
    
    % Name of current variable
    currVar = varNames{j};
     
    % Check that data and smoothing tolerance match
    if ~strcmp(currVar,dataNames{j})
        disp(' The tolerance parameter does not match the data')
        return
    else
    end
    
    % Initial smoothing tolerance for current variable
    spTol = sp.tol.(currVar);
    
    % Initial value of threshold for angular velocity
    spThresh = sp.thresh;
    
    % Time vector
    t = sp.time;
    
    % Current data
    d = data.(dataNames{j});
    
    while ~tolset
        
        % Current smoothing spline fit
        spCurr  = spaps(t,d,spTol);
        Dsp     = fnder(spCurr,1);
        D2sp    = fnder(spCurr,2);
        
        % Roots of second derivative (time of peak/troughs)
        D2roots = fnzeros(D2sp);
        D2roots = D2roots(1,:)';
        
        % Plot raw data and overlay smoothed data
        subplot(2,1,1)
        
        % Plot raw data
        plot(t,d,'.');
        hold on;
        
        % Plot smooth data
        fnplt(spCurr)
        hold off;
        
        grid on
        xlabel('Time (s)')
        ylabel(currVar)
        
        subplot(2,1,2)
        fnplt(Dsp);
        set(gca,'ColorOrderIndex',1)
        hold on
        
        % Plot peaks and troughs
        plot(D2roots,fnval(Dsp,D2roots),'k+')
        
        % Plot the threshold
        line([t(1) t(end)],[spThresh spThresh])
        line([t(1) t(end)],[-spThresh -spThresh])
        hold off;
        
        drawnow
        
        % Store results
        sp.(varNames{j})    = spCurr;
        sp.tol.(currVar)    = spTol;
        
        % Store threshold
        if j==1
            sp.thresh = spThresh;
        else 
        end
        
        % Prompt for input
        [~,~,b] = ginput(1);
        
        % Up arrow
        if b==30
            spTol = spTol*1.25;
            
        % Down arrow
        elseif b==31
            spTol = spTol*0.75;
            
        % Right arrow 
        elseif b==29 & j==1
            spThresh = spThresh + 0.5;
            
        % Left arrow
        elseif b==28 & j==1
            spThresh = spThresh - 0.5;
            
        % esc
        elseif b==27
            disp('Quitting (data not saved)')
            return
            
        % Return
        elseif isempty(b)
            
            tolset = 1;

        end
    end
    
    
    clear Dsp D2sp
end

pause(1)

% Create the field 'tolSet', indicating that the tolerances have been set.
sp.tolSet = 1;

% disp(' ')
% disp(['Smoothing variable: ' num2str(i) ' of ' num2str(length(aRev)) ' completed'])

%         % Save sequence stats
%         save([sPath.data filesep 'Timing data.mat'],'dT');

end


% close