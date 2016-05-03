function D = intervalsGUI(D,numTurns,hdAngle)
%
% intervalsGUI is an interactive GUI for adjusting the time intervals of
% major heading changes of a swimming fish
%
% INPUTS:
%   1) 'D' is a structure with a field 'tInt' of length numTurns. The
%   rows of this subfield contain the endpoints of a turn. 
%   2) 'numTurns' is the total number of turns for the current experiment
%   3) 'hdAngle' is the B-form of the heading angle
%
% OUTPUT:
%   'D' with user adjusted time intervals is the output

%% Interactively choose time intervals

% Instructions
disp('Selecting period before and after a turn')
% disp('     up: more smoothing');
% disp('     down: less smoothing');
disp('     right: move endpoint to the right');
disp('     left: move endpoint to the left');
disp('     return: set value');
disp('     esc: quit');

% Compute Angular velocity 
hdAngle_D1 = fnder(hdAngle);

% Get interval values (size is [2*numTurns, 1])
tInt = D.tInt(:);

% Make figure window
figure

for k = 1:2*numTurns
    
    % Indicator for setting tolerance (0 means tolerance is not yet set)
    intSet = 0;

    % Current time point
    endPnt = tInt(k);
    
    % Endpoints of current interval
    if k <= numTurns
        lEnd = tInt(k);
        rEnd = tInt(k+numTurns);
        
        % Set index for mode (left endpoints)
        iMode = 1;
    else
        lEnd = tInt(k-numTurns);
        rEnd = tInt(k);
        
        % Set index for mode (right endpoints)
        iMode = 2;
    end
    
    while ~intSet
  
        % Plot heading data
        subplot(2,1,1)
        fnplt(hdAngle);
        hold on;
        
        % Plot current time interval
        plot([lEnd rEnd], fnval(hdAngle,[lEnd rEnd]),'*r','MarkerSize',12)
        xlim([lEnd-0.1 rEnd+0.1])
        
        % Text description for interval
        intText = ['Interval = [' num2str([lEnd rEnd]) ']'];
        text(lEnd+0.01,fnval(hdAngle,lEnd),intText, 'FontSize',14);
        
        grid on
        xlabel('Time (s)')
        ylabel('Heading Angle')
        hold off;
        
        % Plot angular velocity 
        subplot(2,1,2)
        fnplt(hdAngle_D1);
        hold on;
        
        % Plot current time interval
        plot([lEnd rEnd], fnval(hdAngle_D1,[lEnd rEnd]),'*r','MarkerSize',12)
        xlim([lEnd-0.1 rEnd+0.1])
        
        % Text description for interval
        intText2 = ['Derivative = ' num2str(fnval(hdAngle_D1,endPnt))];
        text(endPnt,fnval(hdAngle_D1,endPnt),intText2, 'FontSize',14);
        
        grid on
        xlabel('Time (s)')
        ylabel('Heading angular velocity')
        hold off;
        
        drawnow
        
        % Store results
        tInt(k)  = endPnt;
        
        % Update interval endpoint
        if iMode==1
            lEnd = endPnt;
        else
            rEnd = endPnt;
        end
        
        % Prompt for input
        [~,~,b] = ginput(1);
        
        % Up arrow
        if b==30
            endPnt = endPnt + 0.002;
            
        % Down arrow
        elseif b==31
            endPnt = endPnt - 0.002;
            
        % Right arrow
        elseif b==29 
            endPnt = endPnt + 0.002;
            
        % Left arrow
        elseif b==28 
            endPnt = endPnt - 0.002;
            
        % esc
        elseif b==27
            disp('Quitting (data not saved)')
            return
            
        % Return
        elseif isempty(b)
            
            intSet = 1;

        end
    end
    
end

pause(1)

% Create the field 'tolSet', indicating that the tolerances have been set.
D.intSet = 1;

% Reshape intervals and save for output
D.tInt = reshape(tInt,[numTurns, 2]);

% disp(' ') disp(['Smoothing variable: ' num2str(i) ' of '
% num2str(length(aRev)) ' completed'])

%         % Save sequence stats save([sPath.data filesep 'Timing
%         data.mat'],'dT');

end


% close