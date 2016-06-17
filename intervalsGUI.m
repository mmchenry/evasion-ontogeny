function D = intervalsGUI(D,numTurns,hdAngle,mode)
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

switch mode
    
    % Interactively check accuracy of labeled turning events
    case 'turns'
        
        % Instructions
        disp('Verify turning events')
        disp('     Spacebar: keep event');
        disp('     delete: erase event');
        disp('     esc: quit');
        
        % Compute Angular velocity
        hdAngle_D1 = fnder(hdAngle);
        
        % Time values of all turning events from input
        tAV = D.tAV;
        
        for j = 1:D.numTurns
            
            % Indicator for verifying a turn (0 means turn not yet verified)
            eventSet = 0;
            
            % Current time point
            turnT = tAV(j);
        
            while ~eventSet
                
                % Plot heading data
                subplot(2,1,1)
                fnplt(hdAngle);
                hold on;
                
                % Show location of turning event on heading plot
                plot(turnT,fnval(hdAngle,turnT),'*r','MarkerSize',12)
                
                % Get current axes limits
                ax = gca;
                xText = ax.XLim(1) + ax.XLim(1)/5;
                yText = ax.YLim(2) - ax.YLim(2)/3;
                
                % Text description for interval
                if j > 1
                    intText = ['Prev. turn at t= [' num2str(tAV(j-1)) ']'];
                    text(xText,yText,intText, 'FontSize',14);
                end
                
                if j < D.numTurns
                    intText = ['Next turn at t= [' num2str(tAV(j+1)) ']'];
                    text(xText,yText/5,intText, 'FontSize',14);
                end
                
                grid on
                xlabel('Time (s)')
                ylabel('Heading Angle')
                hold off;
                
                % Plot rotational velocity
                subplot(2,1,2)
                fnplt(hdAngle_D1);
                hold on;
                
                % Show location of turning event on rotational velocity plot
                plot(turnT,fnval(hdAngle_D1,turnT),'*r','MarkerSize',12)
                
                grid on
                xlabel('Time (s)')
                ylabel('Heading angular velocity')
                hold off;
                
                drawnow
                
                % Store results
%                 turnT(j)  = endPnt;
                
                % Prompt for input
                [~,~,b] = ginput(1);
                
                % Delete
                if b==8
                    
                    % Delete current turning event (-1 is a placeholder)
                    D.tAV(j) = -1;
                    
                    % Update number of turns
                    numTurns = numTurns - 1;
                    
                    % Turning event has been set
                    eventSet = 1;
                    
                    % esc
                elseif b==27
                    disp('Quitting (data not saved)')
                    return
                    
                    % Spacebar
                elseif b==32
                    
                    % Turning event has been set
                    eventSet = 1;
                    
                end
            end
        end
        
        % Save updated turning event times and total number of turns
        D.tAV(D.tAV < 0) = [];
        D.numTurns = numTurns;
        
        
        
    % Interactively choose time intervals
    case 'intervals'
        
        % Instructions
        disp('Selecting period before and after a turn')
        disp('     right/up: move endpoint to the right');
        disp('     left/down: move endpoint to the left');
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
                    endPnt = endPnt + 0.001;
                    
                    % Down arrow
                elseif b==31
                    endPnt = endPnt - 0.001;
                    
                    % Right arrow
                elseif b==29
                    endPnt = endPnt + 0.001;
                    
                    % Left arrow
                elseif b==28
                    endPnt = endPnt - 0.001;
                    
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