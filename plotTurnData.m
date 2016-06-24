%% Path definitions & Load data

% % Alberto's MacMini (office)
if ~isempty(dir([filesep fullfile('Users','alberto','Documents','GitHub-SourceTree')]))
    path = '/Volumes/VisualPred/ZF_visuomotor/Turn_Data';
        
elseif ~isempty(dir([filesep fullfile('Users','mmchenry','Documents','Projects')]))
    % Matt's computer
    path = '/Users/mmchenry/Dropbox/Labbies/Alberto/predator';
else
    % Alberto's MacBook (laptop)
    path = '/Volumes/BackUp/ZF_visuomotor/Turn_Data';
end

% path = [path filesep 'Spring_2016/Turn_Data'];

% Load the pooled data 
load([path filesep '2016-05-19/2016-05-19-TurnData'])

% Cell array of variable names (used for axes labels)
names = {'\Delta Heading (deg)'; '\Delta Heading_{pre} (deg)';...
    '\Delta Gaze_{pre} (deg)'; '\Delta Bearing_{pre} (deg)';...	
    '\Delta \alpha_{pre} (deg)';...
    '\Delta \theta_E (deg)'; '\Delta \theta_{EG} (deg)';...
    '\Delta d(\theta_E)/dt (deg/s)';'{\theta_E}_{pre} (deg)';...
    '{\theta_EG}_{pre}  (deg)'; ' Pre d(\theta_E)/dt  (deg/s)'; ...	
    'Bearing_{pre} (deg)';'Bearing_{post} (deg)';'bearD1Pre';...
    '\delta_{pre}'; '\delta_{post}'; '\Delta \delta_{pre}'; ...
    'Distance_{pre} (cm)';'Dist_{post} (cm)';... 
    'Turn Duration (s)';'Inter Turn Duration (s)'};

% get a figure window and color order

color = get(gca,'colororder');

%% Plot all data 

% Get change in heading (during turn) in degrees
y = data(:,1)*180/pi;

% iteratively create plots

for k=1:length(names)

h = figure;
    
% Data to plot
if k<18
    % Angle data (convert to degrees)
    x = (data(:,k)*180/pi);
else
    % Dimensionless data
    x = data(:,k);
end

% indx = y > 50;
% y = y(~indx);
% x = x(~indx);

% Coefficients of linear fit
p = polyfit(x,y,1);

% y-values of fit
yfit = polyval(p,x);

% Compute the residual values as a vector of signed numbers:
yresid = y - yfit;

% Square the residuals and total them to obtain the residual sum of
% squares:
SSresid = sum(yresid.^2);

% Compute the total sum of squares of y by multiplying the variance of y by
% the number of observations minus 1:
SStotal = (length(y)-1) * var(y);

% Compute R^2 using the formula:
rsq = 1 - SSresid/SStotal;

% Correlation coefficient matrix
corrMatrix = corrcoef(x,y);

% Correlation coefficient between x & y
coefVal = corrMatrix(2);

% Plot change in heading vs. independent variable
plot(x,y,'o') % 'Color',color(5,:).*[0 0.9 1.05])
xlabel(names{k})
ylabel(names{1})

hold on

% Plot the line of isometry
line([min(ylim), max(ylim)], [min(ylim), max(ylim)],...
    'LineStyle','--','Color',[0.1, 0.1, 0.1])

% Plot linear regression
plot(x,yfit,'-k')

% Figure properties
% axis equal 
grid on

% Get x and y coordinates for placing text
xText = min(xlim) + diff(xlim)/20;
yText = max(ylim) - diff(ylim)/10;

% Show R^2 value
text(xText, yText,['R^2 = ' num2str(rsq)])

% Show correlation coefficient
text(xText,yText-diff(ylim)/20, ['Corr. Coeff. = ' num2str(coefVal)]);

hold off

% Set Fontsize
set(findobj(h, '-property', 'FontSize'),'FontSize',14)

% Save figure
% print(h,['HeadingChange-vs-Var ' num2str(k)],'-dtiff')

pause
% clear x y yfit 

close all

end
