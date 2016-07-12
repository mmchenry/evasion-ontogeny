function [stats,slope,intercept] = reducedMajorAxis(X,Y,b_predicted,alpha,plotting,t_txt)
% Runs  type 2 regression stats
% b_predicted is the predicted slope of the line
% alpha is the p-value at which to accept/reject the null hypothesis (usually .05)
% plotting: 1 makes plots, 0 doesn't
% Code tested using examples from Sokal and Rolf (3rd ed., 1995) p.546, where
% x = [14;17;24;25;27;33;34;37;40;41;42];y=[61;37;65;69;54;93;87;89;100;90;97];

if nargin < 4, alpha   = .05;end
if nargin < 5, plotting = 1; end
if nargin < 6, t_txt = [];   end

if (size(X,2)>size(X,1)) | (size(Y,2)>size(Y,1))
    error('Data should be arranged in columns, not rows');
end

%In the case that you do not have information about the error of your samples, it is appropriate
%to use a Reduced Major Axis (i.e. geometric mean) regression.  This is the case if you don't have repeated
%measures (i.e. the x and y data have only one column).
if size(X,2)==1
    
    %First, find standard linear regression stats:
    x       = X - mean(X);
    y       = Y - mean(Y);
    n       = size(Y,1);
    b       = sum(x.*y)./sum(x.^2); %slope
    a       = mean(Y) - b * mean(X); %intercept
    
    Ypred   = b.*X + a; %predicted y values
    s_b     = ( (sum((Y-Ypred).^2) ./ (n-2)) ./ ...
        sum(x.^2) ).^.5; %standard error of regression
    s_Y     = ( (sum((Y-Ypred).^2) ./ (n-2)) ./ ...
        n ).^.5; %standard error of regression
    
    %calculate confidence intervals for linear regression:
    tValue  = tinv(1-alpha/2,n-2);
    t_s     = tValue.*s_b;
    t_sa    = tValue.*s_Y;
    L1      = b - t_s;
    L2      = b + t_s;
    
    %Now, find RMA stats:
    v       = ( sum(y.^2) / sum(x.^2) ) .^.5;
    %note: the stat v normally cannot find the sign of the slope,
    
    %so we use the least squares regression to find the sign:
    c       = corrcoef(X,Y);
    slope   = (c(2)/abs(c(2))) * v;
    a_v     = mean(Y) - slope * mean(X);
    s_v     = s_b;
    
    %and find confidence intervals for v:
    L1      = v - t_s;
    L2      = v + t_s;
    
    %and find confidence intervals for A_v:
    aL1     = a_v - t_sa;
    aL2     = a_v + t_sa;
    
    %prediction stats:
    a_pred  = mean(Y) - b_predicted * mean(X);
    
    
    intercept           = a_v;
    stats.slope         = slope;
    stats.lowerLimit_b  = L1;
    stats.upperLimit_b  = L2;
    stats.lowerLimit_a  = aL1;
    stats.upperLimit_a  = aL2;
    stats.alpha         = alpha;
    stats.r2            = rSquared(x,y,slope,a_v);
    
    
    if (b_predicted<L2) & (b_predicted>L1)
        stats.hypothesis    = 'accept null';
        title_str = [t_txt ' m = ' num2str(b_predicted,'% 10.2f\n') ...
               ', r^2= ' num2str(stats.r2,'% 10.2f\n')];
    else
        stats.hypothesis    = 'reject null';      
        title_str = [t_txt ' m = ' num2str(stats.slope,'% 10.2f\n') ...
               ', r^2= ' num2str(stats.r2,'% 10.2f\n')];
    end
    
    
    
    if plotting
        %figure;
        pointColor  = 'k';
        h = plot(X,Y,'o');hold on
        yL = ylim;xL = xlim;
        HL = max(abs([yL xL]));
        set(h,'MarkerFaceColor',pointColor)
        set(h,'MarkerEdgeColor',pointColor)
        set(h,'MarkerSize',3)
        h = plot([min(X) max(X)],slope*[min(X) max(X)]+a_v,'k');
        h = plot([-HL HL],b_predicted*[-HL HL]+a_pred,'--');
        %title([t_txt stats.hypothesis, ', r^2= ' num2str(stats.r2,'% 10.2f\n')])
        title(title_str)
        axis square
        axis([-HL HL -HL HL])
        grid on
        hold off
    end
    
    
end


function r2 = rSquared(x,y,slope,intercept)
%Gives r2 for a linear regression

yPred   = slope.*x + intercept;
yMean   = mean(y);

sPred       = sum( (yPred - yMean).^2 );
sTotal      = sPred + sum( (y - yPred).^2 );

r2 = sPred/sTotal;