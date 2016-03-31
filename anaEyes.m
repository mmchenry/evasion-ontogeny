function anaEyes(dPath,vPath,startFrame,redoEyes)

% Indicator to visualize process
showAna = 0;

% Indicator for contrast adjustment of image
adjustOn = 0;

% Tolerance for spline fits
tol.head = 0.005e2;         % default: tol.head = 1e2
tol.rost = 0.250e1;         % default: tol.rost = 0.5e1

% Target number of pixels for eye area
eye_area = 100;

% Set startFrame
if nargin < 3
    startFrame = 1;
end


%% Load data

% Load batch data 'B'
load([dPath filesep 'blob data.mat'])

% Load midline data
load([dPath filesep 'Midline data.mat'])

% % List of frame files
% a = dir([vPath filesep '*.jpg']);

% Alternate to load filenames (use when 'ghost' files present) 
a = dir([vPath filesep 'exp*']);

% Defaults for right eye
pEye.R.xCent = 23;
pEye.R.yCent = 25;
pEye.R.area = 100;
pEye.R.tVal = 40;

% Defaults fo left eye
%pEye.L.xCent = 23;
pEye.L.xCent = 21;
pEye.L.yCent = 9;
pEye.L.area = 100;
pEye.L.tVal = 40;

%% Get Mean image

newMean = 1;

imMean = makeMeanImage(dPath,vPath,[],B,newMean);


%% Fit splines to head and rostrum


% Disable warnings for spline fit
warning off

% intial smoothing of data (use when midline data has a few errors)
% mid.xRost = smooth(mid.xRost,0.005,'rloess');
% mid.yRost = smooth(mid.yRost,0.005,'rloess');
mid.xHead = smooth(mid.xHead,0.020,'rloess');
mid.yHead = smooth(mid.yHead,0.020,'rloess');

% Spline fit the data
sp.xRost = fnval(spaps(mid.t,mid.xRost,tol.rost),mid.t);
sp.yRost = fnval(spaps(mid.t,mid.yRost,tol.rost),mid.t);
sp.xHead = fnval(spaps(mid.t,mid.xHead,tol.head),mid.t);
sp.yHead = fnval(spaps(mid.t,mid.yHead,tol.head),mid.t);

warning on

% Visualize spline fits
if 1
    figure
    subplot(2,1,1)
    plot(mid.xRost,mid.yRost,'.',sp.xRost,sp.yRost,'-')
    axis square
    axis equal
    ylabel('Rost y')
    xlabel('Rost x')
    
    subplot(2,1,2)
    plot(mid.xHead,mid.yHead,'.',sp.xHead,sp.yHead,'-')
    axis square
    axis equal
    ylabel('Head y')
    xlabel('Head x')   
end
pause
% Cranial length
cran_len = mean(sqrt((sp.xHead-sp.xRost).^2 + (sp.yHead-sp.yRost).^2));


%% Acquire the eye mask

if isempty(dir([dPath filesep 'eye data.mat'])) || redoEyes
    
    % Interval for resetting reference frame
    frameResetIntrvl = 100;
    
    % Threshold adjusted to achieve the following area
    eye_area = 100;
    
    % Origin for head FOR in local system
    roi.xL          = -0.35*cran_len;       % default = -0.3
    roi.yL          = -1.1*cran_len;        % default = -1.1
    
    % width and height of bounding box around head
    roi.w           = 2.2*cran_len;
    roi.h           = 2.2*cran_len;
    
    % coordinates for bounding box, local FOR
    roi.xCoordL = [roi.xL roi.xL+roi.w roi.xL+roi.w roi.xL roi.xL]';
    roi.yCoordL = [roi.yL roi.yL roi.yL+roi.h roi.yL+roi.h roi.yL]';
    
    % radius for circle around eyes
    roi.rEye    = cran_len*0.36;     % default = 0.35
    
    % Dimensions of frame around head
    headIMdim = [2*cran_len 1.75*cran_len];
    
    % Initialize image registration
    [optimizer, metric]  = imregconfig('monomodal');
    optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;
    
     % Load image of first video frame
    im = imread([vPath filesep a(startFrame).name]);
    
    % cropped head image, aligned with horizontal
    [imHead0,~] = giveHeadIM(im,sp,roi,startFrame,0);
    
    f = figure;
    imshow(imHead0,'InitialMagnification','fit')
    title('Right eye: select anterior, then posterior margin')
    [x,y,~] = ginput(2);
    hold on
    plot(x(1),y(1),'+r',x,y,'r-')
    
    % Initial right eye angle (relative to body axis)
    % switched sign on y-coord 
    rAng0 = atan2(-(y(2)-y(1)),x(2)-x(1));
    
    title('Left eye: select anterior, then posterior margin')
    [x,y,~] = ginput(2);
    plot(x(1),y(1),'+g',x,y,'g-')
    pause(1)
    
    % Initial left eye angle (relative to body axis)
    % switched sign on y-coord
    lAng0 = atan2(-(y(2)-y(1)),x(2)-x(1));
    
    % Get initial Right eye centroid position
    title('Select center of Right eye')
    [x,y,~] = ginput(1);
    plot(x(1),y(1),'wo')
    
    % save centroid Right eye position
    pEye.R.xCent = x(1);
    pEye.R.yCent = y(1);
    
    % Get initial Left eye centroid position
    title('Select center of Left eye')
    [x,y,~] = ginput(1);
    plot(x(1),y(1),'yo')
    pause(1)
    
    % save centroid Right eye position 
    pEye.L.xCent = x(1);
    pEye.L.yCent = y(1);
    
    % close figure
    close(f)
    
    clear tform tform1 x y 
    
    bk_clr = round(1.5*min(double(im(:))));
    %bk_clr = 255;
    
    % Initial threshold value, based on first frame
    Reye_tVal   = imHead0(round(pEye.R.yCent),round(pEye.R.xCent));
    Leye_tVal   = imHead0(round(pEye.L.yCent),round(pEye.L.xCent));
    tVal0       = mean([Reye_tVal,Leye_tVal]);
    
    % estimate eye area from cropped image
    [eyeArea,tVal] = giveArea(imHead0,pEye,tVal0);
    
    % Get eye coordinate data
    eL = give_eye(imHead0,pEye.L,eyeArea.Leye,bk_clr,roi.rEye);
    eR = give_eye(imHead0,pEye.R,eyeArea.Reye,bk_clr,roi.rEye);
    
    % Eye data for next iteration
    pEye0.R     = eR;
    pEye0.L     = eL;
    pEye0.tVal  = tVal;
    
    % Update pEye
    pEye = pEye0;
    
    clear eL eR
    
    iReset = 1;
    
    % Loop through frames
    for i = startFrame:length(sp.xRost)
        
        % Load image of video frame
        im = imread([vPath filesep a(i).name]);
        
        % if newMean is ON...
        if newMean
            
            % ...Subtract background image (use meanImage2)
            im = imcomplement(imsubtract(imMean,im));
        end
        
        % Get smoothed positions for head coordinates
        head  = [sp.xHead(i) sp.yHead(i)];
        rost  = [sp.xRost(i) sp.yRost(i)];
        
        % Define coord transformation (rostrum as origin)
        tform = local_system(rost,head);
        
        % Origin in global FOR
        [roi.xCoordG,roi.yCoordG] = local_to_global(tform,roi.xCoordL,roi.yCoordL);
        
        % Head image (cropped and aligned to horizontal)
        [imHead,tform1,anglCor] = giveHeadIM(im,sp,roi,i,adjustOn,pEye.tVal);
        
        % NOTE: tform1 does not give the total rotation angle to obtain
        % 'imHead' from 'im'

        % Transformation object to stablize head wrt imHead0
        tform2 = imregtform(imHead,imHead0,'rigid',optimizer,metric);
         
        % Stablize head image
        imStable = imwarp(imHead,tform2,'OutputView',imref2d(size(imHead0)));
        
        % Remove strips of black
        imStable(~im2bw(imStable,1-254/255)) = 255;        
        
        % Get eye coordinate data
        eL = give_eye(imStable,pEye.L,eyeArea.Leye,bk_clr,roi.rEye);
        eR = give_eye(imStable,pEye.R,eyeArea.Reye,bk_clr,roi.rEye);
        
        % Transformation object to stabilize eyes
        if i == startFrame
            
            eR.tform = imregtform(eR.im,pEye0.R.im,'rigid',optimizer,metric);
            eL.tform = imregtform(eL.im,pEye0.L.im,'rigid',optimizer,metric);
            
        else
            
            % Attempt to stabilize Right eye
            try
                eR.tform = imregtform(eR.im,pEye0.R.im,'rigid',...
                    optimizer,metric,'InitialTransformation',pEye.R.tform);
                
           % If error above . . .
            catch
                % Use data from previous iteration
                eR.tform = pEye.R.tform;
            end
             
            % Attempt to stabilize Left eye
            try
                eL.tform = imregtform(eL.im,pEye0.L.im,'rigid',...
                    optimizer,metric,'InitialTransformation',pEye.L.tform);
                
           % If error above . . .
            catch
                % Use data from previous iteration
                eL.tform = pEye.L.tform; 
            end
            
        end
         
        % Stablize eye images
        if showAna
            imStableR = imwarp(pEye.R.im,eR.tform,'OutputView',...
                imref2d(size(pEye0.R.im)));
            imStableL = imwarp(pEye.L.im,eL.tform,'OutputView',...
                imref2d(size(pEye0.L.im)));
        end
             
        % Inverse of first transformation             
%         tmp = invert(tform1);

        % Head angle from first transformation (imHead)
        eyes.angl1(i,1) = atan2(tform1.T(1,2),tform1.T(1,1)) + anglCor;
        
        % Head angle correction from image registration (imStable)
        eyes.angl2(i,1) = atan2(tform2.T(1,2),tform2.T(1,1));
        
        % Head angle in world coordinates
        if sign(eyes.angl1)>=0 
            % if heading is between 0 and 180 ...
            eyes.hdAngle(i,1) = pi - (eyes.angl1(i,1) - eyes.angl2(i,1));
        else
            % ... otherwise
            eyes.hdAngle(i,1) = - pi - (eyes.angl1(i,1) - eyes.angl2(i,1));
        end
                        
        % Eye angles                
        eyes.rAngle(i,1) = atan2(eR.tform.T(1,2),eR.tform.T(1,1)) + rAng0;
        eyes.lAngle(i,1) = atan2(eL.tform.T(1,2),eL.tform.T(1,1)) + lAng0;
        
        % Eye data for next iteration
        pEye.R = eR;
        pEye.L = eL;   
        
        if iReset==frameResetIntrvl
            % reset eye data
            pEye0 = pEye;
            
            % reset reference head image
            imHead0 = imStable;
            
            % reset initial eye angles
            rAng0 = eyes.rAngle(i,1);
            lAng0 = eyes.lAngle(i,1);

            % reset the reset counter
            iReset = 1;
        else
            iReset = iReset + 1;
        end
 
        
        if showAna
            % Border around ROI
            brdr = 20;
            figure(1);
            subplot(1,2,1)
            imshow(im,'InitialMagnification','fit')
            hold on
            plot(roi.xCoordG,roi.yCoordG,'r-')
            %plot(tform3.T(3,1),tform3.T(3,2),'r+')
            hold off
            axis([min(roi.xCoordG)-brdr max(roi.xCoordG)+brdr  ...
                min(roi.yCoordG)-brdr max(roi.yCoordG)+brdr]);

            subplot(1,2,2)
            imshow(imStable,'InitialMagnification','fit')
            hold on
            plot(pEye.R.x,pEye.R.y,'r',pEye.L.x,pEye.L.y,'g')
            set(gcf,'Color',.5.*[1 1 1])
            pause(.01)

        end
        
        if showAna
            figure(2);
            subplot(3,2,1);imshow(imStableR,'InitialMagnification','fit')
            subplot(3,2,3);imshow(imStable,'InitialMagnification','fit')
            hold on; plot(pEye.R.xCirc,pEye.R.yCirc,'r-',...
                          pEye.L.xCirc,pEye.L.yCirc,'g-'); hold off
            subplot(3,2,5);imshow(imStableL,'InitialMagnification','fit')
            subplot(3,2,2);imshow(pEye0.R.im,'InitialMagnification','fit')
            subplot(3,2,4);imshow(imHead0,'InitialMagnification','fit')
            hold on; plot(pEye0.R.xCirc,pEye0.R.yCirc,'r-',...
                          pEye0.L.xCirc,pEye0.L.yCirc,'g-'); hold off
            subplot(3,2,6);imshow(pEye0.L.im,'InitialMagnification','fit')
            set(gcf,'Color',.3.*[1 1 1])
            
            pause(0.01);
        end
        clear tform tform2 tform3
     
        % Update status
        disp(['        Eye masks done for ' num2str(i) ' of ' ...
              num2str(length(sp.xRost))])
        
        % Visualize in local FOR
        if 0
            subplot(2,1,1);
            imshow(imRight,'InitialMagnification','fit')
            hold on
            plot(eyes(i).R.x,eyes(i).R.y,'r-')
            hold off
            title(['Frame ' num2str(i)])
            subplot(2,1,2);
            imshow(imLeft,'InitialMagnification','fit')
            hold on
            plot(eyes(i).L.x,eyes(i).L.y,'r-')
            hold off
            pause(0.1)
        end
        
        % Visualize in global FOR
        if 0
            imshow(im,'InitialMagnification','fit')
            hold on           
            plot(eyes(i).R.x,eyes(i).R.y,'g-',eyes(i).L.x,eyes(i).L.y,'y-')           
            hold off           
            pause(0.1)
        end
        
        clear eR eL
    end
    
    eyes.t = mid.t;
    
    % Save eye data
    save([dPath filesep 'eye data.mat'],'eyes')
    
else
    disp('            Loading eye data')
    load([dPath filesep 'eye data.mat'])
end


%% Visualize the results

gazeR = (unwrap(eyes.hdAngle)+eyes.rAngle+pi/2)./pi*180;
gazeL = (unwrap(eyes.hdAngle)+eyes.lAngle-pi/2)./pi*180;

totFrames = length(gazeR);

figure
subplot(2,1,1);
% plot(eyes.t,unwrap(eyes.hdAngle-eyes.hdAngle(1))./pi*180,'-k');
plot(eyes.t(1:totFrames)*250,unwrap(eyes.hdAngle)./pi*180,'-k');
hold on
% plot(eyes.t,gazeR-gazeR(1),'-',...
%      eyes.t,gazeL-gazeL(1),'-');
 plot(eyes.t(1:totFrames)*250,gazeR,'-',...
     eyes.t(1:totFrames)*250,gazeL,'-');
grid on;
% xlabel('Time (s)')
xlabel('Frame number')
ylabel('Head/Gaze angle (deg)')

subplot(2,1,2);
plot(eyes.t(1:totFrames)*250,eyes.rAngle./pi*180,'-')
hold on
plot(eyes.t(1:totFrames)*250,eyes.lAngle./pi*180,'-')
legend('R','L');
grid on
% xlabel('Time (s)')
xlabel('Frame number')
ylabel('Eye angle (deg)')



% Show eye traces
if 0
    for i = 1:length(eyes)
        
        subplot(1,2,1)
        plot(eyes(i).L.x,eyes(i).L.y,'-')
        hold on
        axis square
        axis equal
        
        subplot(1,2,2)
        plot(eyes(i).R.x,eyes(i).R.y,'-')
        hold on
        axis square
        axis equal
        
        %pause(.1)
    end
end


% Animate the data
if 0  
    figure;
    for i = 1:length(a)
        
        xRost = fnval(sp.xRost,mid.t(i));
        yRost = fnval(sp.yRost,mid.t(i));
        xHead = fnval(sp.xHead,mid.t(i));
        yHead = fnval(sp.yHead,mid.t(i));
        
        im = imread([vPath filesep a(i).name]);
        imshow(im,'InitialMagnification','fit')
        hold on
        plot(xRost,yRost,'+',xHead,yHead,'o')
        hold off
        title([num2str(i) ' of ' num2str(length(a))])
        pause(0.05)
        
    end
end


function [imHead2,tform1,anglCor] = giveHeadIM(im,sp,roi,cFrame,levels,tVal)
% INPUTS: 
%       - im  = original video frame 
%       - sp  = spline fitted head & rostrum data
%       - roi = region of interest
%       - cFrame = current frame
%       - levels = indicator for increasing image contrast
%       - tVal = threshold value
%
% OUTPUT: 
%       - imHead2 = image of cropped and aligned head 
%       - tform1  = transformation matrix to obtain image 

% Set tVal if not given as input
if nargin < 6
    tVal = 35;
end
    
% Get smoothed positions for head coordinate in current frame
head  = [sp.xHead(cFrame) sp.yHead(cFrame)];
rost  = [sp.xRost(cFrame) sp.yRost(cFrame)];

% Define coord transformation (rostrum as origin)
tform = local_system(rost,head);

% Bounding box around head in global FOR 
[roi.xCoordG,roi.yCoordG] = local_to_global(tform,roi.xCoordL,roi.yCoordL);

% Redefine coord transformation (originG as origin)
tform1 = local_system_special(rost,head,[roi.xCoordG(1) roi.yCoordG(1)]);

% Initial head image 
% (rotated and translated, rostrum near top-left, not aligned w/ horiz)
imHead = imwarp(im,invert(tform1),'OutputView',imref2d(size(im)));

pixVal = imHead(1:round(roi.h),round(roi.w));

% Get center point along fish body, posterior to head point
xCent = round(roi.w);
yCent = find(pixVal==min(pixVal),1,'first');

% rostrum point
rostP = [-roi.xL -roi.yL];

% Redefine coord transformation (originG as origin)
% This transformation more closely aligns the fish heading w/ horizontal
tform2 = local_system_special(rostP,[xCent yCent],[0 0]);

% Head image with aligned horz axis
imHead2 = imwarp(imHead,invert(tform2),'OutputView',imref2d(size(im)));

% Remove strips of black
imHead2(~im2bw(imHead2,1-254/255)) = 255;   

% Crop both images
imHead = imHead([1:roi.w],[1:roi.h]);
imHead2 = imHead2([1:roi.w],[1:roi.h]);

% Adjust cropped image contrast
if levels
    imHead2 = imadjust(imHead2,[tVal/255;1],[2/255;1]);
end

% compute correction angle, need to add this to tform1 for total angle
anglCor = atan2(tform2.T(1,2),tform2.T(1,1));

if 0
    figure
    subplot(1,2,1)
    imshow(imHead,'InitialMagnification','fit')
    hold on
    plot(-roi.xL,-roi.yL,'r+',[-roi.xL xCent],[-roi.yL yCent],'r-')
    hold off
    
    subplot(1,2,2)
    imshow(imHead2,'InitialMagnification','fit')
end



function tform3 = sum_tforms(tform1,tform2)
% Sums two coordinate transformations

% Check dimensions
if tform1.Dimensionality~=2 || tform2.Dimensionality~=2
    error('Code only handles 2D transformations')
end

tform1 = invert(tform1);

% Initialize tform3
tform3 = tform2;

% Angles
ang1 = atan2(tform1.T(1,2),tform1.T(1,1));
ang2 = atan2(tform2.T(1,2),tform2.T(1,1));
ang_tot = ang1 + ang2;

% Translations
trans = [tform1.T(3,1)+tform2.T(3,1) tform1.T(3,2)+tform2.T(3,2)];

%Transformation matrix
tform3.T =  [cos(ang_tot) sin(ang_tot) 0; ...
            -sin(ang_tot) cos(ang_tot) 0; ...
            trans(1) trans(2) 1];

%tform3 = invert(tform3);        
%inv_tform = invert(tform2);

%tform3.T = tform1.T*inv_tform.T;

%tform3.T(3,1:2) = [tform1.T(3,1)+tform2.T(3,1) tform1.T(3,2)+tform2.T(3,2)];

function imHead = headImage(im,tform,dims)

% % Angle of head
% imAng = atan2(head(2)-rost(2),head(1)-rost(1));
%     
% % Transformation matrix
% A = [cos(imAng) sin(imAng) 0; ...
%     -sin(imAng) cos(imAng) 0; ...
%     originG(1) originG(2) 1];
% 
% % Formatted
% tform = affine2d(A);
  
% Head image
imHead = imwarp(im,invert(tform),'OutputView',imref2d(size(im)));

% Crop
imHead = imHead([1:dims(1)],[1:dims(2)]);

function [eyeArea,tVal] = giveArea(imHead,pEye,tVal0)
% INPUTS: 
%   - imHead    = cropped head image
%   - pEye      = previous centroid position of eye
%   - r         = radius for circle around eye
%
% OUTPUTS:
%   - eyeArea   = target eye area 

% Indicator for while loop
proceed = 1;

% Set initial threshold value
% tVal1 = mean2(imHead);
tVal = tVal0;

% Maximum threshold value
max_tVal = 90;

while proceed && (tVal < max_tVal)
    
    % threshold image (eye blobs must be white!)
    imBW1  = ~im2bw(imHead,tVal/255);
    
    % Select blob that includes right eye position
    imReye = bwselect(imBW1,pEye.R.xCent,pEye.R.yCent,8);
    
    % Select blob that includes left eye position
    imLeye = bwselect(imBW1,pEye.L.xCent,pEye.L.yCent,8);
    
    % Measure region properties of right eye
    ReyeBlob = regionprops(imReye,'Area','Centroid','Orientation');
    
    % Measure region properties of light eye
    LeyeBlob = regionprops(imLeye,'Area','Centroid','Orientation');
    
%     % keep only the larger connected components, i.e., eye blobs
%     eyeBlobs = eyeBlobs([eyeBlobs.Area] > 10);
    
    % Check that there are only 2 blobs
    if (length(ReyeBlob)+length(LeyeBlob))==2
        
        % if so ... compute the ratio of the blob areas ...
        sizeRatio = min(ReyeBlob.Area,LeyeBlob.Area)/...
            max(ReyeBlob.Area,LeyeBlob.Area);
        
        % ... and check that they are of similar size
        if sizeRatio < 0.5
            % increase threshold value
            tVal = tVal + 1;
        else
            % exit while loop
            proceed = 0;
        end
        
        % otherwise ...
    else
        tVal = tVal + 1;
    end
end

% Store area of right eye
eyeArea.Reye = ReyeBlob.Area;

% Store area of left eye
eyeArea.Leye = LeyeBlob.Area;



function cEye = give_eye(im,pEye,eye_area,bk_clr,r)
% INPUTS: 
%   - im        = cropped head image
%   - pEye      = previous centroid position of eye
%   - eye_area  = target area to achieve (cosider an adaptive area)
%   - bk_clr    = average value of background?
%   - r         = radius for circle around eye
%
% OUTPUTS:
%   - cEye      = coordinates for eye blob & eye mask image

%TODO: Set initial centroid position.

% set initial threshold value
tVal = max([pEye.tVal-3 1]);

% set the stepsize for decreasing tVal
tStep = 1;

max_tVal = 255/2;

% adaptively find threshold value to reach target eye area
while tVal < max_tVal
    
    % Return blob
    [imBW1,cArea,cX,cY] = give_blob(im,tVal,pEye.xCent,pEye.yCent);
    
    % If exceeding target area . . .
    if cArea > eye_area
        % Step back, if area too large
        if cArea > eye_area*1.2
            [imBW1,cArea,cX,cY] = give_blob(im,tVal-1,pEye.xCent,pEye.yCent);         
        end
        
        % Quit loop
        break
        
    % Otherwise, advance threshold
    else
        tVal = tVal + tStep;
    end
    
end

% Get peripheral shapes
bb = bwboundaries(imBW1,'noholes');

if (cArea < .75*eye_area) || (tVal >= max_tVal)
    %cEye.area = pEye.area;
    cEye.x = nan;
    cEye.y = nan;
    cEye.xCent = pEye.xCent;
    cEye.yCent = pEye.yCent;
    cEye.tVal = pEye.tVal;
else       
    
    % Extract coordinates
    cEye.x = bb{1}(:,2);
    cEye.y = bb{1}(:,1);
    cEye.xCent = cX;
    cEye.yCent = cY;
    %cEye.area = cArea;
    cEye.tVal = tVal;
end

if 0
    imshow(im,'InitialMagnification','fit')
    hold on
    plot(cEye.x,cEye.y,'r-')
    hold off
end

% Mask image
%bw = roipoly(im,cEye.x,cEye.y);
[xCirc,yCirc] = makeCircle(cEye.xCent,cEye.yCent,r);
bw = roipoly(im,xCirc,yCirc);

% % Dilate the mask
% se = strel('disk',2);
% bw = imdilate(bw,se);

% Get image of just the eye
imMask = bk_clr * ones(size(im),'uint8');
imMask(bw) = im(bw);

% Cropping rectangle
rect = [round(cEye.xCent-r)-2 round(cEye.yCent-r)-2 2*r+4 2*r+4];

% Crop
imMask = imcrop(imMask,rect);

% Store
cEye.im = imMask;
cEye.xCirc = xCirc;
cEye.yCirc = yCirc;


function [imBW2,cArea,x,y] = give_blob(im,tVal,x,y)

% Threshold image and include only ROI
imBW1  = ~im2bw(im,tVal/255);

% Circle structuring element with radius=1
% seEye = strel('disk',1);

% Perform a morphological open (erosion + dilation) operation on the image.
% This gets rid of small objects in foreground 
% imBW1 = imopen(imBW,seEye);

% Identify blobs
props = regionprops(imBW1,'Centroid','Area');

% If there are blobs . . .
if ~isempty(props)
    
        % Select blob that includes previous centroid
        imBW2 = bwselect(imBW1,x,y,8);
        
        % measure its properties
        props = regionprops(imBW2,'Centroid','Area');
        
        if isempty(props)
            cArea = 0;
            
        elseif length(props)>1
            error('More than one blob!')
            
        else
            cArea = props(1).Area;
            x = props(1).Centroid(1);
            y = props(1).Centroid(2);
        end
    
    imBW2 = imdilate(imBW2,strel('disk',1));

% Otherwise, no area . . .
else
    cArea = 0;
    imBW2 = imBW1;
    x = nan;
    y = nan;
end

function tform = local_system(origin,xPoint)

% Check dimensions
if size(origin,1)~=1 || size(origin,2)~=2 || size(xPoint,1)~=1 || size(xPoint,2)~=2 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = xPoint(1) - origin(1);
xaxis(1,2) = xPoint(2) - origin(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% % Format used for angles:
% imAng = atan2(xaxis(2),xaxis(1));
%  A = [cos(imAng) sin(imAng) 0; ...
%               -sin(imAng) cos(imAng) 0; ...
%               origin(1) origin(2) 1];

% Format for matlab
tform = affine2d(R);

function tform = local_system_special(base,xPoint,origin)
% Defines a local system where the xaxis is defined by base and xPoint, 
% separate from the origin

% Check dimensions
if size(base,1)~=1 || size(base,2)~=2 || size(xPoint,1)~=1 || size(xPoint,2)~=2 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = xPoint(1) - base(1);
xaxis(1,2) = xPoint(2) - base(2);
xaxis(1,3) = 0;

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% Format used for angles:
% A = [cos(imAng) sin(imAng) 0; ...
%              -sin(imAng) cos(imAng) 0; ...
%              xPoint(1) xPoint(2)+0*round(size(im,1)/2) 1];

% Format for matlab
tform = affine2d(R);


function [xT,yT] = global_to_local(tform,x,y)
% Assumes column vectors for coordinates

% Check dimensions
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end

pts = [x y];

% Translate
pts(:,1) = pts(:,1) - tform.T(3,1);
pts(:,2) = pts(:,2) - tform.T(3,2);

% Rotate points
ptsT = [tform.T(1:2,1:2) * pts']';

% Extract columns of points
xT = ptsT(:,1);
yT = ptsT(:,2);


function [xT,yT] = local_to_global(tform,x,y)
% Assumes columns vectors for coordinates

% Check dimensions
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end

% Loop thru columns of coordinates
for i = 1:size(x,2)
    
    pts = [x(:,i) y(:,i)];
    
    % Rotate points
    ptsT = (tform.T(1:2,1:2) \ pts')';
    
    % Translate global coordinates wrt origin
    ptsT(:,1) = ptsT(:,1) + tform.T(3,1);
    ptsT(:,2) = ptsT(:,2) + tform.T(3,2);
    
    % Extract columns of points
    xT(:,i) = ptsT(:,1);
    yT(:,i) = ptsT(:,2);
    
    clear ptsT pts
end

function [x,y] = makeCircle(xC,yC,r)
theta = 0:pi/500:2*pi;
x = r.*cos(theta) + xC;
y = r.*sin(theta) + yC;



