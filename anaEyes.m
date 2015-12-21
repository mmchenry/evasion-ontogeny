function anaEyes(dPath,vPath)

showAna = 1;

% Tolerance for spline fits
tol.head = 1e2;
tol.rost = 0.5e1;

% Target number of pixels for eye area
eye_area = 100;


%% Load data

% Load batch data 'B'
load([dPath filesep 'blob data.mat'])

% Load midline data
load([dPath filesep 'Midline data.mat'])

% List of frame files
a = dir([vPath filesep '*.jpg']);

% Defaults for right eye
pEye.R.xCent = 23;
pEye.R.yCent = 25;
pEye.R.area = 100;
pEye.R.tVal = 1;

% Defaults fo left eye
%pEye.L.xCent = 23;
pEye.L.xCent = 21;
pEye.L.yCent = 9;
pEye.L.area = 100;
pEye.L.tVal = 1;


%% Fit splines to head and rostrum


% Disable warnings for spline fit
warning off

% Spline fit the data
sp.xRost = fnval(spaps(mid.t,mid.xRost,tol.rost),mid.t);
sp.yRost = fnval(spaps(mid.t,mid.yRost,tol.rost),mid.t);
sp.xHead = fnval(spaps(mid.t,mid.xHead,tol.head),mid.t);
sp.yHead = fnval(spaps(mid.t,mid.yHead,tol.head),mid.t);

warning on

% Visualize spline fits
if 0
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

% Cranial length
cran_len = mean(sqrt((sp.xHead-sp.xRost).^2 + (sp.yHead-sp.yRost).^2));


%% Acquire the eye mask

if isempty(dir([dPath filesep 'eye data.mat']))   
    
    startFrame = 1;
    
    % Interval for resetting reference frame
    frameResetIntrvl = 100;
    
    % Threshold adjusted to achieve the following area
    eye_area = 100;
    
    % Origin for head FOR in local system
    roi.xL          = -.3*cran_len;
    roi.yL          = -1.1*cran_len;
    roi.w           = 2.2*cran_len;
    roi.h           = 2.2*cran_len;
    roi.xCoordL = [roi.xL roi.xL+roi.w roi.xL+roi.w roi.xL roi.xL]';
    roi.yCoordL = [roi.yL roi.yL roi.yL+roi.h roi.yL+roi.h roi.yL]';
    roi.rEye    = cran_len*.35;
    
    % Dimensions of frame around head
    headIMdim = [2*cran_len 1.75*cran_len];
    
    % Initialize image registration
    [optimizer, metric]  = imregconfig('monomodal');
    optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;
    
     % Load image of first video frame
    im = imread([vPath filesep a(startFrame).name]);
    
    [imHead0,tform1] = giveHeadIM(im,sp,roi,startFrame);
    
    f = figure;
    imshow(imHead0,'InitialMagnification','fit')
    title('Right eye: select anterior, then posterior margin')
    [x,y,b] = ginput(2);
    hold on
    plot(x(1),y(1),'+r',x,y,'r-')
    
    % Initial right eye angle
    rAng0 = atan2(y(2)-y(1),x(2)-x(1));
    
    title('Left eye: select anterior, then posterior margin')
    [x,y,b] = ginput(2);
    plot(x(1),y(1),'+g',x,y,'g-')
    pause(2)
    close(f)
    
    % Initial left eye angle
    lAng0 = atan2(y(2)-y(1),x(2)-x(1));
    
    clear tform tform1 x y b
    
    bk_clr = round(1.5*min(double(im(:))));
    %bk_clr = 255;
    
    % Approximate centroid positions of eye
    pEye.R.xCent = 0.35*roi.w-roi.xL;
    pEye.R.yCent = abs(roi.yL)/2+roi.h*0.12;
    pEye.L.xCent = 0.35*roi.w-roi.xL;
    pEye.L.yCent = abs(roi.yL)/2+roi.h*.33;
    
    % Get eye coordinate data
    eL = give_eye(imHead0,pEye.L,eye_area,bk_clr,roi.rEye);
    eR = give_eye(imHead0,pEye.R,eye_area,bk_clr,roi.rEye);
    
    % Eye data for next iteration
    pEye0.R = eR;
    pEye0.L = eL;
    
    % Update pEye
    pEye = pEye0;
    
    clear eL eR
    
    iReset = 1;
    
    % Loop through frames
    for i = startFrame:length(sp.xRost)
        
        % Load image of video frame
        im = imread([vPath filesep a(i).name]);
        
        % Get smoothed positions for head coordinates
        head  = [sp.xHead(i) sp.yHead(i)];
        rost  = [sp.xRost(i) sp.yRost(i)];
        
        % Define coord transformation (rostrum as origin)
        tform = local_system(rost,head);
        
        % Origin in global FOR
        [roi.xCoordG,roi.yCoordG] = local_to_global(tform,roi.xCoordL,roi.yCoordL);
        
        % Head image
        [imHead,tform1] = giveHeadIM(im,sp,roi,i);
        

        % Transformation object to stablize head wrt imHead0
        tform2 = imregtform(imHead,imHead0,'rigid',optimizer, metric);
         
        % Stablize
        imStable = imwarp(imHead,tform2,'OutputView',imref2d(size(imHead0)));
        
        % Remove strips of black
        imStable(~im2bw(imStable,1-254/255)) = 255;        
        
        % Get eye coordinate data
        eL = give_eye(imStable,pEye.L,eye_area,bk_clr,roi.rEye);
        eR = give_eye(imStable,pEye.R,eye_area,bk_clr,roi.rEye);
        
        % Transformation object to stabalize head
        if i == startFrame
            eR.tform = imregtform(eR.im,pEye0.R.im,'rigid',optimizer, metric);
            eL.tform = imregtform(eL.im,pEye0.L.im,'rigid',optimizer, metric);
        else
            eR.tform = imregtform(eR.im,pEye0.R.im,'rigid',optimizer, metric,...
                            'InitialTransformation',pEye.R.tform);
            eL.tform = imregtform(eL.im,pEye0.L.im,'rigid',optimizer, metric,...
                            'InitialTransformation',pEye.L.tform);
        end
         
        % Stablize image
        imStableR = imwarp(pEye.R.im,eR.tform,'OutputView',...
                           imref2d(size(pEye0.R.im)));
        imStableL = imwarp(pEye.L.im,eL.tform,'OutputView',...
                           imref2d(size(pEye0.L.im)));
             
        % Inverse of first transformation             
        tmp = invert(tform1);

        % Head angle, from the inverse head angle and stablization
        eyes.hdAngle(i,1) = atan2(tmp.T(1,2),tmp.T(1,1)) + atan2(tform.T(1,2),tform2.T(1,1));
        eyes.rAngle(i,1) = atan2(eR.tform.T(1,2),eR.tform.T(1,1)) + rAng0;
        eyes.lAngle(i,1) = atan2(eL.tform.T(1,2),eL.tform.T(1,1)) + lAng0;
        
        % Eye data for next iteration
        pEye.R = eR;
        pEye.L = eL;   
        
        if iReset==frameResetIntrvl
            pEye0 = pEye;
            imHead0 = imStable;
            iReset = 1;
        else
            iReset = iReset + 1;
        end
 
        
        if 0
            % Border around ROI
            brdr = 20;
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

figure
subplot(2,1,1);
plot(eyes.t,unwrap(eyes.hdAngle-eyes.hdAngle(1))./pi*180,'-k');
hold on
plot(eyes.t,gazeR-gazeR(1),'-',...
     eyes.t,gazeL-gazeL(1),'-');
grid on;
xlabel('Time (s)')
ylabel('Head/Gaze angle (deg)')
%ylim([-135 135])

subplot(2,1,2);
plot(eyes.t,eyes.rAngle./pi*180,'-',eyes.t,eyes.lAngle./pi*180,'-');
legend('R','L');
grid on
xlabel('Time (s)')
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


function [imHead2,tform1] = giveHeadIM(im,sp,roi,cFrame)
    
% Get smoothed positions for head coordinate in first frame
head  = [sp.xHead(cFrame) sp.yHead(cFrame)];
rost  = [sp.xRost(cFrame) sp.yRost(cFrame)];

% Define coord transformation (rostrum as origin)
tform = local_system(rost,head);

% Origin in global FOR
[roi.xCoordG,roi.yCoordG] = local_to_global(tform,roi.xCoordL,roi.yCoordL);


%[xRost,yRost] = global_to_local(tform,rost(1),rost(2));

% Redefine coord transformation (originG as origin)
tform1 = local_system_special(rost,head,[roi.xCoordG(1) roi.yCoordG(1)]);

% Initial head image
imHead = imwarp(im,invert(tform1),'OutputView',imref2d(size(im)));

pixVal = imHead(1:round(roi.h),round(roi.w));

xCent = round(roi.w);
yCent = find(pixVal==min(pixVal),1,'first');

% Redefine coord transformation (originG as origin)
tform2 = local_system_special([-roi.xL -roi.yL],[xCent yCent],[0 0]);

% Head image with aligned horz axis
imHead2 = imwarp(imHead,invert(tform2),'OutputView',imref2d(size(im)));

% Remove strips of black
imHead2(~im2bw(imHead2,1-254/255)) = 255;   

% Crop both images
imHead = imHead([1:roi.w],[1:roi.h]);
imHead2 = imHead2([1:roi.w],[1:roi.h]);

%TODO: Calculate the total transformation

% itform1 = invert(tform1);
% itform2 = invert(tform2);
% 
% tform_tot = tform1;
% 
% tform_tot.T(1:2,1:2) = itform1.T(1:2,1:2)*itform2.T(1:2,1:2);
% tform_tot.T(3,1:2) = [itform1.T(3,1)+itform2.T(3,1) itform1.T(3,2)+itform2.T(3,2)];
% 
% imHead3 = imwarp(im,tform_tot,'OutputView',imref2d(size(im)));
% imHead3 = imHead3([1:roi.w],[1:roi.h]);
% 
% if 1
%     subplot(2,2,1)
%     imshow(imHead2)
%     subplot(2,2,2)
%     imshow(imHead3)
%     subplot(2,2,3:4)
%     imshowpair(imHead2,imHead3)
% end

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


function cEye = give_eye(im,pEye,eye_area,bk_clr,r)

% Adjust image contrast
%im = imadjust(im);

%TODO: Set initial centroid position.

% set initial threshold value
tVal = max([pEye.tVal-3 1]);

% set the stepsize for decreasing tVal
tStep = 1;

max_tVal = 255/2;

% adaptively find threshold value s.t. there are only two blobs
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
imMask = cast(zeros(size(im))+bk_clr,class(im));
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
imBW  = ~im2bw(im,tVal/255);

% Circle structuring element with radius=2
seEye = strel('disk',2);

% Perform a morphological close (dilation + erosion) operation on the image.
imBW = imerode(imBW,strel('disk',1));

%imBW = imopen(imBW,seEye);
% figure, imshow(imBW,'InitialMagnification','fit')
% title('Threshold + Open')

% find convex hulls of objects in imBW, with 4-connected neighborhood
imBW1 = bwconvhull(imBW,'objects',8);
% figure, imshow(CH_objects,'InitialMagnification','fit');
% title('Convex hulls')

% Identify blobs
LL    = bwlabel(imBW1);
props = regionprops(LL,'Centroid','Area');

% If there are blobs . . .
if ~isempty(props)
    
        % Select blob that includes previous centroid
        imBW2 = bwselect(imBW1,x,y,8);
        LL    = bwlabel(imBW2);
        props = regionprops(LL,'Centroid','Area');
        
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



