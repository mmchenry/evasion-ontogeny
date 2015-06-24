function skel = remove_branches(skel,B,E,xHead,yHead,xTail,yTail)
% Removes branches from a skeleton


% Define center cross im, peripheral im, diagonal im
im_cntr = [0 1 0;1 0 1; 0 1 0];
im_peri = [1 1 1;1 0 1;1 1 1];
im_diag = [1 0 1; 0 0 0;1 0 1];


%% Prep images

% Coordinates of endpoints
[yEnd,xEnd] = find(E);

% Remove head and tail endpoints
E(yTail,xTail) = false;
E(yHead,xHead) = false;

% Make clubs out of head and tail
ptsX = [-1 0 1 1 0 -1 -1];
ptsY = [-1 -1 -1 0 1 1 1 0];

skel(ptsY+yTail,ptsX+xTail) = true;
skel(ptsY+yHead,ptsX+xHead) = true;

clear iVal dists

% Pad images
skel = pad_im(skel);
B = pad_im(B);
E = pad_im(E);

% New coordinates of endpoints
[yEnd,xEnd] = find(E);






% 
% %% Eliminate branches, start
% 
% % Step thru skeletal points, until all branches removed
% for i = 1:length(xEnd)
%     
%     % Start with endpoints
%     xLast = xEnd(i);
%     yLast = yEnd(i);
%     
%     % Black out endpoint
%     %skel_hist(yEnd(i),xEnd(i)) = false;
%     
%     %     xLast = max([xEnd 2]);
%     %     xLast = min([xLast size(skel_hist,2)-1]);
%     %
%     %     yLast = max([yLast 2]);
%     %     yLast = min([yLast size(skel_hist,1)-1]);
%     while true       
%         
%          % Stop if branchpoint reached
%         if B(yLast,xLast)==1
%             break
%         end
%         
%         % Black out current
%         skel(yLast,xLast) = false;
%         
%         % Define little roi
%         imLocal = skel((yLast-1):(yLast+1),(xLast-1):(xLast+1));
%         
%         % Local branchpoints
%         %B_local = B((yLast-1):(yLast+1),(xLast-1):(xLast+1));
%             
%         % If next point on vertical or horizontal
%         if max(max(imLocal & im_cntr))
%             
%             % Get next vertical or horizontal
%             [yNext,xNext] = find(imLocal & im_cntr,1,'first');
%             
%             % Otherwise, get next diagonal
%         else
%             % Find next pixel in local FOR
%             [yNext,xNext] = find(imLocal & im_diag,1,'first');
%         end
%         
%         % Check that we have a point
%         if isempty(xNext)
%             error('Oops');
%         end
%         
%         % Visual check
%         if 0
%             warning off
%             subplot(1,2,1)
%             imshow(skel);hold on
%             plot(xLast,yLast,'r+');hold off
%             subplot(1,2,2)
%             imshow(imLocal)
%             pause(0.01);
%         end
%         
%         % Displace for next iteration
%         xLast = xLast + xNext - 2;
%         yLast = yLast + yNext - 2;
%         
%     end
% end
% 
% % Unpad images
% skel = unpad_im(skel);



function im = pad_im(im)
% Pad images with single border of black pixels

% Top row
im = [false(1,size(im,2)); im];

% Bottom row
im = [im; false(1,size(im,2))];

% Left column
im = [false(size(im,1),1) im];

% Right column
im = [im false(size(im,1))];


function im = unpad_im(im)
% UnPad images with single border of black pixels

% Top row
im = im(2:end,:);

% Bottom row
im = im(1:end-1,:);

% Left column
im = im(:,2:end);

% Right column
im = im(:,1:end-1);   
    
