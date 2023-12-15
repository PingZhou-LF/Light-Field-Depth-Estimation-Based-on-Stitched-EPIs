function label=meanshsegm(im,hs,hr) 
% Input:
%     im   :[height,length,1],
%            image to be labeled
%     hs   :int,size of window for x-coordinate and y-coordinate
%     hr   :int,size of window for color domain
% Output:
%     label:[height,length,1],
%           label for different connected regions

%% Prepare the parameters
im                = double(im);
[height,length,~] = size(im);
h                 = [hs hs hr];% size of window for x-coordinate¡¢y-coordinate and color domain
clusterCenter                 = zeros( height, length, 3 );

for ix = 1:length
  for iy = 1:height
    %% generation of fw
    xmin = max(ix-hs,1);          % lower limit of window in the x direction
    xmax = min(ix+hs,length);     % upper limit of window in the x direction
    xlen = xmax-xmin+1;           % length of window in the x direction
    ymin = max(iy-hs,1);          % lower limit of window in the y direction
    ymax = min(iy+hs,height);     % upper limit of window in the y direction
    ylen = ymax-ymin+1;           % length of window in the y direction
    len  = xlen*ylen;             % length of window in the y direction
    iw   = (0:(len-1))';          % number of elements in the window
    fw   = [fix(iw/ylen+xmin) ... % x-coordinate of the element in the window
           mod(iw,ylen)+ymin  ... % y-coordinate of the element in the window
           reshape(im(ymin:ymax,xmin:xmax),len,1)]; % grayscale of the elements in the window
    %% mean shift
    center = double([ix iy im(iy,ix,:)]);% x-coordinate¡¢y-coordinate and grayscale of the current center point 
    
    if(sum(isnan(center))>0)
        label=ones(height,length);
        return
    end
    while true
      dis = (fw-repmat(center,len,1)) ./ repmat(h,len,1);
      dis = sum( dis.*dis, 2 );% distance between elements in the window and the current center point 
      index = dis<1.0;% find elements which are close to the current center point 
      temp = center;
      center = mean( fw(index,:), 1 );% find the new the new current point 
      if sqrt(sum((center-temp).^2))<1e-3*max(hs,hr)% condition of convergence 
          break; 
      end
    end
    clusterCenter(iy,ix,:) = center;    % z(:,:,1):convergence of x-coordinate
                            % z(:,:,2):convergence of y-coordinate
                            % z(:,:,3):convergence of z-coordinate

  end
end 

%% find connected regions
% set the pixel corresponds to the img to 1,between which are zero-pixels
% corresponding to the state of connection(0:disconnected
s = ones( 2*height+1, 2*length+1, 'int8' );
s(1:2:(2*height+1),:) = zeros( height+1, 2*length+1, 'int8' );
s(:,1:2:(2*length+1)) = zeros( 2*height+1, length+1, 'int8' );

% set the pixel satisfying the hs/hr constraints to 1(1:connected
% horizontal edges
s(2:2:2*height,3:2:(2*length-1)) = all(cat(3, ...             
  abs(clusterCenter(:,2:end,1:2)-clusterCenter(:,1:(end-1),1:2)) < hs, ...     % spatial distance between two adjacent elements < hs?
  abs(clusterCenter(:,2:end,3)-clusterCenter(:,1:(end-1),3)) < hr ),3);        % gray difference between two adjacent elements < hr?
% vertical edges
s(3:2:(2*height-1),2:2:2*length) = all(cat(3, ...         
  abs(clusterCenter(2:end,:,1:2)-clusterCenter(1:(end-1),:,1:2)) < hs, ...
  abs(clusterCenter(2:end,:,3:end)-clusterCenter(1:(end-1),:,3:end)) < hr ),3);

% find connected regions(4 connected domain
label = bwlabel( s, 4 );                                   
label = label( 2:2:2*height, 2:2:2*length ); % extract labeling                            


