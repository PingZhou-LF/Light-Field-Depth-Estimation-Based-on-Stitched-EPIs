function bw = adaptivethresh(im, fsize, t, filterType, thresholdMode)

    % Set up default parameter values as needed
    if nargin < 2
	fsize = fix(length(im)/20);
    end
        
    if nargin < 3
	t = 15;
    end
    
    if nargin < 4
	filterType = 'gaussian';
    end    
    
    if nargin < 5
	thresholdMode = 'relative';
    end
    
    % Apply Gaussian or median smoothing
    if strncmpi(filterType, 'gaussian', 3)
	g = fspecial('gaussian', 6*fsize, fsize);
	fim = filter2(g, im);
    elseif strncmpi(filterType, 'median', 3)
	fim = medfilt2(im, [fsize fsize], 'symmetric');	
    else
	error('Filtertype must be ''gaussian'' or ''median'' ');
    end
    
    % Finally apply the threshold
    if strncmpi(thresholdMode,'relative',3)
	bw = im > fim*(1-t/100);
    elseif  strncmpi(thresholdMode,'fixed',3)
	bw = im > fim-t;
    else
	error('mode must be ''relative'' or ''fixed'' ');
    end