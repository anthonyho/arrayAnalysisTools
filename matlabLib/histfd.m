function [nelements,centers] = histfd(varargin)
% HISTFD  Histogram using the Freedman-Diaconis rule.
%    N = HISTFD(Y) bins the elements of Y into N bins calculated by the 
%    Freedman-Diaconis rule.
%    
%    N = HISTFD(Y,K), where K is a scalar, shows a normalized histogram if
%    K~=0. Default K=0. 
%
%    N = HISTFD(Y,K,LP,UP), where LP and UP are scalar between 0 and 100 
%    with UP>LP, shows a histogram ranging from LB to UB given by the 
%    percentile specified by LP and UP where LB = prctile(Y,LP) and 
%    UB = prctile(Y,UP).
%
%    N = HISTFD(Y,K,LP,UP,SF), where SF is a scalar, shows a histogram with
%    the bin width scaled by the scaling factor SF. 
%
%    [N,X] = HISTFD(...) also returns the position of the bin centers in X.
% 
%    HISTFD(...) without output arguments produces a histogram bar plot of
%    the results. 
%
%    HISTFD(AX,...) plots into AX instead of GCA.
% 
% v1, Anthony Ho, 10/22/2014


    %% Reading arguments
    
    % Parse possible Axes input
    error(nargchk(1,5,nargin,'struct'));
    [cax,args,nargs] = axescheck(varargin{:});
    
    % Reading data
    data = args{1};
    
    % Reading the rest of the arguments
    if nargs<2
        normalized = 0;
    else
        normalized = args{2};
    end
    if nargs<3
        lowerPercentile = 0;
    else
        lowerPercentile = args{3};
    end
    if nargs<4
        upperPercentile = 100;
    else
        upperPercentile = args{4};
    end
    if nargs<5
        scalingFactor = 1;
    else
        scalingFactor = args{5};
    end
    
    
    %% Defining constants
    
    percentOffset = 0.01;
    
    
    %% Processing data
    
    N = length(data);                                   % Number of data points
    Nprime = N*(upperPercentile-lowerPercentile)/100;   % Effective number of data points within the percentile range
    h = scalingFactor*2*iqr(data)*(Nprime^(-1/3));      % Calculate bin width
    upperBound = prctile(data,upperPercentile);         % Upper bound of histogram data
    lowerBound = prctile(data,lowerPercentile);         % Lower bound of histogram data
    bins = lowerBound:h:upperBound;                     % Defining bins

    [nelements,centers] = hist(data,bins);              % Calculating histogram data
    
    
    %% Plotting
    
    if nargout==0
        
        if isempty(cax)
            cax=gca;
        end
        
        if normalized
            bar(cax,centers,nelements/N,'hist');
        else
            bar(cax,centers,nelements,'hist');
        end
        
        offset = (upperBound-lowerBound)*percentOffset;  % Space at the ends of the histogram
        xlim(cax,[lowerBound-offset, upperBound+offset]);
    
    end

    makepretty;
    
    
end