function [nelements,centers] = perfhist(varargin)
% PERFHIST  Histogram with user-specific range and binsize
%    N = PERFHIST(Y,X), where X is a vector, returns the distribution of Y
%    among bins with centers specified by X. The first bin includes
%    data between -inf and the first center and the last bin
%    includes data between the last bin and inf.
%    
%    N = PERFHIST(Y,X,K), where K is a scalar, shows a normalized histogram
%    if K~=0. Default K=0. 
%
%    [N,X] = PERFHIST(...) also returns the position of the bin centers in 
%    X.
% 
%    PERFHIST(...) without output arguments produces a histogram bar plot
%    of the results. 
%
%    PERFHIST(AX,...) plots into AX instead of GCA.
%
% v1, Anthony Ho, 10/22/2014


    %% Reading arguments
    
    % Parse possible Axes input
    error(nargchk(2,4,nargin,'struct'));
    [cax,args,nargs] = axescheck(varargin{:});
    
    % Reading arguments
    data = args{1};
    bins = args{2};
    
    if nargs<3
        normalized = 0;
    else
        normalized = args{3};
    end
    
    
    %% Defining constants
    
    percentOffset = 0.01;
    
    
    %% Processing data
    
    N = length(data);                           % Number of data points
    
    [nelements,centers] = hist(data,bins);      % Calculating histogram data
    
    
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
        
        offset = (max(bins)-min(bins))*percentOffset;  % Space at the ends of the histogram
        xlim(cax,[min(bins)-offset, max(bins)+offset]);
    
    end

    makepretty;
    
    
end