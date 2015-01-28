function pColor = colorgradpercent(color0,color1,percentile,scalingFactor)
% COLORGRADPERCENT Returns a color vector by the percentile
%    pColor = COLORGRADPERCENT(color0,color1,percentile), where
%    color0 and color1 are RGB color vectors and percentile is a scalar
%    between 0 and 1, returns pColor, a color vector given by the
%    percentile of a gradient formed by color0 and color1. 
% 
% v1, Anthony Ho, 10/24/2014


    %% Reading arguments

    % Checking the number of arguments
    error(nargchk(3,4,nargin));
    
    if nargin<4
        scalingFactor = 10;
    end
    
    % Making sure percentile is between 0 and 1
    if percentile>1 
        percentile = 1;
    elseif percentile<0
        percentile = 0;
    end
    
    
    %% Computing the new color vector
    
    % Defining the transform function (logistic function in this case)
    transform = @(x) 1/(1+exp(-scalingFactor*(x-0.5)));
    
    % Return the color vector by percentile
    pColor = (color1-color0)*transform(percentile)+color0;
    

end