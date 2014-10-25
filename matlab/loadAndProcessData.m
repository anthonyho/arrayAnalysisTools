function [data,avgData,medData] = loadAndProcessData(filepath)
% Load, filter, normalize, and compute the average/median of the time series
%
% Anthony Ho, 10/24/2014


    %% Constants
    
    threshold = 0.05;
    

    %% Import data
    
    disp(sprintf('   Processing %s',filepath));
    
    data = importdata(filepath);
    
    N0 = size(data,1);
    nTimePoint = size(data,2);
    avgData = zeros(nTimePoint,1);
    medData = zeros(nTimePoint,1);
    
    
    %% Normalization
    
    for j = 1:size(data,1)
        if data(j,1)>threshold
            amp = data(j,1);
            data(j,:) = data(j,:)/amp;
        end
    end

    
    %% Cleaning up data
    
    rowsToRemove = any(data(:,1)<threshold,2);
    data(rowsToRemove,:) = [];
    
    rowsToRemove = any(data==0,2);
    data(rowsToRemove,:) = [];

    N1 = size(data,1);
    
    %% Averaging
    
    for i=1:nTimePoint
        avgData(i) = mean(data(:,i));
        medData(i) = median(data(:,i));
    end
      
    
    %% Summary
    
    disp(sprintf('   Processed %d clusters.',N0));
    disp(sprintf('   Retained %d clusters (%d%%).\n',N1,round(N1/N0*100)));
    
    
end
