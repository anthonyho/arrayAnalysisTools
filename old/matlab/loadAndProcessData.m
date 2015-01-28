function [nData,rData,avgData,medData] = loadAndProcessData(filepath)
% Load, filter, normalize, and compute the average/median of the time series
%
% Anthony Ho, 10/24/2014


    %% Constants
    
    thresholdPercentile = 5;
    

    %% Import data
    
    disp(sprintf('   Processing %s',filepath));
    
    rData = importdata(filepath);
    nData = rData;
    
    N0 = size(rData,1);
    nTimePoint = size(rData,2);
    avgData = zeros(nTimePoint,1);
    medData = zeros(nTimePoint,1);
    
    if N0 == 0
        return
    end
    
    threshold = prctile(rData(:,1),thresholdPercentile);
    %threshold = 0.05;
    
    
    %% Normalization
    
    for j = 1:N0
        
        if rData(j,1)>threshold
            amp = rData(j,1);
            nData(j,:) = rData(j,:)/amp;
        end
        
    end

    
    %% Cleaning up data
    
    rowsToRemove = any(rData(:,1)<threshold,2);
    rData(rowsToRemove,:) = [];
    nData(rowsToRemove,:) = [];
    
    rowsToRemove = any(rData==0,2);
    rData(rowsToRemove,:) = [];
    nData(rowsToRemove,:) = [];

    N1 = size(rData,1);
    
    %% Averaging
    
    for i = 1:nTimePoint
        
        avgData(i) = mean(nData(:,i));
        medData(i) = median(nData(:,i));
        
    end
      
    
    %% Summary
    
    disp(sprintf('   Processed %d clusters.',N0));
    disp(sprintf('   Retained %d clusters (%d%%).\n',N1,round(N1/N0*100)));
    
    
end
