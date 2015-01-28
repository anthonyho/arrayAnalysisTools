classdef singleClustersData
    
    properties
 
       name
       ID
       
       nClusters
       nTimePoint
       
       data
       rData
       nData
       
       avgRData
       medRData
       
       avgNData
       medNData
       
       xAxis
       
       thresholdPercentile = 5;
       
    end % End of properties
    
    methods
        
        %% Constructor method for the singleClustersData class
        function SCD = singleClustersData(name,ID,filepath,thresholdPercentile)
            
            % Defining threshold percentile
            if nargin>=4
                SCD.thresholdPercentile = thresholdPercentile;
            end

            % Defining names and ID of the class
            SCD.name = name;
            SCD.ID = ID;
            
            % Import data from filepath
            disp(sprintf('Processing %s...',filepath));

            SCD.data = importdata(filepath);
            SCD.rData = SCD.data;
            SCD.nData = SCD.data;

            % Size of dataset
            nClusters0 = size(SCD.data,1);
            SCD.nTimePoint = size(SCD.data,2);
            
            % Initialize avgeraged and median data 
            SCD.avgRData = zeros(1,SCD.nTimePoint);
            SCD.avgNData = zeros(1,SCD.nTimePoint);
            SCD.medRData = zeros(1,SCD.nTimePoint);
            SCD.medNData = zeros(1,SCD.nTimePoint);

            % Computing threshold for the first time point
            threshold = prctile(SCD.rData(:,1),SCD.thresholdPercentile);

            % Normalizing to the first time point of each cluster
            for j = 1:nClusters0
                if SCD.rData(j,1)>threshold
                    amp = SCD.rData(j,1);
                    SCD.nData(j,:) = SCD.rData(j,:)/amp;
                end
            end

            % Removing clusters below threshold
            rowsToRemove = any(SCD.rData(:,1)<threshold,2);
            SCD.rData(rowsToRemove,:) = [];
            SCD.nData(rowsToRemove,:) = [];

            % Removing clusters with entrances equal to 0
            rowsToRemove = any(SCD.rData==0,2);
            SCD.rData(rowsToRemove,:) = [];
            SCD.nData(rowsToRemove,:) = [];

            % Computing new number of clusters
            SCD.nClusters = size(SCD.rData,1);

            % Computing the average/median raw/normalized data
            for i = 1:SCD.nTimePoint
                SCD.avgRData(i) = mean(SCD.rData(:,i));
                SCD.medRData(i) = median(SCD.rData(:,i));
                SCD.avgNData(i) = mean(SCD.nData(:,i));
                SCD.medNData(i) = median(SCD.nData(:,i));
            end

            % Printing summary
            disp(sprintf('   Processed %d clusters.',nClusters0));
            disp(sprintf('   Retained %d clusters (%d%%).\n', ...
                SCD.nClusters,round(SCD.nClusters/nClusters0*100)));

        end
        
        %% Method to save median data
        
        

        
        
    end % End of methods
    
end % End of class 
