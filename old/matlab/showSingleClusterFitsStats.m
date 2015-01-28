function showSingleClusterFitsStats(fitParams,fitParamsStruct)
% SHOWSINGLECLUSTERFITSSTATS Shows the histograms of the fit parameters of
% all clusters 
%
% Anthony Ho, 10/27/2014


    %% Transforming fit parameters
    
    nClusters = size(fitParams,1);
    nFitParams = size(fitParamsStruct,1);
    
    transformedFitParams = zeros(nClusters,nFitParams);
    
    for i = 1:nClusters
        for j = 1:nFitParams
            transformedFitParams(i,j) = fitParamsStruct{j,2}(fitParams(i,:));
        end
    end
    
    
    %% Plotting

    disp('Showing statistics of single-cluster fits...');
    
    hFig = figure;
    set(hFig,'Position',[5 1200 3000 350])
    %set(hFig,'Position',[5 1200 1200 350])

    for j = 1:nFitParams
        subplot(1,nFitParams,j);
        %histfd(transformedFitParams(:,j),0,fitParamsStruct{j,3},fitParamsStruct{j,4},fitParamsStruct{j,5});
        perfhist(transformedFitParams(:,j),fitParamsStruct{j,3}:fitParamsStruct{j,5}:fitParamsStruct{j,4});
        title(fitParamsStruct{j,1});
    end
    
    
end