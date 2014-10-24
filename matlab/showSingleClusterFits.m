function showSingleClusterFits(dataset,xAxis,fitParams)%,functionHandle)
% SHOWSINGLECLUSTERFITS Shows all clusters signals and their fits individually
%
% Anthony Ho, 10/23/2014


    %% Defining constants and such
    
    displayWidth = 8;
    displayHeight = 5;
    
    N = size(dataset,1);
    
    nPlotsInAFig = displayWidth*displayHeight;
    nFrames = floor(N/nPlotsInAFig);
    nRemainders = mod(N,nPlotsInAFig);
    
    xUB = (max(xAxis)-min(xAxis))*1.1;
    
    
    %% Plotting

    figure; 
    
    for i=1:nFrames
        
        clf;
        
        % Make sure to show the right number of plots in the last frame
        if i==nFrames
            nDisplayed = nRemainders;
        else
            nDisplayed = nPlotsInAFig;
        end
        
        % Plot all clusters at each frame
        for j=1:nDisplayed
            
            % Compute current index
            currentIndex = (i-1)*nPlotsInAFig+j;
            
            h = subplot(displayHeight,displayWidth,j);
            hold on;
            
            plot(xAxis,dataset(currentIndex,:),'o');
            plot(xAxis,doubleExponential4Distinct(fitParams(currentIndex,:),xAxis),'-');
            
            hold off;
            
            xlim([min(xAxis), xUB]);
            title(sprintf('eof = %0.2g',fitParams(currentIndex,end)));
            
            % Makes subplot more compact
            ti = get(h,'TightInset');
            op = get(h,'OuterPosition');
            set(h,'Position',[op(1)+ti(1) op(2)+ti(2) op(3)-ti(3)-ti(1) op(4)-ti(4)-ti(2)]);
            
        end
        
        makepretty(2,8,15,15,2,'w');
        
        k = waitforbuttonpress;
        
    end

    
end