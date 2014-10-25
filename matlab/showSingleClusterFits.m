function showSingleClusterFits(dataset,xAxis,fitParams,fitParamsStruct, ...
    functionHandle,colorByField,sortByField)
% SHOWSINGLECLUSTERFITS Shows all clusters raw signals and their fits 
% individually
%
% Anthony Ho, 10/23/2014


    %% Reading arguments
    
    if nargin<6
        colorByField = size(fitParams,2);
    end
    
    if nargin<7
        sortByField = 0;
    end
    

    %% Defining constants and such
    
    % Number of subplots per figure
    displayWidth = 8;
    displayHeight = 5;
    
    % Number of cliusters
    N = size(dataset,1);
    
    % Number of figures
    nPlotsInAFig = displayWidth*displayHeight;
    nFrames = floor(N/nPlotsInAFig);
    nRemainders = mod(N,nPlotsInAFig);
    
    % Upperbound of the x-axis
    xUB = (max(xAxis)-min(xAxis))*1.1;
    
    % Color gradient parameters
    color0 = [0 0 1];
    color1 = [1 0 0];
    colorLowerThreshold = 0;
    colorUpperThreshold = 90;
    colorParam0 = prctile(fitParams(:,colorByField),colorLowerThreshold);
    colorParam1 = prctile(fitParams(:,colorByField),colorUpperThreshold);
    
    % Sorting by fit parameter
    if sortByField~=0
        [fitParams,sortIndex] = sortrows(fitParams,sortByField);
        dataset = dataset(sortIndex,:);
    end
    
    % Textbox positions
    xpos = 0.7;
    ypos = 0.85;
        
    
    %% Plotting

    disp('Showing individual single cluster fits...');
    
    hFig = figure; 
    set(hFig,'Position',[5 1200 3000 1500])
    
    for i = 1:nFrames
        
        clf;
        
        % Make sure to show the right number of plots in the last frame
        if i==nFrames
            nDisplayed = nRemainders;
        else
            nDisplayed = nPlotsInAFig;
        end
        
        % Plot all clusters at each frame
        for j = 1:nDisplayed
            
            % Compute current index
            index = (i-1)*nPlotsInAFig+j;
            
            % Show subplot
            h0 = subplot(displayHeight,displayWidth,j);
            
            % Makes subplot more compact
            ti = get(h0,'TightInset');
            op = get(h0,'OuterPosition');
            set(h0,'Position',[ ...
                op(1)+ti(1) ...
                op(2)+ti(2) ...
                op(3)-ti(3)-ti(1) ...
                op(4)-ti(4)-ti(2)]);
            
            % Plot
            hold on;
            h1 = plot(xAxis,dataset(index,:),'o');
            h2 = plot(xAxis,functionHandle(fitParams(index,:),xAxis),'-');
            hold off;
            
            % X-axis
            xlim([min(xAxis), xUB]);
            
            % Show fit parameters in plot
            fitParamsText = '';
            for k = 1:size(fitParamsStruct,1)
                textLine = sprintf(' = %0.2g',fitParamsStruct{k,2}(fitParams(index,:)));
                fitParamsText = strcat(fitParamsText,fitParamsStruct{k,1},textLine,'\n');
            end
            text(xpos,ypos,sprintf(fitParamsText),'Units','normalized');
            
            % Color of the line and markers
            colorPercentile = (fitParams(index,colorByField)-colorParam0)/(colorParam1-colorParam0);
            set([h1, h2],'color',colorgradpercent(color0,color1,colorPercentile));
            
        end
        
        makepretty(2,8,16,13,2,'w');
        
        k = waitforbuttonpress;
        
    end

    
end