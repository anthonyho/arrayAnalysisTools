function makepretty(varargin)
% MAKEPRETTY Makes figures look prettier
%    MAKEPRETTY(lineWidthPlot,markerSize,fontSizeAxes,fontSizeOther,
%    lineWidthAxes,backgroundColor) sets the line width and marker size of 
%    the plotted line, the font size of the axes labels, the font size of 
%    other text within the figure, the line width of the axes, and the 
%    background color of a figure according to the user-specific parameters. 
%
%    Default parameters are:
%       lineWidthPlot = 2
%       markerSize = 8
%       fontSizeAxes = 22
%       fontSizeOther = 22
%       lineWidthAxes = 2
%       backgroundColor = 'w'
%
%    MAKEPRETTY(AX,...) resets the parent figure of AX instead of GCF.
%
% v1, Anthony Ho, 10/22/2014


    %% Reading arguments
    
    % Parse possible Axes input
    error(nargchk(0,inf,nargin,'struct'));
    [cax,args,nargs] = axescheck(varargin{:});
    
    if isempty(cax)
        figureHandle=gcf;
    else
        figureHandle=ancestor(cax,'figure');
    end
    
    % Reading arguments
    if nargs<1
        lineWidthPlot = 2;
    else
        lineWidthPlot = args{1};
    end
    if nargs<2
        markerSize = 8;
    else
        markerSize = args{2};
    end
    if nargs<3
        fontSizeAxes = 22;
    else
        fontSizeAxes = args{3};
    end
    if nargs<4
        fontSizeOther = 22;
    else
        fontSizeOther = args{4};
    end
    if nargs<5
        lineWidthAxes = 2;
    else
        lineWidthAxes = args{5};
    end
    if nargs<6
        backgroundColor = 'w';
    else
        backgroundColor = args{6};
    end
    

    %% Getting the right handles
    
    allAxes = findall(figureHandle,'type','axes');
    allLines = findall(figureHandle,'type','line');
    allText = findall(figureHandle,'type','text');

    
    %% Setting the properties of the figure
    
    % Set line width of the axes
    set(allAxes,'lineWidth',lineWidthAxes);
    
    % Set line width of the plotted lines
    set(allLines,'lineWidth',lineWidthPlot);
    
    % Set markers size of the plotted lines
    set(allLines,'markerSize',markerSize);
    
    % Set font size of the axis labels
    set(allAxes,'fontSize',fontSizeAxes);
    
    % Set font size of xlabels, ylabels, titles, etc
    set(allText,'fontSize',fontSizeOther);
    
    % Set background background
    set(figureHandle,'color',backgroundColor);
    
     
end