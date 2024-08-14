function readTDcsvLog(fileNames, colors,  maxFrames,h_in, plotNoContact)
    if nargin < 3
        maxFrames = -1;
    end
    if nargin < 4
        h=1;%print time in step count
    else
        h = h_in;
    end
    if nargin < 5
        plotNoContact = false;
    end

    fontSize = 8;
    fontName = 'Linux Biolinum O';
    hfig = figure('Renderer', 'painters', 'Units','Inches', 'Position', [0,0,3.3125,2]);%
    ax1 = axes('Parent', hfig);
    set(gcf,'color','w');
    hold(ax1,"off");
    for i = 1:numel(fileNames)
       myTable1 = readtable(fileNames(i));
       if maxFrames > 0
           myTable1 = myTable1(1:maxFrames, :);
       end
       %col 10 is the total sim time and 11 is the time to compute the
       %energy
       times1 = table2array(myTable1(:,10)) - table2array(myTable1(:,11));
       timeNoContact1 = times1 - table2array(myTable1(:,15));
       timesSum = sum(times1);
       timesNCSum = sum(timeNoContact1);
       disp(fileNames(i));
       disp("sum time:"+string(timesSum));
       disp("sum time NC:"+string(timesNCSum));

       endFrame = (size(times1,1)-1);
       if ~plotNoContact
           semilogy(ax1,[0:endFrame]'.*h, times1, "Color", colors(i));
       else
           semilogy(ax1,[0:endFrame]'.*h, timeNoContact1, "Color", colors(i));
       end
       hold(ax1,"on");
       ylabel('Log of Wall Clock Time');
       if h == 1
            xlabel('Time Step');
       else
           xlabel('Time (seconds)');
       end

    end
    
    handle = gca;
    handle.Color = [0.95, 0.95, 0.95];
    handle.XLim = [0, size(times1,1)*h];
    handle.TickLength = [0,0];
    handle.YGrid = 'on';
    handle.GridColor = [1 1 1];
    handle.XColor = [0 0 0];
    set(gca, 'fontSize', fontSize, 'fontName', fontName);
%     set(hfig,'Units','Inches');
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperPosition',[0,0,3.3125,2]);%
    printSvg = strrep(fileNames(1),'.csv','.svg');
%     exportgraphics(gca,"benchmark/"+printPdf,'ContentType','vector');
    saveas(gca,"benchmark/"+printSvg)
    %legend(fileNames);
    
    
end