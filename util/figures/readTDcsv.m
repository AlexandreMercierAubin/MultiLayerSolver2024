%readTDcsv(fileNames, maxFrames,h_in, plotNoContact)
function readTDcsv(fileNames, maxFrames,h_in, plotNoContact)
    if nargin < 2
        maxFrames = -1;
    end
    if nargin < 3
        h=1;%print time in step count
    else
        h = h_in;
    end

    if nargin < 4
        plotNoContact = false;
    end

    fontSize = 8;
    fontName = 'Linux Biolinum O';
    hfig = figure('Renderer', 'painters', 'Units','Inches', 'Position', [0,0,3.3125,2]);%
    set(gcf,'color','w');
    hold on;
    %for i = 1:numel(fileNames)
       myTable1 = readtable(fileNames(1));
       myTable2 = readtable(fileNames(2));
       if maxFrames > 0
           myTable1 = myTable1(1:maxFrames, :);
           myTable2 = myTable2(1:maxFrames, :);
       end
       %col 10 is the total sim time and 11 is the time to compute the
       %energy
       times1 = table2array(myTable1(:,10)) - table2array(myTable1(:,11));
       sumTime1 = sum(times1)
       times2 = table2array(myTable2(:,10)) - table2array(myTable2(:,11));
       sumTime2 = sum(times2)
       timeNoContact1 = times1 - table2array(myTable1(:,15));
       sumTimeNoContact1 = sum(timeNoContact1)
       timeNoContact2 = times2 - table2array(myTable2(:,15));
       sumTimeNoContact2 = sum(timeNoContact2)

       totSpeedup = sumTime2/sumTime1
       totSpeedupNoContact = sumTimeNoContact2/sumTimeNoContact1

       speedup = times2./times1;
       meanspeed = mean(speedup)

       speedupNoContact = timeNoContact2./timeNoContact1;
       meanspeedNoContact = mean(speedupNoContact)

       endFrame = (size(times1,1)-1);
       if ~plotNoContact
            plot([0:endFrame]'.*h, speedup);
            plot([0, (endFrame*h)], [meanspeed, meanspeed], 'g');
       else
            plot([0:endFrame]'.*h, speedupNoContact);
            plot([0, (endFrame*h)], [meanspeedNoContact, meanspeedNoContact], 'g');
       end
       string(num2str(speedup));
       ylabel('Speedup Factor');
       if h == 1
            xlabel('Time Step');
       else
           xlabel('Time (seconds)');
       end

    %end
    
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