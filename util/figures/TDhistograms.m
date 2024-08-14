function TDHistograms(scene, path, tests, residualHistBinSize, stackHistograms)
close all;
figure;

if nargin < 4
    residualHistBinSize = -1; %automatic
end
if nargin < 5
    stackHistograms = false; %automatic
end

% scene = "PachinkoPill_04-16-2024_04-13";
comparisonNames = ["XPBD", "Residual Velocity Layers"];

if ~stackHistograms
    comparisonTimes = zeros(numel(comparisonNames),0);
    for i = 1:numel(tests)
        filename = path +string(tests(i)) +scene+".mat";
        % filename = path +scene+".mat";
        [td] = loadTDcellarray(filename);
        entries = numel(td{1}.log(6,:));
    
        %allocating space
        comparisonTimes = [comparisonTimes, zeros(numel(comparisonNames),entries)];
        for plotnum = 1:numel(td)
            comparisonTimes(plotnum,(end-entries+1):end) = td{plotnum}.log(6,:);
        end
    end
    
    k=numel(comparisonNames);
    hold("on");
    if residualHistBinSize > 0
        for plotnum = 1:k
            fig = histogram(comparisonTimes(plotnum,:), 'BinWidth', residualHistBinSize );
        end
    else
        hist1 = histogram(comparisonTimes(1,:));
        residualHistBinSize = hist1.BinWidth;
        for plotnum = 2:k
            histogram(comparisonTimes(plotnum,:), 'BinWidth', residualHistBinSize);
        end
    end
    title("Histograms of Runtime for "+string(cat(tests))+"% Residual Error Improvement");
    fontSize = 12;
    fontName = 'Linux Biolinum O';
    set(gca, 'fontSize', fontSize, 'fontName', fontName);
else
    tileHnandle = tiledlayout(numel(tests),1);
    for i = 1:numel(tests)
        nexttile;
        filename = path +string(tests(i)) +scene+".mat";
        [td] = loadTDcellarray(filename);
        for plotnum = 1:numel(comparisonNames)
            histogram(td{plotnum}.log(6,:), 'BinWidth', residualHistBinSize);
            hold on;
            xlim([0.15,0.6]);
            ylim([0,4000]);
            fontSize = 12;
            fontName = 'Linux Biolinum O';
            set(gca, 'fontSize', fontSize, 'fontName', fontName);
        end
    end
end


w = warning ('off','all');
[~] = legend(comparisonNames([1,2]), 'Location','northeast');
w = warning ('on','all');

xlabel("Time (s)");
ylabel("Frames");


saveas(gca,path+scene+".svg");
end
