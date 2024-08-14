function figurePlotResidual(filename, iterationsPerLayer, cropRigidLayers, skipLineDraw)
    err = readmatrix(filename);
    figure;
    hold on;

    if cropRigidLayers
        xvals = [1,iterationsPerLayer(end)];
        xlim(xvals);
        semilogy( xvals(1):xvals(2), err(1,1:iterationsPerLayer(end)));
        for plotnum = 2:size(err,1)
            semilogy( xvals(1):xvals(2), err(plotnum,end-(iterationsPerLayer(end)-1):end));
        end
    else
        xvals = [1,size(err,2)];
        xlim(xvals);
        for plotnum = 1:size(err,1)
            semilogy( xvals(1):xvals(2), err(plotnum,:));
        end
        if ~skipLineDraw
            for layer = 1:numel(iterationsPerLayer)
                iterationsAtLayer = sum(iterationsPerLayer(1:layer)); 
                xline(iterationsAtLayer);  
                % text(iterationsAtLayer,10^residualRange(2),string(integrators{layerLines}.layers(layer)))
            end
        end
    end
    w = warning ('off','all');
    % [~] = legend(["XPBD", "Residual Velocity Layers"], 'Location','northeast');
    % [~] = legend(["Elastic","Strain rate", "MinvK eigs", "Random", "VStripes", "HStripes"], 'Location','northeast');
    [~] = legend(["elastic","d50", "d33", "d20","d10","d5"], 'Location','northeast');
    w = warning ('on','all');
    xlabel("Iterations");
    ylabel("|r|");
end

