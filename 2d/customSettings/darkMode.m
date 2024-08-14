function settings = darkMode(settings)
    %rendering options
        settings.elementLineColor = [1,1,1];
        settings.boundaryColor = 'white';
        settings.elementLineAlpha = 0.1;
        settings.elementFaceAlpha = 0.5; %TODO: maybe move this to the material properties with the rest of the color so we have per element alphas
        settings.backgroundColor = 'black';
        settings.plotRigidification = 0;
end

