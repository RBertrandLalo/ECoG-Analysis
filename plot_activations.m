
    addpath([cd,'/activeBrain']);

    load([cd,'/brain_models/','AMC078.mat']);
    
    viewstruct.viewvect     = [90, 0]; % view
    viewstruct.lightpos     = [150, 0, 0]; % light
    
    viewstruct.what2view = {'brain', 'activations'}; % remove 'activations' to only visualize the brain
        
    tala.activations = randn(size(tala.electrodes, 1),1); % tala is a structure with the activations and electrodes coordinates
            
    % sets min and max of colorbar
    cmapstruct.cmin = -max(abs(tala.activations));
    cmapstruct.cmax = max(abs(tala.activations));
    cmapstruct.enablecolorbar = 1;
    cmapstruct.enablecolormap = 1;
    cmapstruct.fading         = 0;
    c_steps                   = 256;
    cmapstruct.cmap           = jet(c_steps);
    
    figure

    % plots viewstruct.what2view (i.e. brain w/wo activations)
    activateBrain(cortex, vcontribs, tala, ix, cmapstruct, viewstruct); 
   
    % to view the electrodes
    for el = 1:size(tala.electrodes, 1)
        
        % size by default: 
        % size_el = 0.7;
        
        % size proportional to the activation:
        size_el = double(tala.activations(el));
        
        plotBalls(tala.electrodes(el, :), 'k', size_el) % be sure not to have any negative values for the activations
        
    end
   