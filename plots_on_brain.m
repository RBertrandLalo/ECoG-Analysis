% load brain model
load('AMC078.mat'); 

% make sure activateBrain script is in path
addpath(genpath('activeBrain')); 

% change activateBrain settings
viewstruct.what2view = {'brain'};
cmapstruct.enablecolorbar = 0;

% plot brain
figure('Position',[100,100,1000,800]);
activateBrain( cortex, vcontribs, tala, ix, cmapstruct, viewstruct );

% plot for each electrode

inset_width = .035;
inset_height = .035;
xl = xlim;
yl = ylim;
zl = zlim;

channel2plot0 = [3 4 52];
channel2plot1 = [12:14];
channel2plot2 = [20:24];
channel2plot3 = [28:30];
channel2plot4 = [34:38];
channel2plot5 = [71:73];
channel2plot6 = [86:88];


for idx_ch = channel2plot,
 
    
    % this isn't perfect but it'll do - ljc
    % inset axes need positions based on 0:1 for x and y
    axes('Position',[(0.125+0.75*(tala.trielectrodes(idx_ch,2)-yl(1))/(yl(2)-yl(1)))...
        0.2+0.6*(tala.trielectrodes(idx_ch,3)-zl(1))/(zl(2)-zl(1)) inset_width inset_height]) % x, y, width, height
    
    % plot something
    imagesc(times, freqs(dispf), ersp{idx_ch});

    set(gca,'xtick',[],'ytick',[]);
    
end


plotSpheres(tala.electrodes(channel2plot0, :), 'w');
plotSpheres(tala.electrodes(channel2plot1, :), 'y');
plotSpheres(tala.electrodes(channel2plot2, :), 'g');
plotSpheres(tala.electrodes(channel2plot3, :), 'k');
plotSpheres(tala.electrodes(channel2plot4, :), 'r');
plotSpheres(tala.electrodes(channel2plot5, :), 'm');
plotSpheres(tala.electrodes(channel2plot6, :), 'b');

plotSpheres(tala.electrodes([21 23 28 29 35 36 37 ], :), 'r');



plotSpheres(tala.electrodes(35, :), 'y');

