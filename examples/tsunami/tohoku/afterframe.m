% Patches
hidegridlines(1:6);
showpatchborders(1:15);
% setpatchborderprops('linewidth',1);

colormap(parula);
colorbar;
tol = -0.8;
c1 = -0.1;
c2 = 0.1;
caxis([c1,c2]);

set(gca,'zlim',[-10,1])

% Adds gauges (doesn't replace the figure) 
add_gauges('geoclaw');
add_regions(t,'geoclaw');

fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);


% Axes
axis([132 210 9 53])
daspect([1 1 1]);
set(gca,'fontsize',16);

% Use these zoom regions to see finer resolution grids. 
% Zoom Frame = 2
% axis([134.6060,  164.9044,   28.2413,   45.3327])

% Zoom 1 (Frame = 17)
% axis([201.3305,  206.7376,   18.8768,  21.9269]);

% Zoom 2 (Frame 17)
% axis([202.2093, 204.2401, 20.3161, 21.4617])

% Zoom 3 (Frame 17)
% axis([202.9598,  203.7538,   20.6753,   21.1232]);

% Zoom 4 (Frame 17)
% axis([203.4689,  203.6265,   20.8651,   20.9540]);

% Zoom 5 (Frame 18)
% axis([203.5126,  203.5443,   20.8883,   20.9062]);


title(sprintf('Tohoku (%d) : t = %.2f (%.2f,%.2f)',Frame,t,qmin,qmax),'fontsize',18);

NoQuery = 0;
prt = false;
MaxFrames = 26;
if (prt)
    filename = sprintf('tohoku_%04d_fc.png',Frame);
    fprintf('Print file %s\n',filename);
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
