function plot_gauges(plot_var)

if nargin < 1
    plot_var = 'v';
end

close all

% ------------------------------------------
geo_outdir = './_output';

plot_geo = true;
plot_obs_data = true;


% Scaling
tscale = 1/3600;     % Convert time to hours
if (plot_obs_data)
    tshift = 10*60;
end
vel_scale = 100;       % Convert velocity to cm/sec


% ----------------------------------
gdata = read_gauge_data();

t_idx = 2;
h_idx = 3;
hu_idx = 4;
hv_idx = 5;
eta_idx = 6;

switch plot_var
    case 'eta'
        pvidx = 1;
        ystr = 'Surface height';
        yla = [-3,3];    % Y axis limits
    case 'u'
        pvidx = 2;
        ystr = 'u-velocity';
        yla = [-275,275];
    case 'v'
        pvidx = 3;
        ystr = 'v-velocity';
        yla  = [-275,275];
    case 'speed'
        pvidx = 4;
        ystr = 'Speed';
        yla  = [-275,275];
    otherwise
        error('No valid plot_var was specified');
end

ph = [];
lstr = {};
num_gauges = length(gdata);
for i = 1:num_gauges
    clear ph lstr;
    g = gdata(i);
    
    if (g.id == 1123 && strcmpi(plot_var,'eta'))
        continue
    elseif (g.id == 5680 && ~strcmpi(plot_var,'eta'))
        continue
    end
    
    figure(100+i);
    clf;
    hold on;
    
    k = 1;
                    
    if (plot_geo)
        gname_geo = sprintf('%s/gauge%05d.txt',geo_outdir,g.id);
        if (exist(gname_geo,'file'))
            tseries_geo = importdata(gname_geo,' ',3);    
            t_geo = tscale*(tseries_geo.data(:,2) + tshift);  % Shift by 10 minutes
            eta_geo = tseries_geo.data(:,eta_idx);
            h_geo = tseries_geo.data(:,h_idx);
            u_geo = vel_scale*tseries_geo.data(:,hu_idx)./h_geo;
            v_geo = vel_scale*tseries_geo.data(:,hv_idx)./h_geo;
            speed_geo = sqrt(u_geo.^2 + v_geo.^2);
            pvars_geo = {eta_geo, u_geo, v_geo, speed_geo};
            
            pv_geo = pvars_geo{pvidx};
            
            hold on;    
            ph(k) = plot(t_geo,pv_geo,'r.-','linewidth',2,'markersize',8);
            lstr{k} = 'GeoClaw';
            k = k + 1;
        else
            fprintf('File %s does not exist\n',gname_geo);
        end
    end
                        
    if (plot_obs_data)
        hold on;
        pout = plot_obs(g.id,plot_var);
        if (pout ~= 0)
            ph(k) = pout;
            lstr{k} = 'Observations';
            k = k + 1;
        end
    end
           
    xl = xlim;
    ph(k) = plot(xl,0*xl,'k');
    lstr{k} = 'Sea level';
    
    title(sprintf('Gauge %d',g.id),'fontsize',18);
    xlabel('t (hours)','fontsize',16);
    ylabel(ystr,'fontsize',16);
    set(gca,'fontsize',16);
    legend(ph,lstr);
    set(gca,'box','on');
        
%{    
    if (~plot_obs_data)
        yl = [min([pv; pv_geo]), max([pv; pv_geo])];
        ym = mean(yl);
        ys = 1.1*diff(yl)/2;
        ylim(ym + ys*[-1,1]);
    else    
        set(gca,'xlim',[7.5,13]);
        set(gca,'ylim',xla);
        % set(gca,'ylim',[-3,3]);
        set(gcf,'position',[173   358   986   420]);
    end
%}    

    set(gca,'xlim',[7.5,13]);
    set(gca,'ylim',yla);
    % set(gca,'ylim',[-3,3]);
    set(gcf,'position',[173   358   986   420]);

    
    hold off;
    shg
       
end



end