function pout = plot_obs(id,plot_var)

if (id == 1123 && ~strcmpi(plot_var,'eta'))
    data = load('HAI1123_Kahului_harbor_detided.txt');
    t = data(:,1);
    m = find(t >= 0 & t <= 13);
    tobs = t(m);        % Time in hours
    uobs = data(m,2);   % detided velocities
    vobs = data(m,3);
    switch plot_var
        case 'u'
            pvobs = uobs;
        case 'v'
            pvobs = vobs;
        case 'speed'
            pvobs = sqrt(uobs.^2 + vobs.^2);
    end
elseif (id == 5680 && strcmpi(plot_var,'eta'))
    data = load('1615680_detided.txt');
    t = data(:,1);
    m = find(t >= 0 & t <= 13*3600);
    tobs = t(m);
    pvobs = data(m,2);
else
    pout = 0;
    return
end
pout = plot(tobs,pvobs,'k.-','linewidth',2,'markersize',20);
end