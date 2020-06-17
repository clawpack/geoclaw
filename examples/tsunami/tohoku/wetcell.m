function eta = wetcell(data)
% read q data:
mbathy = 1;
mheight = 1;
h = data(:,mheight);

eta = data(:,4);

dry_tol = 0.001;
eta(h <= dry_tol) = nan;
end