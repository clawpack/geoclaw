function setcolors(p,x,y,z,q)
%
% SETCOLORS specifies user-defined mapping between q data and color map
%
%      User-defined routine which maps q values to the colormap.  To
%      indicate that this routine should be called, the user should set
%      the flag 'UserColorMapping' to 1 in SETPLOT2 or SETPLOT3 files.
%
%      The syntax for this routine is
%
%              setcolors(p,x,y,z,q);
%
%      where p is the handle to the patch whose colors are being set,
%      (x,y,z) are the Cartesian locations of cell centers with corresponding
%      data values in q.  If a mapped grid or manifold is being plotted,
%      (x,y,z) are given in the Cartesian locations, NOT the physical
%      locations of cell centers.
%
%      Possible uses for this routine include :
%
%             -- Use different colormaps in different domain regions
%             -- Mask out embedded boundary regions.
%             -- Visualize parallel partitions
%
%      To visualize under and over shoots, see UNDEROVER.
%
%      Example : Use the 'spring' colormap for the left half of a region
%      and the 'winter' color map for the right half.
%
%              c1 = spring;
%              c2 = winter;
%
%              cm_new = [c1; c2];
%
%              n1 = length(c1);
%              n2 = length(c2);
%
%              % Map q(m1) to indices [1,n1] and q(m2) to indices [n1+1,n1+n2]
%              idx = 0*q + nan;   % To make sure we don't leave any q values out.
%
%              m1 = x > 0.5;      % Left and right halves of computational domains
%              m2 = x <= 0.5;
%
%              qmin = -0.01;      % global min and max values for data
%              qmax = 1.01;
%
%              q(q == qmax) = qmax*(1-1e-8);     % So floor works below.
%              slope = (q - qmin)./(qmax-qmin);  % slopes must be in [0,1]
%
%              idx(m1) = floor(1 + slope(m1)*(n1 - 1));
%              idx(m2) = n1 + floor(1 + slope(m2)*(n2 - 1));
%
%              set(p,'cdata',idx);
%              fv = get(p,'FaceVertexCData');
%
%              if (sum(isnan(fv)) > 0)
%                error('setcolors : Not all values in q have been mapped')
%              end;
%
%              set(p,'FaceVertexCData',cm_new(fv,:));
%              set(p,'facecolor','flat');
%
%
%      See also PATCH, COLORMAP, UNDEROVER, UNDEROVER_COLORBAR,
%      MULTICOLORMAP, MULTICOLORMAP_COLORBAR.
%

cm = colormap;

% Map nan values to brown (land)
cnew = [[101,67,33]/255; cm];

qmin = min(q(:));
qmax = max(q(:));

qmin = -0.1;
qmax = 0.1;

q(q < qmin) = qmin;
q(q > qmax) = qmax;

slope = (q - qmin)./(qmax-qmin);

m1 = isnan(q);
m2 = ~m1;

N = length(cm);

idx = 0*q + nan;
idx(m1) = 1;
idx(m2) = round(2 + slope(m2)*(N-1));

set(p,'cdata',idx);
fv = get(p,'FaceVertexCData');

if (sum(isnan(fv)) > 0)
    fprintf('setcolors : Not all values in q have been mapped\n')
end

if (min(fv(:)) < 1)
    fprintf('Error in min index : \n');
end


if (max(fv(:)) > length(cnew))
    fprintf('Error in max index : \n');
end

set(p,'FaceVertexCData',cnew(fv,:));
set(p,'facecolor','flat');

set(p,'cdatamapping','direct');

end