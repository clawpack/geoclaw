"""
Marching front algorithm to select points from a topo file based on
various criteria.  See select_by_flooding docstring 

"""

from __future__ import print_function
from pylab import *
import sys
from clawpack.geoclaw import topotools



def select_by_flooding(Ztopo, mask=None, prev_pts_chosen=None,
                       Z1=-5., Z2=0., max_iters=None, verbose=False):
    """
    Uses Ztopo as the topography DEM.

    Returns pts_chosen, a 1/0 array of the same shape as Ztopo, indicating
    whether each point is chosen (1) or not (0), based on various criteria.

    If Z1 <= Z2 then points where Ztopo < Z1 are first selected and 
    then a marching algorithm is used to select neighboring points that 
    satisfy Ztopo < Z2. 

    Think of chosen points as "wet" and unchosen points as "dry" 
    and new points are flooded if at least one neighbor is "wet" and 
    the topography is low enough. Starting from deep water (e.g. Z1 = -5)
    this allows flooding up to MHW (Z2 = 0) without going past dikes that
    protect dry land with lower elevation.
     
    However, this can also be called with Z1 > Z2, in which case points
    where Ztopo > Z1 are first selected and then the marching algorithm is 
    used to select neighboring points that satisfy Ztopo > Z2.
    This is useful to select offshore points where the water is shallow,
    (and that are connected to the shore by a path of shallow points) 
    e.g. on the continental shelf. 

    If max_iters=None then iterate to convergence.

    Otherwise, at most max_iters iterations are taken, so setting
    this to a small value only selects points within this many grid
    points of where Ztopo < Z1, useful for buffering.

    If prev_pts_chosen is None we are starting from scratch, otherwise
    we possibly add additional chosen points to an existing array.
    Points where prev_pts_chosen[i,j]==1 won't change but those ==0 may be
    changed to 1 based on the new criteria.

    If mask==None or mask==False then the entire array is subject to 
    selection.  If mask is an array with mask[i,j]==True at some points,
    then either:
        - these masked points are marked as not selected (pts_chosen[i,j]=0)
          if prev_pts_chosen==None, or
        - these masked points are not touched if prev_pts_chosen is an
          array (so previous 0/1 values are preserved).

    Basic strategy (assuming Z1 <= Z2):
        Outside the masked region, if any, initialize pt_chosen to
            pt_chosen=1 where Ztopo < Z1
            pt_chosen=0 where Ztopo >= Z2
            pt_chosen=-5 for intermediate points to indicate unknown        
        Let water intrude to level Z2
            (setting pt_chosen=1 where Ztopo < Z2 and a neighbor is wet)
            This is done by a fast marching iteration in which we keep track 
            of the cells at the front to check next, those for which pt_chosen==1
            and at least one neighbor has pt_chosen=-5 (so sum of neighbors < 0)  
        Do this iteration for at most max_iters or until no new points are added
        
        At the end, any point that still has pt_chosen==-5 is set to 0 
        (not chosen).  (e.g. these correspond to dry points below MHW not
        connected to the coast.)
        
    Returns: pt_chosen array of same shape as Ztopo with values 1 (where chosen)
        or 0 (where not chosen).  

    To use as a mask to plot only the chosen points, use e.g.
        Zmasked = ma.masked_array(Ztopo, mask=logical_not(pt_chosen))
        
    Sample applications:
     1. With the default parameters, choose all points with Ztopo < 0 and that 
        can be connected by a path of grid points below this elevation to 
        "deep water" where Ztopo < -5.
    """
    
    # determine if we should first select points < Z1 and march up to Z2,
    # of if we should initially select points > Z1 and march down to Z2:
    increasing = (Z1 <= Z2)
    
    if mask is None:
        mask = zeros(Ztopo.shape)
    
    if prev_pts_chosen is not None:
        # reset previous unchosen points to unset for further consideration,
        # but only in regions not masked out
        reset = logical_and(prev_pts_chosen==0, logical_not(mask))
        pt_chosen = where(reset, -5, prev_pts_chosen)
        if verbose:
            print('Resetting %i points to unset' % reset.sum())    
    else:
        # set pt_chosen=1 if Z<Z1 and =0 if Z>Z2, else =-5 (unset)
        if increasing:
            cond1 = logical_and(Ztopo < Z1, logical_not(mask))
            cond0 = logical_and(Ztopo >= Z2, logical_not(mask))
        else:
            cond1 = logical_and(Ztopo > Z1, logical_not(mask))
            cond0 = logical_and(Ztopo <= Z2, logical_not(mask))
        pt_chosen = where(cond1, 1, -5)  # mark as chosen or unknown
        pt_chosen = where(cond0, 0, pt_chosen) # mark as unchosen
    
    pt_chosen_neighbor_sum = zeros(pt_chosen.shape)
    pt_chosen_neighbor_sum[:,1:] += pt_chosen[:,:-1]
    pt_chosen_neighbor_sum[:,:-1] += pt_chosen[:,1:]
    pt_chosen_neighbor_sum[1:,:] += pt_chosen[:-1,:]
    pt_chosen_neighbor_sum[:-1,:] += pt_chosen[1:,:]

    
    # Determine points at front to use in next iteration:
    # Points that have been chosen and for which at least one neighbor is
    # unmarked and needs to be determined.
    set_next_i, set_next_j = where(logical_and(pt_chosen == 1, \
                               pt_chosen_neighbor_sum < 0))

    nfront = len(set_next_i)
    if verbose:
        print('Initially: %i cells out of %i cells on frontier' \
                % (nfront, prod(Z.shape)))

    if max_iters is None:
        max_iters = prod(Ztopo.shape)  # upper bound, should need much fewer
    
    print('Selecting points with Z1 = %g, Z2 = %g, max_iters=%i' \
          % (Z1,Z2,max_iters))

    iter_count = 0
    maxi = Ztopo.shape[0] - 1
    maxj = Ztopo.shape[1] - 1

    for k in range(max_iters):
        if len(set_next_i) == 0:
            # done
            break
        iter_count += 1
        if verbose==2:
            print('Iteration %3i: there are %7i points on frontier' \
                   % (k, len(set_next_i)))
        new_set_next_i = []  # for the next iteration
        new_set_next_j = []  # for the next iteration
        for (i,j) in zip(set_next_i,set_next_j):
            for (ii,jj) in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if (ii>=0) and (ii<=maxi) and (jj>=0) and (jj<=maxj):
                    if increasing:
                        cond1 = (Ztopo[ii,jj] < Z2)
                    else:
                        cond1 = (Ztopo[ii,jj] > Z2)
                    if (pt_chosen[ii,jj] < 0) and cond1 \
                                and (not mask[ii,jj]):
                        pt_chosen[ii,jj] = 1
                        new_set_next_i.append(ii)
                        new_set_next_j.append(jj)
                else:
                    # (ii,jj) outside grid
                    pass
        set_next_i = new_set_next_i
        set_next_j = new_set_next_j
        
            
    # any remaining unset points are dry points below Z2
    pt_chosen = where(pt_chosen<0, 0, pt_chosen)

    print('Done after %i iterations with %i points chosen' \
          % (iter_count,pt_chosen.sum()))

    return pt_chosen

    
    
