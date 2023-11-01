"""
Useful things for plotting 1D GeoClaw results.

Ultimately this should be combined with clawpack.visclaw.geoplot
and the functionality improved to handle either 1d or 2d more seamlessly.
"""

drytol_default = 1e-3

def topo(current_data):
   """
   Return topography = eta - h.
   Surface eta is assumed to be output as 4th column of fort.q files.
   """
   #aux = current_data.aux
   #topo = aux[0,:]
   h = current_data.q[0,:]
   eta = current_data.q[2,:]
   topo = eta - h
   return topo

def land(current_data):
   """
   Return a masked array containing the surface elevation only in dry cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
   q = current_data.q
   aux = current_data.aux
   h = q[0,:]
   eta = q[2,:]
   land = ma.masked_where(h>drytol, eta)
   return land


def depth(current_data):
   """
   Return a masked array containing the depth of fluid only in wet cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
   q = current_data.q
   h = q[0,:]
   depth = ma.masked_where(h<=drytol, h)
   return depth

def surface(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo.
    """
    from numpy import ma
    drytol = current_data.user.get('dry_tolerance', drytol_default)
    q = current_data.q
    aux = current_data.aux
    h = q[0,:]
    #eta = aux[0,:] + h
    eta = q[2,:]
    water = ma.masked_where(h<=drytol, eta)
    return water

def velocity(current_data):
    """
    Return a masked array containing the water velocity only in wet cells.
    Surface is eta = h+topo.
    """
    from numpy import ma
    drytol = current_data.user.get('dry_tolerance', drytol_default)
    q = current_data.q
    h = q[0,:]
    h_wet = ma.masked_where(h<=drytol, h)
    u_wet = q[1,:] / h_wet
    return u_wet

