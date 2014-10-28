"""
Useful things for plotting 1D GeoClaw results.
"""
import setrun
rundata=setrun.setrun()
drytol_default = rundata.geo_data.dry_tolerance

def topo(current_data):
   """
   topo is assumed to be aux[0,:]
   """
   aux = current_data.aux
   topo = aux[0,:]
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
   eta = aux[0,:] + h
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
    drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
    q = current_data.q
    aux = current_data.aux
    h = q[0,:]
    eta = aux[0,:] + h
    water = ma.masked_where(h<=drytol, eta)
    return water


