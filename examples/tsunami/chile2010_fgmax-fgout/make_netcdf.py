
"""
Test script to read in a sequence of fgout frames, and write a netCDF file
containing only eta from each fgout frame, along with the final topography B.

Then print info about it, and read it back in.
With datatype = 'f8' there should be no difference from the original.
With datatype = 'f4' (32-bit floats) the netCDF file is only half as large.

This example should be turned into a notebook with more discussin...
"""

from clawpack.geoclaw import fgout_tools

try:
    import netCDF4
except:
    'You must install netCDF4 in order to write a netCDF file'
    raise

fgno = 1  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output

fgframenos = range(1,5)  # frames of fgout solution to write to netCDF file

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 

fgout_frames = []
for fgframeno in fgframenos:
    fgout = fgout_grid.read_frame(fgframeno)
    fgout_frames.append(fgout)
    
datatype = 'f4'  # 'f4' or 'f8' 

fgout_tools.write_netcdf(fgout_frames,
                        fname_nc='fgout_frames.nc',
                        qois=['eta'],
                        datatype=datatype,
                        include_B0=True,
                        include_Bfinal=False,
                        description='Chile 2010 test problem',
                        verbose=True)
                
# ================================================
# Print out info about this file:
fgout_tools.print_netcdf_info('fgout_frames.nc') 

# ================================================
# Read back the saved fgout frames as arrays:

x,y,t,qoi_arrays = fgout_tools.read_netcdf_arrays('fgout_frames.nc',
                                                  ['eta','B0'], verbose=True)
                                                  
B0 = qoi_arrays['B0']
eta = qoi_arrays['eta']

# Compare h in first frame with the original data:
h0_original = fgout_frames[0].h
h0 = eta[0,:,:] - B0

dh0 = abs(h0 - h0_original).max()
print("With datatype='%s', abs(h0 - h0_original).max() = %.3e" % (datatype,dh0))  

    