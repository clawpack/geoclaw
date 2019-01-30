  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:   1025.0000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           2
    sea_level:   0.0000000000000000     
  
    coriolis_forcing: F
    theta_0:   0.0000000000000000     
    friction_forcing: T
    manning_coefficient:   2.5000000000000001E-002
    friction_depth:   1000000.0000000000     
  
    dry_tolerance:   1.0000000000000000E-003
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:  0.10000000000000001     
    speed_tolerance:   1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000     
    maxleveldeep:           3
    depthdeep:   100.00000000000000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            1
    fname:/Users/hudaqureshi/clawpack/geoclaw/scratch/dtopo_usgs100227.tt3                                                                                      
    topo type:           3
    minlevel, maxlevel:
           3           3
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    /Users/hudaqureshi/clawpack/geoclaw/scratch/etopo10min120W60W60S0S.asc                                                                                
   itopotype =            2
   minlevel, maxlevel =            1           3
   tlow, thi =    0.0000000000000000        10000000000.000000     
   mx =          361   x = (  -120.00000000000000      ,  -59.999999999879996      )
   my =          361   y = (  -60.000000000000000      ,   1.2000356264252332E-010 )
   dx, dy (meters/degrees) =   0.16666666666700000       0.16666666666700000     
  
   Ranking of topography files  finest to coarsest:            2           1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
   qinit_type = 0, no perturbation
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:  0.94999999999999996     
    eigen_method:           4
    inundation_method:           2
    dry_tolerance:   1.0000000000000000E-003
