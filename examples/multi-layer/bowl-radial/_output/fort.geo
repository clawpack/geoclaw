  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:  0.90000000000000002        1.0000000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           1
    sea_level:   0.0000000000000000     
  
    coriolis_forcing: F
    theta_0:   0.0000000000000000     
    friction_forcing: T
    manning_coefficient:   2.5000000000000001E-002
    friction_depth:   20.000000000000000     
  
    dry_tolerance:  0.10000000000000001     
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:   1.0000000000000000E-002
    speed_tolerance:   1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000     
    maxleveldeep:           3
    depthdeep:   100.00000000000000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            0
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    /Users/hudaqureshi/clawpack/geoclaw/examples/multi-layer/bowl-radial/bowl.topotype2                                                                   
   itopotype =            2
   minlevel, maxlevel =            1           3
   tlow, thi =    0.0000000000000000        10000000000.000000     
   mx =          201   x = (  -1000.0000000000000      ,   1000.0000000000000      )
   my =          201   y = (  -1000.0000000000000      ,   1000.0000000000000      )
   dx, dy (meters/degrees) =    10.000000000000000        10.000000000000000     
  
   Ranking of topography files  finest to coarsest:            1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
    min_level, max_level, qinit_fname:
           0           1 /Users/hudaqureshi/clawpack/geoclaw/examples/multi-layer/bowl-radial/hump.xyz                                                                         
   
 Reading qinit data from
 /Users/hudaqureshi/clawpack/geoclaw/examples/multi-layer/bowl-radial/hump.xyz                                                                         
   
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:                  Infinity
    eigen_method:           2
    inundation_method:           2
    dry_tolerance:  0.10000000000000001       0.10000000000000001     
