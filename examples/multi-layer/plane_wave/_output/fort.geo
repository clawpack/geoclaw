  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:   922.50000000000000        1025.0000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           1
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
    num dtopo files =            0
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    /Users/hudaqureshi/clawpack/geoclaw/examples/multi-layer/plane_wave/topo.tt2                                                                          
   itopotype =            2
   minlevel, maxlevel =            1           5
   tlow, thi =    0.0000000000000000        10000000000.000000     
   mx =           88   x = (  -2.0000000000000000      ,   3.9999999999999734      )
   my =           88   y = (  -2.0000000000000000      ,   3.9999999999999734      )
   dx, dy (meters/degrees) =    6.8965517241379004E-002   6.8965517241379004E-002
  
   Ranking of topography files  finest to coarsest:            1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
  epsilon =    2.0000000000000000E-002
  init_location =  -0.10000000000000001        0.0000000000000000     
  wave_family =            4
  angle =    0.0000000000000000     
  sigma =    2.0000000000000000E-002
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:  0.94999999999999996     
    eigen_method:           2
    inundation_method:           2
    dry_tolerance:   1.0000000000000000E-003   1.0000000000000000E-003
