  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:   1025.0000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           1
    sea_level:  -10.000000000000000     
  
    coriolis_forcing: F
    theta_0:   0.0000000000000000     
    friction_forcing: F
    manning_coefficient: not used
    friction_depth: not used
  
    dry_tolerance:   1.0000000000000000E-003
  
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
    
    /Users/hudaqureshi/clawpack/geoclaw/examples/tsunami/bowl-slosh/bowl.topotype2                                                                        
   itopotype =            2
   minlevel, maxlevel =            1          10
   tlow, thi =    0.0000000000000000        10000000000.000000     
   mx =          200   x = (  -2.0000000000000000      ,   2.0000000000000373      )
   my =          200   y = (  -2.0000000000000000      ,   2.0000000000000373      )
   dx, dy (meters/degrees) =    2.0100502512563002E-002   2.0100502512563002E-002
  
   Ranking of topography files  finest to coarsest:            1
  
  
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
