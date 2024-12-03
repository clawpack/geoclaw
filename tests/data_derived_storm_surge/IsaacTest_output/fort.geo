  
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
  
    coriolis_forcing: T
    theta_0:   0.0000000000000000     
    friction_forcing: T
    manning_coefficient:   2.5000000000000001E-002
    friction_depth:   10000000000.000000     
  
    dry_tolerance:   1.0000000000000000E-002
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:   1.0000000000000000     
    speed_tolerance:   1.0000000000000000        2.0000000000000000        3.0000000000000000        4.0000000000000000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            0
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    /home/catherinej/clawpack/geoclaw/scratch/gulf_caribbean.tt3                                                                                          
   itopotype =            3
   mx =         1471   x = (  -99.000000000000000      ,  -49.999990199999999      )
   my =          721   y = (   8.0000000000000000      ,   32.000004799999999      )
   dx, dy (meters/degrees) =    3.3333340000000003E-002   3.3333340000000003E-002
  
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
    dry_tolerance:   1.0000000000000000E-002
