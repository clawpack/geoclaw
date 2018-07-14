  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.81000000000000     
    density water:   1025.00000000000     
    density air:   1.15000000000000     
    ambient pressure:   101300.000000000     
    earth_radius:   6367500.00000000     
    coordinate_system:           2
    sea_level:  0.280000000000000     
  
    coriolis_forcing: T
    theta_0:  0.000000000000000E+000
    friction_forcing: T
    manning_coefficient:  2.500000000000000E-002
    friction_depth:   10000000000.0000     
  
    dry_tolerance:  1.000000000000000E-002
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:   1.00000000000000     
    speed_tolerance:   1.00000000000000        2.00000000000000     
   3.00000000000000        4.00000000000000     
    maxleveldeep:           4
    depthdeep:   300.000000000000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            0
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    
 /rigel/apam/users/hq2152/clawpack/geoclaw/scratch/gulf_caribbean.tt3           
                                                                        
   itopotype =            3
   minlevel, maxlevel =            1           5
   tlow, thi =   -259200.000000000        86400.0000000000     
   mx =         1471   x = (  -99.0000000000000      ,  -49.9999902000000      )
   my =          721   y = (   8.00000000000000      ,   32.0000048000000      )
   dx, dy (meters/degrees) =   3.333334000000000E-002  3.333334000000000E-002
  
   Ranking of topography files  finest to coarsest:            1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
   qinit_type = 0, no perturbation
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:  0.950000000000000     
    eigen_method:           4
    inundation_method:           2
    dry_tolerance:  1.000000000000000E-002
