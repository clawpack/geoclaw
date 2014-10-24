subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use gauges_module
   use geoclaw_module

   implicit none

   call set_gauges()
   call set_geo()

end subroutine setprob
