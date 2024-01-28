module amber

use constants,  only :  dp

use input,      only :  Var_Input

implicit none
save

type, public :: Var_Amber
   integer,     dimension(:)      :: num
   real(kind=dp), dimension(:,:) :: 


end type Var_Amber

logical, private :: load_amber = .false.

private ::

contains

 function amb_load(vinp) result(vamb)
   implicit none
   
   if (load_amber) return

   



   load_amber = .true.

 end function amb_load


end module amber
