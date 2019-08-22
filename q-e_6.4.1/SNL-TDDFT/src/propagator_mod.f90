MODULE propagator_mod
  ! 
  ! ... Class that provides a time propagator for TDDFT
  ! 
  USE kinds,  ONLY : dp
  ! 
  SAVE
  ! 
  PRIVATE 
  PUBLIC :: &
      set_dt_el,                    &! sets the size of the time step (attoseconds)
      set_implicit_flag,            &! sets the flag that indicates we're using an implicit solver
      set_nstages, 		    &! sets the number of stages
      electron_step		     ! advances the electrons by a single time step
  
  TYPE propagator_type
      INTEGER :: &
          nstages 		     ! number of stages in the time integrator
      LOGICAL :: &
          implicit_flag	             ! logical flag for implicit solvers (.TRUE.=>implicit, .FALSE.=>explicit)   
      REAL(dp) :: &
          dt_el                      ! the time step associated with the electrons   
      TYPE(it_solver_type), POINTER :: &
          implicit_solver	     ! used for implicit time stepping (e.g., Crank-Nicolson) 

  END TYPE propagator_type

CONTAINS

SUBROUTINE set_dt_el(this, dt_el)
    ! 
    ! ... interface for setting the size of the electronic time step (in attoseconds)
    ! 
    IMPLICIT NONE
    ! input variables
    CLASS(propagator_type), INTENT(INOUT) :: this
    REAL(dp), INTENT(IN) :: dt_el

    this%dt_el = dt_el

    RETURN

END SUBROUTINE set_dt_el

SUBROUTINE set_implicit_flag

END SUBROUTINE set_implicit_flag

END MODULE propagator_mod
