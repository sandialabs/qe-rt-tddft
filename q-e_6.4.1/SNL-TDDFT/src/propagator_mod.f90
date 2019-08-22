MODULE propagator_mod
  ! 
  ! ... Class that provides a time propagator for TDDFT
  ! 
  USE kinds,  ONLY : dp
  USE it_solver_mod
  ! 
  SAVE
  ! 
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

SUBROUTINE set_implicit_flag(this, implicit_flag)
    !
    ! ... interface for setting the implicit propagation flag's value
    !
    IMPLICIT NONE
    ! input variables
    CLASS(propagator_type), INTENT(INOUT) :: this
    LOGICAL :: implicit_flag
   
    this%implicit_flag = implicit_flag

    RETURN

END SUBROUTINE set_implicit_flag

SUBROUTINE set_nstages(this, nstages)
    !
    ! ... interface for setting the number of stages for a propagator
    ! 
    IMPLICIT NONE
    ! input variables
    CLASS(propagator_type), INTENT(INOUT) :: this
    INTEGER :: nstages

    this%nstages = nstages
  
    RETURN

END SUBROUTINE set_nstages

SUBROUTINE electron_step(this)
    ! 
    ! ... interface for marching the electrons forward by one step
    ! 
    IMPLICIT NONE
    ! input variables
    CLASS(propagator_type), INTENT(INOUT) :: this

    RETURN

END SUBROUTINE electron_step

END MODULE propagator_mod
