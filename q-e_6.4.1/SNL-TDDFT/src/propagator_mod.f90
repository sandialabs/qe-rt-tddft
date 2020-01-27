MODULE propagator_mod
  ! 
  ! ... Class that provides a time propagator for TDDFT
  ! 
  USE kinds,  ONLY : dp
  USE it_solver_mod
  ! 
  SAVE
  ! 
  TYPE, PUBLIC :: propagator_type
     
      INTEGER :: nstages                  ! number of stages in the time integrator
      LOGICAL :: limplicit                ! logical flag for implicit solvers (.TRUE.=>implicit, .FALSE.=>explicit)   
      REAL(dp) :: dt                      ! the time step by which the propagator increments
      TYPE(it_solver_type), POINTER :: &  ! the object that actually does linear algebra
          implicit_solver          

      CONTAINS
        PROCEDURE :: read_settings_file => read_propagator_settings
	PROCEDURE :: print_summary => print_propagator_summary
	PROCEDURE :: broadcast_settings => broadcast_propagator_settings
	PROCEDURE :: propagate => propagator_propagate

  END TYPE propagator_type

  CONTAINS

    SUBROUTINE read_propagator_settings(this)
        ! 
	! ... Reads in the namelist associated with the propagator
	! 
        IMPLICIT NONE
	! input variable
	CLASS(propagator_type), INTENT(INOUT) :: this 

        RETURN

    END SUBROUTINE read_propagator_settings

    SUBROUTINE print_propagator_summary(this, io_unit)
        ! 
	! ... Prints a summary of the propagator's settings
	! 
	IMPLICIT NONE
	! input variables
	CLASS(propagator_type), INTENT(INOUT) :: this
	INTEGER, INTENT(IN) :: io_unit

	WRITE(io_unit,'(5x,"Summary!")')

        RETURN

    END SUBROUTINE print_propagator_summary

#ifdef __MPI
    SUBROUTINE broadcast_propagator_settings(this)
        ! 
	! ... Broadcast input read in on IO node to all nodes
	!   
	USE mp,          ONLY : mp_bcast
	USE mp_world,    ONLY : world_comm 

	IMPLICIT NONE
	! input variable
	CLASS(propagator_type), INTENT(INOUT) :: this

        RETURN

    END SUBROUTINE broadcast_propagator_settings
#endif

    SUBROUTINE propagator_propagate(this, io_unit)
        ! 
        ! ... interface for marching the electrons forward by one step
        ! 
        IMPLICIT NONE
        ! input variables
        CLASS(propagator_type), INTENT(INOUT) :: this
	INTEGER :: io_unit

        WRITE(io_unit,'(6x,"step")')

        RETURN

    END SUBROUTINE propagator_propagate

END MODULE propagator_mod
