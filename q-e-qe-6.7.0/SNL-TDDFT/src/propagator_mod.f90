! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

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
 
    COMPLEX(dp) :: cn_factor                 ! scalar factor of idt/2 for Crank-Nicolson
    INTEGER :: nstages                       ! number of stages in the time integrator
    LOGICAL :: limplicit                     ! logical flag for implicit solvers (.TRUE.=>implicit, .FALSE.=>explicit)
    TYPE(it_solver_type) :: implicit_solver  ! the object that actually does linear algebra

  CONTAINS
    PROCEDURE :: read_settings_file => read_propagator_settings
    PROCEDURE :: print_summary => print_propagator_summary
    PROCEDURE :: broadcast_settings => broadcast_propagator_settings

  END TYPE propagator_type

CONTAINS

  SUBROUTINE read_propagator_settings(this)
    !
    ! ... Reads in the namelist associated with the propagator
    !
    IMPLICIT NONE
    ! input variable
    CLASS(propagator_type), INTENT(INOUT) :: this
    ! internal variables
    INTEGER :: ierr
    INTEGER :: nstages
    INTEGER :: max_kry_dim, max_restarts
    LOGICAL :: limplicit
    REAL(dp) :: tol

    NAMELIST /propagator/ limplicit, nstages, &
                          max_kry_dim, max_restarts, tol
 
    ! set default values
    limplicit = .TRUE.
    nstages = 2
 
    max_kry_dim = 50
    max_restarts = 10
    tol = 1.e-12_DP

    READ(5, propagator, err = 201, iostat = ierr)
    201 CALL errore('read_propagator_settings', 'reading propagator namelist', ierr)

    this%limplicit = limplicit
    this%nstages = nstages

    IF(this%limplicit)THEN
      this%implicit_solver%max_kry_dim = max_kry_dim
      this%implicit_solver%max_restarts = max_restarts
      this%implicit_solver%tol = tol
    ELSE
      CALL errore('read_propagator_settings', 'no non-implicit propagators...yet', ierr)
    ENDIF

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
    ! internal variable
    INTEGER :: ierr

    IF(this%limplicit)THEN
      WRITE(io_unit,'(5x,"Implicit propagator")')
      WRITE(io_unit,'(5x," Number of stages         =",I12)') this%nstages
      WRITE(io_unit,'(5x," Max Krylov dimension     =",I12)') this%implicit_solver%max_kry_dim
      WRITE(io_unit,'(5x," Max restarts             =",I12)') this%implicit_solver%max_restarts
      WRITE(io_unit,'(5x," Tolerance for Ax=b       =",E12.4)') this%implicit_solver%tol
    ELSE
      CALL errore('print_propagatory_summary', 'no non-implicit propagators...yet', ierr)
    ENDIF
    
    WRITE(io_unit,*)

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
    ! internal variables
    INTEGER, parameter :: root = 0

    CALL mp_bcast(this%limplicit, root, world_comm)
    CALL mp_bcast(this%nstages, root, world_comm)
    CALL mp_bcast(this%implicit_solver%max_kry_dim, root, world_comm)
    CALL mp_bcast(this%implicit_solver%max_restarts, root, world_comm)
    CALL mp_bcast(this%implicit_solver%tol, root, world_comm)

    RETURN

  END SUBROUTINE broadcast_propagator_settings
#endif

END MODULE propagator_mod
