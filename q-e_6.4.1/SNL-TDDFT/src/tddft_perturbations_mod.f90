MODULE tddft_perturbations_mod
  ! 
  ! ... This module contains the variables used for various perturbations applied in TDDFT
  !
  USE kinds,			ONLY : dp
  IMPLICIT NONE
  ! 
  SAVE
  !
  TYPE, PUBLIC :: scalar_perturbation_type

      INTEGER :: dummy
      CONTAINS
        PROCEDURE :: read_settings_file => read_scalar_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_scalar_inputs
#endif __MPI		
      
  END TYPE scalar_perturbation_type

  TYPE stopping_perturbation_type

      INTEGER :: dummy
      CONTAINS
        PROCEDURE :: read_settings_file => read_stopping_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_stopping_inputs
#endif		
	 
  END TYPE stopping_perturbation_type

  TYPE vector_perturbation_type

      INTEGER :: dummy
      CONTAINS
        PROCEDURE :: read_settings_file => read_vector_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_vector_inputs
#endif		
      
  END TYPE vector_perturbation_type

  TYPE xray_perturbation_type

      INTEGER :: dummy
      CONTAINS 
        PROCEDURE :: read_settings_file => read_xray_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_xray_inputs
#endif		

  END TYPE xray_perturbation_type

  CONTAINS

    SUBROUTINE read_scalar_settings(this)

        IMPLICIT NONE

        ! input variables
        CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /scalar/ dummy

        READ(5, scalar, err = 201, iostat = ierr)
201 CALL errore('read_settings_file', 'reading scalar namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_scalar_settings

    SUBROUTINE read_stopping_settings(this)

        IMPLICIT NONE
        ! input variables
        CLASS(stopping_perturbation_type), INTENT(INOUT) :: this 
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /stopping/ dummy

        READ(5, stopping, err = 202, iostat = ierr)
202 CALL errore('read_settings_file', 'reading stopping namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_stopping_settings

    SUBROUTINE read_vector_settings(this)

        IMPLICIT NONE
        ! input variables
        CLASS(vector_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /vector/ dummy

        READ(5, vector, err = 203, iostat = ierr)
203 CALL errore('read_settings_file', 'reading vector namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_vector_settings

    SUBROUTINE read_xray_settings(this)

        IMPLICIT NONE 
        ! input variables
        CLASS(xray_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /xray/ dummy

        READ(5, xray, err = 204, iostat = ierr)
204 CALL errore('read_settings_file', 'reading xray namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_xray_settings

#ifdef __MPI
    SUBROUTINE broadcast_scalar_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%dummy, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_scalar_inputs

    SUBROUTINE broadcast_stopping_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(stopping_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%dummy, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_stopping_inputs

    SUBROUTINE broadcast_vector_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(vector_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%dummy, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_vector_inputs

    SUBROUTINE broadcast_xray_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(xray_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%dummy, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_xray_inputs
#endif

END MODULE tddft_perturbations_mod
