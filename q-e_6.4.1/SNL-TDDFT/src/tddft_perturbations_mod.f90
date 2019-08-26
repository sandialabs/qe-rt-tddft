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
      
  END TYPE scalar_perturbation_type

  TYPE stopping_perturbation_type

      INTEGER :: dummy
      CONTAINS
        PROCEDURE :: read_settings_file => read_stopping_settings
	 
  END TYPE stopping_perturbation_type

  TYPE vector_perturbation_type

      INTEGER :: dummy
      CONTAINS
        PROCEDURE :: read_settings_file => read_vector_settings
      
  END TYPE vector_perturbation_type

  TYPE xray_perturbation_type

      INTEGER :: dummy
      CONTAINS 
        PROCEDURE :: read_settings_file => read_xray_settings

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

END MODULE tddft_perturbations_mod
