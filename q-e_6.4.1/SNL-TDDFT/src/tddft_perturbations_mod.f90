MODULE tddft_perturbations_mod
  ! 
  ! ... This module contains the variables used for various perturbations applied in TDDFT
  !
  USE kinds,			ONLY : dp
  IMPLICIT NONE
  ! 
  SAVE
  !
  PUBLIC :: &
      read_settings_file_scalar,	&
      read_settings_file_stopping,      &
      read_settings_file_vector,        &
      read_settings_file_xray
  
  TYPE scalar_perturbation_type

      INTEGER :: dummy
      
  END TYPE scalar_perturbation_type

  TYPE stopping_perturbation_type

      INTEGER :: dummy
      
  END TYPE stopping_perturbation_type

  TYPE vector_perturbation_type

      INTEGER :: dummy
      
  END TYPE vector_perturbation_type

  TYPE xray_perturbation_type

      INTEGER :: dummy

  END TYPE xray_perturbation_type

  CONTAINS

    SUBROUTINE read_settings_file_scalar(this)

        IMPLICIT NONE

        ! input variables
        TYPE(scalar_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /scalar/ dummy

        READ(5, scalar, err = 201, iostat = ierr)
201 CALL errore('read_settings_file', 'reading scalar namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_settings_file_scalar

    SUBROUTINE read_settings_file_stopping(this)

        IMPLICIT NONE
        ! input variables
        TYPE(stopping_perturbation_type), INTENT(INOUT) :: this 
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /stopping/ dummy

        READ(5, stopping, err = 202, iostat = ierr)
202 CALL errore('read_settings_file', 'reading stopping namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_settings_file_stopping

    SUBROUTINE read_settings_file_vector(this)

        IMPLICIT NONE
        ! input variables
        TYPE(vector_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /vector/ dummy

        READ(5, vector, err = 203, iostat = ierr)
203 CALL errore('read_settings_file', 'reading vector namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_settings_file_vector

    SUBROUTINE read_settings_file_xray(this)

        IMPLICIT NONE 
        ! input variables
        TYPE(xray_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        INTEGER :: dummy
        INTEGER :: ierr

        NAMELIST /xray/ dummy

        READ(5, xray, err = 204, iostat = ierr)
204 CALL errore('read_settings_file', 'reading xray namelist', ierr)

        this%dummy = dummy

    END SUBROUTINE read_settings_file_xray

END MODULE tddft_perturbations_mod
