MODULE tddft_perturbations_mod
  ! 
  ! ... This module contains the variables used for various perturbations applied in TDDFT
  !
  USE kinds,			ONLY : dp
  IMPLICIT NONE
  ! 
  SAVE
  !
  TYPE scalar_perturbation_type

      INTEGER :: dummy
      PUBLIC :: read_settings_file
      
      CONTAINS

      SUBROUTINE read_settings_file

          IMPLICIT NONE

          INTEGER :: dummy
          INTEGER :: ierr

          READ(5, scalar, err = 201, iostat = ierr)
201 CALL errore('read_settings_file', 'reading scalar namelist', ierr)

          this%dummy = dummy

      END SUBROUTINE read_settings_file

  END TYPE scalar_perturbation_type
  
  TYPE stopping_perturbation_type

      INTEGER :: dummy
      PUBLIC :: read_settings_file
      
      CONTAINS

      SUBROUTINE read_settings_file

          IMPLICIT NONE

          INTEGER :: dummy
          INTEGER :: ierr

          READ(5, scalar, err = 202, iostat = ierr)
202 CALL errore('read_settings_file', 'reading stopping namelist', ierr)

          this%dummy = dummy

      END SUBROUTINE read_settings_file

  END TYPE stopping_perturbation_type
  
  TYPE vector_perturbation_type

      INTEGER :: dummy
      PUBLIC :: read_settings_file
      
      CONTAINS

      SUBROUTINE read_settings_file

          IMPLICIT NONE

          INTEGER :: dummy
          INTEGER :: ierr

          READ(5, scalar, err = 203, iostat = ierr)
203 CALL errore('read_settings_file', 'reading vector namelist', ierr)

          this%dummy = dummy

      END SUBROUTINE read_settings_file

  END TYPE vector_perturbation_type
    
  TYPE xray_perturbation_type

      INTEGER :: dummy
      PUBLIC :: read_settings_file

      CONTAINS

      SUBROUTINE read_settings_file

          IMPLICIT NONE

          INTEGER :: dummy
          INTEGER :: ierr

          READ(5, scalar, err = 204, iostat = ierr)
204 CALL errore('read_settings_file', 'reading xray namelist', ierr)

          this%dummy = dummy

      END SUBROUTINE read_settings_file
 
  END TYPE xray_perturbation_type 

END MODULE tddft_perturbations_mod
