MODULE tddft_perturbations_mod
  ! 
  ! ... This module contains the variables used for various perturbations applied in TDDFT
  !
  USE kinds,			ONLY : dp
  USE tddft_envelope_mod

  IMPLICIT NONE
  ! 
  SAVE
  !
  TYPE projectile_perturbation_type

      CHARACTER(len=1) :: projectile_direction          ! direction of motion of the projectile, 'a', 'b', 'c', or 'r' = a lattice vector or a random vector
      INTEGER :: projectile_index                       ! the number of the atom that is designated as projectile
      REAL(dp) :: projectile_kinetic_energy             ! kinetic energy of the projectile atom in Rydberg units
      TYPE(tddft_envelope_type) :: projectile_envelope  ! defines the time-dependence of the motion of the projectile

      CONTAINS
        PROCEDURE :: read_settings_file => read_projectile_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_projectile_inputs
#endif		
	 
  END TYPE projectile_perturbation_type

  TYPE, PUBLIC :: scalar_perturbation_type

      REAL(dp) :: efield_strength(3)                    ! electric field strength in Rydberg units (energy / length)
      TYPE(tddft_envelope_type) :: scalar_envelope      ! defines the time-dependence of the scalar kick

      CONTAINS
        PROCEDURE :: read_settings_file => read_scalar_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_scalar_inputs
#endif __MPI		
      
  END TYPE scalar_perturbation_type

  TYPE vector_perturbation_type

      REAL(dp) :: afield_strength(3)                    ! vector potential strength in Rydberg units (energy time / length)
      TYPE(tddft_envelope_type) :: vector_envelope      ! defines the time-dependence of the vector kick

      CONTAINS
        PROCEDURE :: read_settings_file => read_vector_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_vector_inputs
#endif		
      
  END TYPE vector_perturbation_type

  TYPE xray_perturbation_type

      CHARACTER(len=1) :: cos_or_sin                    ! character indicating whether the cosine or sine part of exp(i*q*r) is used as a perturbation
      INTEGER :: xray_pert(3)                           ! integer indices corresponding to q in a basis of reciprocal lattice vectors
      REAL(dp) :: xray_strength
      TYPE(tddft_envelope_type) :: xray_envelope        ! defines the time-dependence of the x-ray kick
      
      CONTAINS 
        PROCEDURE :: read_settings_file => read_xray_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_xray_inputs
#endif		

  END TYPE xray_perturbation_type

  CONTAINS

    SUBROUTINE read_projectile_settings(this)

        IMPLICIT NONE
        ! input variables
        CLASS(projectile_perturbation_type), INTENT(INOUT) :: this 
        ! internal variables
        CHARACTER(len=1) :: projectile_direction
	INTEGER :: projectile_index
	REAL(dp) :: projectile_kinetic_energy
	TYPE(tddft_envelope_type) :: projectile_envelope
	
        INTEGER :: ierr

        NAMELIST /projectile/ projectile_direction, projectile_index, projectile_kinetic_energy

        ! set default values
	projectile_direction = 'r'
	projectile_index = 1
	projectile_kinetic_energy = 1.0_dp	! note: add conversion for units

        ! read values from file
        READ(5, projectile, err = 201, iostat = ierr)
201 CALL errore('read_settings_file', 'reading projectile namelist', ierr)

        ! load values from file into calling instance
        this%projectile_direction = projectile_direction
	this%projectile_index = projectile_index
	this%projectile_kinetic_energy = projectile_kinetic_energy
	this%projectile_envelope = projectile_envelope

    END SUBROUTINE read_projectile_settings

    SUBROUTINE read_scalar_settings(this)

        IMPLICIT NONE

        ! input variables
        CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        REAL(dp) :: efield_strength(3)
	TYPE(tddft_envelope_type) :: scalar_envelope
        INTEGER :: ierr

        NAMELIST /scalar/ efield_strength

        ! set default values
	efield_strength(:) = (/0.001_dp, 0.0_dp, 0.0_dp/)

        ! read values from file
        READ(5, scalar, err = 202, iostat = ierr)
202 CALL errore('read_settings_file', 'reading scalar namelist', ierr)
        
	! load values from the file into the calling instance
        this%efield_strength(:) = efield_strength(:)

    END SUBROUTINE read_scalar_settings

    SUBROUTINE read_vector_settings(this)

        IMPLICIT NONE
        ! input variables
        CLASS(vector_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        REAL(dp) :: afield_strength(3)
	TYPE(tddft_envelope_type) :: vector_envelope
        INTEGER :: ierr

        NAMELIST /vector/ afield_strength
   
        ! set default values
	afield_strength(:) = (/0.001_dp, 0.0_dp, 0.0_dp/)
   
        ! read values from file 
        READ(5, vector, err = 203, iostat = ierr)
203 CALL errore('read_settings_file', 'reading vector namelist', ierr)

        ! load values from file into the calling instance
        this%afield_strength(:) = afield_strength(:)

    END SUBROUTINE read_vector_settings

    SUBROUTINE read_xray_settings(this)

        IMPLICIT NONE 
        ! input variables
        CLASS(xray_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        CHARACTER(len=1) :: cos_or_sin
	INTEGER :: xray_pert(3)
	REAL(dp) :: xray_strength
	TYPE(tddft_envelope_type) :: xray_envelope
        INTEGER :: ierr

        NAMELIST /xray/ cos_or_sin, xray_pert, xray_strength

        ! set default values
        cos_or_sin = 'c'
	xray_pert = (/1, 0, 0/)
	xray_strength = 0.001_dp

        ! read values from file
        READ(5, xray, err = 204, iostat = ierr)
204 CALL errore('read_settings_file', 'reading xray namelist', ierr)

        ! load values from file into the calling instance
        this%cos_or_sin = cos_or_sin
	this%xray_pert(:) = xray_pert(:)
	this%xray_strength = xray_strength

    END SUBROUTINE read_xray_settings

#ifdef __MPI

    SUBROUTINE broadcast_projectile_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(projectile_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%projectile_direction, root, world_comm)
	CALL mp_bcast(this%projectile_index, root, world_comm)
	CALL mp_bcast(this%projectile_kinetic_energy, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_projectile_inputs

    SUBROUTINE broadcast_scalar_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%efield_strength, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_scalar_inputs

    SUBROUTINE broadcast_vector_inputs(this)
         
        USE mp,		ONLY : mp_bcast
	USE mp_world, 	ONLY : world_comm

	IMPLICIT NONE
	! input variables
	CLASS(vector_perturbation_type), INTENT(INOUT) :: this
	! internal variables
	INTEGER, parameter :: root = 0

        CALL mp_bcast(this%afield_strength, root, world_comm)

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

        CALL mp_bcast(this%cos_or_sin, root, world_comm)
	CALL mp_bcast(this%xray_pert, root, world_comm)
	CALL mp_bcast(this%xray_strength, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_xray_inputs
#endif

END MODULE tddft_perturbations_mod
