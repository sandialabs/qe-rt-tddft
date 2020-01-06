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
      REAL(dp) :: projectile_kinetic_energy             ! kinetic energy of the projectile in units of eV
      TYPE(tddft_envelope_type) :: projectile_envelope  ! defines the time-dependence of the motion of the projectile

      CONTAINS
        PROCEDURE :: print_summary => print_projectile_summary
        PROCEDURE :: read_settings_file => read_projectile_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_projectile_inputs
#endif		
	 
  END TYPE projectile_perturbation_type

  TYPE, PUBLIC :: scalar_perturbation_type

      REAL(dp) :: efield_strength(3)                    ! strength in units of eV/Angstrom
      TYPE(tddft_envelope_type) :: scalar_envelope      ! defines the time-dependence of the scalar kick

      CONTAINS
        PROCEDURE :: print_summary => print_scalar_summary
        PROCEDURE :: read_settings_file => read_scalar_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_scalar_inputs
#endif __MPI		
      
  END TYPE scalar_perturbation_type

  TYPE vector_perturbation_type

      REAL(dp) :: afield_strength(3)                    ! strength in units of eV/Angstrom
      TYPE(tddft_envelope_type) :: vector_envelope      ! defines the time-dependence of the vector kick

      CONTAINS
        PROCEDURE :: print_summary => print_vector_summary
        PROCEDURE :: read_settings_file => read_vector_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_vector_inputs
#endif		
      
  END TYPE vector_perturbation_type

  TYPE xray_perturbation_type

      CHARACTER(len=1) :: cos_or_sin                    ! character indicating whether the cosine or sine part of exp(i*q*r) is used as a perturbation
      INTEGER :: xray_pert(3)                           ! integer indices corresponding to q in a basis of reciprocal lattice vectors
      REAL(dp) :: xray_strength                         ! strength in units of eV
      TYPE(tddft_envelope_type) :: xray_envelope        ! defines the time-dependence of the x-ray kick
      
      CONTAINS 
        PROCEDURE :: print_summary => print_xray_summary
        PROCEDURE :: read_settings_file => read_xray_settings
#ifdef __MPI
	PROCEDURE :: broadcast_inputs => broadcast_xray_inputs
#endif		

  END TYPE xray_perturbation_type

  CONTAINS

    SUBROUTINE print_projectile_summary(this, io_unit) 
        !
	! ... Prints a summary of the settings in this instance of projectile perturbation
        !
        IMPLICIT NONE
	! input variables
	CLASS(projectile_perturbation_type), INTENT(INOUT) :: this
	INTEGER, INTENT(IN) :: io_unit

	WRITE(io_unit,'(5x,"Projectile perturbation active")')
	WRITE(io_unit,'(5x,"Direction                  =",A)') this%projectile_direction
	WRITE(io_unit,'(5x,"Index                      =",I12)') this%projectile_index
	WRITE(io_unit,'(5x,"Kinetic energy             =",F12.4, " eV ")') this%projectile_kinetic_energy
	WRITE(io_unit,'(5x,"Envelope")')
	CALL this%projectile_envelope%print_summary(io_unit)

        RETURN 

    END SUBROUTINE print_projectile_summary

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

	! envelope variables
	INTEGER :: envelope_index
	REAL(dp) :: linear_slope, amplitude, carrier_frequency, delay, width
	
        NAMELIST /projectile/ projectile_direction, projectile_index, projectile_kinetic_energy, &
                              envelope_index, linear_slope, amplitude, carrier_frequency, delay, width

        ! set default values
	projectile_direction = 'r'
	projectile_index = 1

        ! default values for the envelope variables
	envelope_index = 0
	linear_slope = 0.0_dp
	amplitude = 0.0_dp
	carrier_frequency = 0.0_dp
	delay = 0.0_dp
	width = 0.0_dp
	
        ! read values from file
        READ(5, projectile, err = 201, iostat = ierr)
201 CALL errore('read_settings_file', 'reading projectile namelist', ierr)

        ! load values from file into calling instance
        this%projectile_direction = projectile_direction
	this%projectile_index = projectile_index
	this%projectile_kinetic_energy = projectile_kinetic_energy
	! envelope values
	this%projectile_envelope%envelope_index = envelope_index
	this%projectile_envelope%linear_slope = linear_slope
        this%projectile_envelope%amplitude = amplitude
	this%projectile_envelope%carrier_frequency = carrier_frequency
	this%projectile_envelope%delay = delay
	this%projectile_envelope%width = width

    END SUBROUTINE read_projectile_settings

    SUBROUTINE print_scalar_summary(this, io_unit) 
        !
	! ... Prints a summary of the settings in this instance of scalar perturbation
        !
        IMPLICIT NONE
	! input variables
	CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
	INTEGER, INTENT(IN) :: io_unit

	WRITE(io_unit,'(5x,"Scalar perturbation active")')
	WRITE(io_unit,'(5x,"E-field strength           =",3F12.4," eV/A ")') this%efield_strength(1:3)
	WRITE(io_unit,'(5x,"Envelope")')
	CALL this%scalar_envelope%print_summary(io_unit)

        RETURN 

    END SUBROUTINE print_scalar_summary

    SUBROUTINE read_scalar_settings(this)

        IMPLICIT NONE

        ! input variables
        CLASS(scalar_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        REAL(dp) :: efield_strength(3)
	TYPE(tddft_envelope_type) :: scalar_envelope
        INTEGER :: ierr

	! envelope variables
	INTEGER :: envelope_index
	REAL(dp) :: linear_slope, amplitude, carrier_frequency, delay, width

        NAMELIST /scalar/ efield_strength, &
                          envelope_index, linear_slope, amplitude, carrier_frequency, delay, width

        ! set default values
	efield_strength(:) = (/0.001_dp, 0.0_dp, 0.0_dp/)

        ! default values for the envelope variables
	envelope_index = 0
	linear_slope = 0.0_dp
	amplitude = 0.0_dp
	carrier_frequency = 0.0_dp
	delay = 0.0_dp
	width = 0.0_dp

        ! read values from file
        READ(5, scalar, err = 202, iostat = ierr)
202 CALL errore('read_settings_file', 'reading scalar namelist', ierr)
        
	! load values from the file into the calling instance
        this%efield_strength(:) = efield_strength(:)
	! envelope values
	this%scalar_envelope%envelope_index = envelope_index
	this%scalar_envelope%linear_slope = linear_slope
        this%scalar_envelope%amplitude = amplitude
	this%scalar_envelope%carrier_frequency = carrier_frequency
	this%scalar_envelope%delay = delay
	this%scalar_envelope%width = width

    END SUBROUTINE read_scalar_settings

    SUBROUTINE print_vector_summary(this, io_unit) 
        !
	! ... Prints a summary of the settings in this instance of vector perturbation
        !
        IMPLICIT NONE
	! input variables
	CLASS(vector_perturbation_type), INTENT(INOUT) :: this
	INTEGER, INTENT(IN) :: io_unit

	WRITE(io_unit,'(5x,"Vector perturbation active")')

        RETURN 

    END SUBROUTINE print_vector_summary

    SUBROUTINE read_vector_settings(this)

        IMPLICIT NONE
        ! input variables
        CLASS(vector_perturbation_type), INTENT(INOUT) :: this
        ! internal variables
        REAL(dp) :: afield_strength(3)
	TYPE(tddft_envelope_type) :: vector_envelope
        INTEGER :: ierr

	! envelope variables
	INTEGER :: envelope_index
	REAL(dp) :: linear_slope, amplitude, carrier_frequency, delay, width

        NAMELIST /vector/ afield_strength, &
                          envelope_index, linear_slope, amplitude, carrier_frequency, delay, width
   
        ! set default values
	afield_strength(:) = (/0.001_dp, 0.0_dp, 0.0_dp/)

        ! default values for the envelope variables
	envelope_index = 0
	linear_slope = 0.0_dp
	amplitude = 0.0_dp
	carrier_frequency = 0.0_dp
	delay = 0.0_dp
	width = 0.0_dp
	   
        ! read values from file 
        READ(5, vector, err = 203, iostat = ierr)
203 CALL errore('read_settings_file', 'reading vector namelist', ierr)

        ! load values from file into the calling instance
        this%afield_strength(:) = afield_strength(:)
	! envelope values
	this%vector_envelope%envelope_index = envelope_index
	this%vector_envelope%linear_slope = linear_slope
        this%vector_envelope%amplitude = amplitude
	this%vector_envelope%carrier_frequency = carrier_frequency
	this%vector_envelope%delay = delay
	this%vector_envelope%width = width

    END SUBROUTINE read_vector_settings

    SUBROUTINE print_xray_summary(this, io_unit) 
        !
	! ... Prints a summary of the settings in this instance of xray perturbation
        !
        IMPLICIT NONE
	! input variables
	CLASS(xray_perturbation_type), INTENT(INOUT) :: this
	INTEGER, INTENT(IN) :: io_unit

	WRITE(io_unit,'(5x,"X-ray perturbation active")')

        RETURN 

    END SUBROUTINE print_xray_summary

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

	! envelope variables
	INTEGER :: envelope_index
	REAL(dp) :: linear_slope, amplitude, carrier_frequency, delay, width

        NAMELIST /xray/ cos_or_sin, xray_pert, xray_strength, &
                        envelope_index, linear_slope, amplitude, carrier_frequency, delay, width

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
	! envelope values
	this%xray_envelope%envelope_index = envelope_index
	this%xray_envelope%linear_slope = linear_slope
        this%xray_envelope%amplitude = amplitude
	this%xray_envelope%carrier_frequency = carrier_frequency
	this%xray_envelope%delay = delay
	this%xray_envelope%width = width

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
	CALL mp_bcast(this%projectile_envelope%envelope_index, root, world_comm)
	CALL mp_bcast(this%projectile_envelope%linear_slope, root, world_comm)
	CALL mp_bcast(this%projectile_envelope%amplitude, root, world_comm)
	CALL mp_bcast(this%projectile_envelope%carrier_frequency, root, world_comm)
	CALL mp_bcast(this%projectile_envelope%delay, root, world_comm)
	CALL mp_bcast(this%projectile_envelope%width, root, world_comm)

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
	CALL mp_bcast(this%scalar_envelope%envelope_index, root, world_comm)
	CALL mp_bcast(this%scalar_envelope%linear_slope, root, world_comm)
	CALL mp_bcast(this%scalar_envelope%amplitude, root, world_comm)
	CALL mp_bcast(this%scalar_envelope%carrier_frequency, root, world_comm)
	CALL mp_bcast(this%scalar_envelope%delay, root, world_comm)
	CALL mp_bcast(this%scalar_envelope%width, root, world_comm)

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
	CALL mp_bcast(this%vector_envelope%envelope_index, root, world_comm)
	CALL mp_bcast(this%vector_envelope%linear_slope, root, world_comm)
	CALL mp_bcast(this%vector_envelope%amplitude, root, world_comm)
	CALL mp_bcast(this%vector_envelope%carrier_frequency, root, world_comm)
	CALL mp_bcast(this%vector_envelope%delay, root, world_comm)
	CALL mp_bcast(this%vector_envelope%width, root, world_comm)

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
	CALL mp_bcast(this%xray_envelope%envelope_index, root, world_comm)
	CALL mp_bcast(this%xray_envelope%linear_slope, root, world_comm)
	CALL mp_bcast(this%xray_envelope%amplitude, root, world_comm)
	CALL mp_bcast(this%xray_envelope%carrier_frequency, root, world_comm)
	CALL mp_bcast(this%xray_envelope%delay, root, world_comm)
	CALL mp_bcast(this%xray_envelope%width, root, world_comm)

	RETURN

    END SUBROUTINE broadcast_xray_inputs
#endif

END MODULE tddft_perturbations_mod
