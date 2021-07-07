! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

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

    INTEGER :: projectile_index                       ! the number of the atom that is designated as projectile  
    REAL(dp) :: projectile_kinetic_energy             ! kinetic energy of the projectile in units of eV
    REAL(dp) :: projectile_velocity(3)         ! vector that orients the projectile's velocity (normalized...or forced to be normalized :))
    TYPE(tddft_envelope_type) :: projectile_envelope  ! defines the time-dependence of the motion of the projectile

  CONTAINS
    PROCEDURE :: print_summary => print_projectile_summary
    PROCEDURE :: read_settings_file => read_projectile_settings
#ifdef __MPI
    PROCEDURE :: broadcast_settings => broadcast_projectile_settings
#endif

  END TYPE projectile_perturbation_type

  TYPE, PUBLIC :: scalar_perturbation_type

    REAL(dp) :: efield_strength(3)                    ! strength in units of eV/Angstrom
    TYPE(tddft_envelope_type) :: scalar_envelope      ! defines the time-dependence of the scalar kick

  CONTAINS
    PROCEDURE :: print_summary => print_scalar_summary
    PROCEDURE :: read_settings_file => read_scalar_settings
#ifdef __MPI
    PROCEDURE :: broadcast_settings => broadcast_scalar_settings
#endif __MPI

  END TYPE scalar_perturbation_type

  TYPE vector_perturbation_type

    REAL(dp) :: afield_strength(3)                    ! strength in units of eV/Angstrom
    TYPE(tddft_envelope_type) :: vector_envelope      ! defines the time-dependence of the vector kick

  CONTAINS
    PROCEDURE :: print_summary => print_vector_summary
    PROCEDURE :: read_settings_file => read_vector_settings
#ifdef __MPI
    PROCEDURE :: broadcast_settings => broadcast_vector_settings
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
    PROCEDURE :: broadcast_settings => broadcast_xray_settings
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
    WRITE(io_unit,'(5x,"Velocity vector            =",3F12.4," unitless ")') this%projectile_velocity(1:3)
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
    INTEGER :: projectile_index
    REAL(dp) :: projectile_kinetic_energy
    REAL(dp) :: projectile_velocity(3)
    REAL(dp) :: scratch_real
    TYPE(tddft_envelope_type) :: projectile_envelope
    INTEGER :: ierr

    ! envelope variables
    INTEGER :: envelope_index
    REAL(dp) :: linear_slope, amplitude, carrier_frequency, delay, width

    NAMELIST /projectile/ projectile_velocity, projectile_index, projectile_kinetic_energy, &
    envelope_index, linear_slope, amplitude, carrier_frequency, delay, width

    ! set default values
    projectile_velocity = (/0.0_dp, 0.0_dp, 1.0_dp/)
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
    201 CALL errore('read_projectile_settings', 'reading projectile namelist', ierr)

    ! force normalization of the projectile velocity vector
    scratch_real = SQRT(projectile_velocity(1)**2 + projectile_velocity(2)**2 + projectile_velocity(3)**2)

    ! load values from file into calling instance
    this%projectile_velocity(:) = projectile_velocity(:)/scratch_real
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
202 CALL errore('read_scalar_settings', 'reading scalar namelist', ierr)

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
203 CALL errore('read_vector_settings', 'reading vector namelist', ierr)

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
204 CALL errore('read_xray_settings', 'reading xray namelist', ierr)

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

  SUBROUTINE broadcast_projectile_settings(this)

    USE mp,		ONLY : mp_bcast
    USE mp_world, 	ONLY : world_comm

    IMPLICIT NONE
    ! input variables
    CLASS(projectile_perturbation_type), INTENT(INOUT) :: this
    ! internal variables
    INTEGER, parameter :: root = 0

    CALL mp_bcast(this%projectile_velocity, root, world_comm)
    CALL mp_bcast(this%projectile_index, root, world_comm)
    CALL mp_bcast(this%projectile_kinetic_energy, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%envelope_index, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%linear_slope, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%amplitude, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%carrier_frequency, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%delay, root, world_comm)
    CALL mp_bcast(this%projectile_envelope%width, root, world_comm)

    RETURN

  END SUBROUTINE broadcast_projectile_settings

  SUBROUTINE broadcast_scalar_settings(this)

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

  END SUBROUTINE broadcast_scalar_settings

  SUBROUTINE broadcast_vector_settings(this)

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

  END SUBROUTINE broadcast_vector_settings

  SUBROUTINE broadcast_xray_settings(this)

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

  END SUBROUTINE broadcast_xray_settings
#endif

END MODULE tddft_perturbations_mod
