! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

MODULE tddft_mod
  !
  ! ... This module contains the variables used for TDDFT calculations
  !
  USE kinds,                   ONLY : dp
  USE propagator_mod
  USE tddft_perturbations_mod
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE, PUBLIC :: tddft_type

    INTEGER :: &
    iuntdorbs = 42,           &  ! fixed unit of the file associated with old time-dependent orbitals
    iverbosity,		    &  ! integer indicating the level of verbosity
    nbands_occupied_max,		& ! integer indicating the largest number of occupied bands at any given k point
    nsteps_el,          	    &  ! total number of electronic steps
    nsteps_el_per_nsteps_ion, &  ! number of electronic steps per ionic steps
    nsteps_ion		       ! total number of ionic steps
    LOGICAL :: &
    lcorrect_ehrenfest_forces,&  ! flag = .TRUE. => compute the correct Ehrenfest forces (USPP/PAW feature)
    lcorrect_moving_ions,     &  ! flag = .TRUE. => compute the "gauge" correction for moving ions (USPP/PAW feature)
    lprojectile_perturbation, &  ! flag = .TRUE. => applies a perturbation by moving an ion at a fixed velocity
    lscalar_perturbation,     &  ! flag = .TRUE. => applies a perturbation through a homogeneous scalar potential
    lvector_perturbation,     &  ! flag = .TRUE. => applies a perturbation through a homogeneous vector potential
    lxray_perturbation           ! flag = .TRUE. => applies a perturbation through an inhomogeneous scalar potential
    REAL(dp) :: &
    band_threshold,           &  ! threshold for considering a band as occupied
    dt_el,                    &  ! electronic time step
    dt_ion,                   &  ! ionic time step
    duration                     ! total duration in attoseconds
    TYPE(projectile_perturbation_type) :: projectile_perturbation
    TYPE(scalar_perturbation_type) :: scalar_perturbation
    TYPE(vector_perturbation_type) :: vector_perturbation
    TYPE(xray_perturbation_type) :: xray_perturbation
    TYPE(propagator_type) :: propagator

    INTEGER, ALLOCATABLE :: nbands_occupied(:)	! array consisting of the number of occupied bands at each k point
    COMPLEX(dp), ALLOCATABLE :: psi(:,:,:)  ! orbitals across all stages of the propagator
    COMPLEX(dp), ALLOCATABLE :: Hpsi(:,:), Spsi(:,:)  ! H and S acting on the orbitals
    COMPLEX(dp), ALLOCATABLE :: rhs(:,:)  ! RHS for implicit propagators

  CONTAINS
    PROCEDURE :: read_settings_file => read_tddft_settings  ! reads the settings file for a TDDFT calculation in from a file
    PROCEDURE :: print_summary_clock => print_tddft_summary_clock  ! prints out clock information about this calculation on stdout
    PROCEDURE :: print_summary => print_tddft_summary  ! prints out information about this calculation on stdout    
    PROCEDURE :: perfunctory_business => perfunctory_tddft_business ! sets things up so the rest of the code doesn't throw a fit
    PROCEDURE :: initialize_calculation => initialize_tddft_calculation ! initializes various quantities before heading into the main loop
#ifdef __MPI
    PROCEDURE :: broadcast_settings => broadcast_tddft_settings  ! broadcasts settings to all tasks after reading in settings on the IO node
    PROCEDURE :: stop_calculation => stop_tddft_calculation  ! synchronizes processes before stopping
#endif
    PROCEDURE :: open_files => open_tddft_files ! opens file containing the Kohn-Sham orbitals
    PROCEDURE :: close_files => close_tddft_files ! closes file containing the Kohn-Sham orbitals

    PROCEDURE :: set_hamiltonian => set_tddft_hamiltonian ! sets the Hamiltonian to its value at the passed time step
    PROCEDURE :: allocate_preloop => allocate_tddft_preloop  ! allocates space for orbitals and stuff before the main loop
    PROCEDURE :: deallocate_postloop => deallocate_tddft_postloop ! deallocates things after the main loop

    PROCEDURE :: propagate_all_orbitals => propagate_tddft_all_orbitals ! updates all of the orbitals in the calculation...

  END TYPE tddft_type

CONTAINS

  SUBROUTINE read_tddft_settings(this)
    !
    ! ... Reads in the tddft input file, which is just a single namelist for now
    !
    USE io_files,         ONLY : prefix, tmp_dir
    USE io_global,        ONLY : ionode
    USE mp_images,        ONLY : my_image_id

    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this

    INTEGER :: ierr
    CHARACTER(len=256), EXTERNAL :: trimcheck
    CHARACTER(len=256) :: verbosity
    INTEGER :: iverbosity, nsteps_ion, nsteps_el, nsteps_el_per_nsteps_ion
    LOGICAL :: &
    lcorrect_ehrenfest_forces,&
    lcorrect_moving_ions,     &
    lprojectile_perturbation,   &
    lscalar_perturbation,     &
    lvector_perturbation,     &
    lxray_perturbation
    REAL(dp) :: band_threshold, dt_el, dt_ion, duration

    NAMELIST /tddft/ prefix, tmp_dir, verbosity, &
    nsteps_el, nsteps_ion, &
    band_threshold, &
    dt_el, dt_ion, duration, &
    lcorrect_ehrenfest_forces, lcorrect_moving_ions, &
    lprojectile_perturbation, lscalar_perturbation, &
    lvector_perturbation, lxray_perturbation

    IF(ionode .or. my_image_id == 0)THEN

      ! attaches unit 5 to the file specified as a command line argument
      CALL input_from_file()

      ! define input default values
      CALL get_environment_variable( 'ESPRESSO_TMPDIR', tmp_dir )
      IF(trim(tmp_dir) == ' ') tmp_dir = './scratch/'
      tmp_dir = trimcheck(tmp_dir)
      prefix       = 'pwscf'
      tmp_dir      = './scratch/'
      verbosity    = 'low'

      lcorrect_ehrenfest_forces = .FALSE.
      lcorrect_moving_ions = .FALSE.
      lprojectile_perturbation = .FALSE.
      lscalar_perturbation = .FALSE.
      lvector_perturbation = .FALSE.
      lxray_perturbation = .FALSE.

      ! default band threshold is the usual rule from Mike Desjarlais
      band_threshold = 1.0e-5_dp

      ! default propagation is one step that is 1 as long
      dt_el = 0.0_dp
      dt_ion = 0.0_dp
      duration = 1.0_dp
      nsteps_el = 1
      nsteps_ion = 1
      nsteps_el_per_nsteps_ion = 1

      ! read tddft namelist from the input file
      READ( 5, tddft, err = 200, iostat = ierr )
      200 CALL errore('read_tddft_settings', 'reading tddft namelist', ierr)

      this%lcorrect_ehrenfest_forces = lcorrect_ehrenfest_forces
      this%lcorrect_moving_ions = lcorrect_moving_ions
      this%lprojectile_perturbation = lprojectile_perturbation
      this%lscalar_perturbation = lscalar_perturbation
      this%lvector_perturbation = lvector_perturbation
      this%lxray_perturbation = lxray_perturbation

      this%band_threshold = band_threshold

      ! simulation duration is set by the specified duration
      this%duration = duration
      IF(dt_ion == 0.0_dp)THEN  ! if the ionic time step isn't set, then set the number of steps...
        this%nsteps_ion = nsteps_ion
        this%dt_ion = this%duration/this%nsteps_ion  ! ...and use that to compute the ionic time step
      ELSE  ! otherwise set the ionic time step and...
        this%dt_ion = dt_ion
        this%nsteps_ion = CEILING(this%duration/this%dt_ion)  ! ...use it to compute the number of ionic time steps
      ENDIF

      IF(dt_el == 0.0_dp)THEN  ! if the electronic time step isn't set, then set the number of steps...
        this%nsteps_el = nsteps_el
        this%dt_el = this%duration/this%nsteps_el  ! ...and use it to compute the number of electronic time steps
      ELSE  ! otherwise set the electronic time step and...
        this%dt_el = dt_el
        this%nsteps_el = CEILING(this%duration/this%dt_el)  ! ...use it to compute the number of electronic time steps
      ENDIF

      this%nsteps_el_per_nsteps_ion = nsteps_el/nsteps_ion

      ! set the integer verbosity flag
      SELECT CASE (verbosity)
      CASE('minimal')
        this%iverbosity = -1
      CASE('low')
        this%iverbosity = 0
      CASE('medium')
        this%iverbosity = 1
      CASE('high')
        this%iverbosity = 2
      CASE('debug')
        this%iverbosity = 3
      CASE DEFAULT
        CALL errore('read_tddft_settings', 'verbosity can be ''minimal'', ''low'', ''medium'', ''high'', or ''debug''', 1)
      END SELECT

      IF(this%lprojectile_perturbation)THEN
        CALL this%projectile_perturbation%read_settings_file()
      ENDIF

      IF(this%lscalar_perturbation)THEN
        CALL this%scalar_perturbation%read_settings_file()
      ENDIF

      IF(this%lvector_perturbation)THEN
        CALL this%vector_perturbation%read_settings_file()
      ENDIF

      IF(this%lxray_perturbation)THEN
        CALL this%xray_perturbation%read_settings_file()
      ENDIF

      CALL this%propagator%read_settings_file()

    ENDIF

#ifdef __MPI
    ! broadcast input variables
    CALL this%broadcast_settings
#endif

    RETURN

  END SUBROUTINE read_tddft_settings

  SUBROUTINE print_tddft_summary_clock(this, io_unit)
    ! 
    ! ... Prints the timing information for TDDFT routines called by this 
    !     instance to the io_unit passed as an argument
    ! 
    IMPLICIT NONE
    ! input variables
    CLASS(tddft_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit

    WRITE(io_unit,'(5x,"Summary of TDDFT timing information")')
    WRITE(io_unit,'(5x,"-----------------------------------")')
    WRITE(io_unit,'(5x,"Initialization")')
    CALL print_clock('init_calc')
    WRITE(io_unit,'(5x,"Time propagation")')
    CALL print_clock('set_hamiltonian') 
    RETURN
    
  END SUBROUTINE print_tddft_summary_clock

  SUBROUTINE print_tddft_summary(this, io_unit)
    !
    ! ... Prints a summary of the settings in this instance to the io_unit
    !     passed as an argument
    !
    IMPLICIT NONE
    ! input variables
    CLASS(tddft_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit

    WRITE(io_unit,'(5x,"Summary of TDDFT calculation")')
    WRITE(io_unit,'(5x,"----------------------------")')
    WRITE(io_unit,'(5x,"Duration                   = ",F12.4," as ",/)') this%duration

    WRITE(io_unit,'(5x,"Electronic time step       = ",F12.4," as ")') this%dt_el
    WRITE(io_unit,'(5x,"Total electronic steps     = ",I12,/)') this%nsteps_el

    WRITE(io_unit,'(5x,"Ionic time step            = ",F12.4," as ")') this%dt_ion
    WRITE(io_unit,'(5x,"Total ionic steps          = ",I12)') this%nsteps_ion
    WRITE(io_unit,'(5x,"Electron-ion step ratio    = ",I12,/)') this%nsteps_el_per_nsteps_ion

    IF ( this%lprojectile_perturbation ) THEN
      CALL this%projectile_perturbation%print_summary(io_unit)
    ENDIF

    IF ( this%lscalar_perturbation ) THEN
      CALL this%scalar_perturbation%print_summary(io_unit)
    ENDIF

    IF ( this%lvector_perturbation ) THEN
      CALL this%vector_perturbation%print_summary(io_unit)
    ENDIF

    IF ( this%lxray_perturbation ) THEN
      CALL this%xray_perturbation%print_summary(io_unit)
    ENDIF

    CALL this%propagator%print_summary(io_unit)

    RETURN

  END SUBROUTINE print_tddft_summary

  SUBROUTINE perfunctory_tddft_business(this)
    !
    ! ... Conducts business that might not be needed for TDDFT, per se, but *is* needed for the rest of QE to be happy
    ! 
    USE klist,  ONLY : nkstot
    USE wvfct,  ONLY : btype, nbndx
     
    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this
    ! internal variables
    INTEGER :: ierr   ! error flag  

    ! btype is allocated because sum_band needs it...
    ! really, it is used in diagonalization routines for determining which bands need to be fully converged 
    ! but sum_bands will be mad at us if this isn't allocated and we need it for constructing charge densities
    ALLOCATE(btype(nbndx, nkstot), stat=ierr)
    IF(ierr/=0) CALL errore('perfunctory_tddft_business','error allocating btype',ierr)

  END SUBROUTINE perfunctory_tddft_business

  SUBROUTINE initialize_tddft_calculation(this, io_unit)
    ! 
    ! ...Initializes all of the necessary (non-perfunctory) quantities before entering the main loop
    ! 
    USE constants,      ONLY : rytoev
    USE dfunct,         ONLY : newd
    USE fft_base,       ONLY : dfftp
    USE gvecs,          ONLY : doublegrid
    USE klist,          ONLY : nks, wk
    USE lsda_mod,       ONLY : nspin
    USE scf,            ONLY : v, vrs, vltot, kedtau
    USE wvfct,          ONLY : nbnd, et, wg  

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit
    ! internal
    INTEGER :: band_counter, kpt_counter
    INTEGER :: ierr   ! error flag

    ! start the clock on initialization
    CALL start_clock('init_calc')

    ! initialize pseudopotentials/projectors
    CALL init_us_1
    CALL init_at_1
    ! set the Hamiltonian at the 0th time step... 
    CALL this%set_hamiltonian(io_unit, 0, 0)

    ! Davide Ceresoli's CE-TDDFT code issues a bunch of warnings at this stage...might be good, someday

    ! allocate an array for tracking the number of occupied bands
    ALLOCATE(this%nbands_occupied(nks), stat=ierr)
    IF(ierr/=0) CALL errore('perfunctory_tddft_business','error allocating btype',ierr)
    this%nbands_occupied(:) = 0
    this%nbands_occupied_max = 0
    ! loop over each k point and evaluate the number of occupied bands
    DO kpt_counter = 1, nks
       ! loop over each band
       DO band_counter = 1, nbnd
          ! make sure that the k point that we're on is actually contributing to the density...
          IF(wk(kpt_counter) > 0.d0)THEN
             ! check the Fermi-Dirac weight, or at least something like it...
             IF(wg(band_counter, kpt_counter)/wk(kpt_counter) > this%band_threshold) this%nbands_occupied(kpt_counter) = band_counter
          ENDIF
       ENDDO
       ! if this is the largest number of occupied bands, so far, then update the maximum number of occupied bands
       IF(this%nbands_occupied(kpt_counter) > this%nbands_occupied_max) this%nbands_occupied_max = this%nbands_occupied(kpt_counter)
    ENDDO

    ! Davide Ceresoli's CE-TDDFT code computes alpha_pv here, used for shifting bands... I don't think we need it  
    
    ! stop the clock on initialization
    CALL stop_clock('init_calc')

  END SUBROUTINE initialize_tddft_calculation

#ifdef __MPI
  SUBROUTINE broadcast_tddft_settings(this)
    !
    ! ... Broadcast input read in on IO node to all nodes
    !
    USE io_files,	ONLY : prefix, tmp_dir
    USE mp,		ONLY : mp_bcast
    USE mp_world,	ONLY : world_comm

    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this
    ! internal variable
    INTEGER, PARAMETER :: root = 0

    CALL mp_bcast(prefix, root, world_comm)
    CALL mp_bcast(tmp_dir, root, world_comm)
    CALL mp_bcast(this%band_threshold, root, world_comm)
    CALL mp_bcast(this%dt_el, root, world_comm)
    CALL mp_bcast(this%dt_ion, root, world_comm)
    CALL mp_bcast(this%iverbosity, root, world_comm)
    CALL mp_bcast(this%lcorrect_ehrenfest_forces, root, world_comm)
    CALL mp_bcast(this%lcorrect_moving_ions, root, world_comm)
    CALL mp_bcast(this%lprojectile_perturbation, root, world_comm)
    CALL mp_bcast(this%lscalar_perturbation, root, world_comm)
    CALL mp_bcast(this%lvector_perturbation, root, world_comm)
    CALL mp_bcast(this%lxray_perturbation, root, world_comm)
    CALL mp_bcast(this%nsteps_el, root, world_comm)
    CALL mp_bcast(this%nsteps_el_per_nsteps_ion, root, world_comm)
    CALL mp_bcast(this%nsteps_ion, root, world_comm)
    CALL mp_bcast(this%dt_el, root, world_comm)
    CALL mp_bcast(this%dt_ion, root, world_comm)
    IF(this%lprojectile_perturbation) CALL this%projectile_perturbation%broadcast_settings
    IF(this%lscalar_perturbation) CALL this%scalar_perturbation%broadcast_settings
    IF(this%lvector_perturbation) CALL this%vector_perturbation%broadcast_settings
    IF(this%lxray_perturbation) CALL this%xray_perturbation%broadcast_settings
    CALL this%propagator%broadcast_settings

    RETURN

  END SUBROUTINE broadcast_tddft_settings

  SUBROUTINE stop_tddft_calculation(this, lclean_stop)
    !
    ! ... Synchronizes before stopping
    !
    USE mp_global,	ONLY : mp_global_end
    USE parallel_include

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: lclean_stop

    CALL mp_global_end()

    IF(lclean_stop)THEN
      STOP
    ELSE
      STOP 1
    ENDIF

    RETURN

  END SUBROUTINE stop_tddft_calculation
#endif

  SUBROUTINE open_tddft_files(this)
    !
    ! ... Open big files for TDDFT
    !
    USE wvfct,  ONLY : nbnd, npwx
    USE noncollin_module, ONLY : npol
    USE io_files,  ONLY : iunwfc, nwordwfc
    USE buffers, ONLY : open_buffer
    USE control_flags, ONLY : io_level

    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this
    ! internal variable
    LOGICAL :: extant

    ! record length (in real words) of the file containing the KS orbitals
    ! io_level is set to 1 in tddft.f90, but you could change this in principle
    nwordwfc = nbnd*npwx*npol
    CALL open_buffer(iunwfc, 'wfc', nwordwfc, io_level, extant)
    CALL open_buffer(this%iuntdorbs, 'tddft', nwordwfc, io_level, extant)

  END SUBROUTINE open_tddft_files

  SUBROUTINE close_tddft_files(this)
    !
    ! ... Close big files for TDDFT
    !
    USE io_files,  ONLY : iunwfc
    USE buffers,  ONLY : close_buffer

    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this

    CALL close_buffer(iunwfc, 'keep')
    CALL close_buffer(this%iuntdorbs, 'keep')

  END SUBROUTINE close_tddft_files

  SUBROUTINE set_tddft_hamiltonian(this, io_unit, electron_step_counter, ion_step_counter)
    !
    ! ... Sets the Hamiltonian to its value at the specified time step 
    !     This is a part of tddft_mod (for now) because it will rely on perturbations (and eventually, maybe the propagator)
    !
    USE kinds,         ONLY : dp
    USE becmod,        ONLY : becp, is_allocated_bec_type, deallocate_bec_type
    USE cell_base,     ONLY : alat, omega, at, bg
    USE control_flags, ONLY : gamma_only
    USE dfunct,        ONLY : newd
    USE fft_base,      ONLY : dfftp
    USE gvecs,         ONLY : doublegrid
    USE gvect,         ONLY : g, gg, gcutm, gstart, ngm
    USE ions_base,     ONLY : nsp, zv, nat, tau, ityp
    USE ldaU,          ONLY : lda_plus_U
    USE lsda_mod,      ONLY : nspin
    USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, kedtau, vrs
    USE uspp,          ONLY : okvan
    USE pwcom

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit
    INTEGER, INTENT(IN) :: electron_step_counter, ion_step_counter
    ! internal variables
    REAL(dp) :: charge
    REAL(dp) :: eth  ! Hubbard energy that we will never use, but interfaces...
    REAL(dp) :: etotefield  ! energy due to coupling to electric field that we will use...
    REAL(dp), EXTERNAL :: ewald, get_one_electron_shift

    ! start the clock on setting the Hamiltonian
    CALL start_clock('set_hamiltonian')
   
    ! begin by computing the charge density
    ! first, null out real and reciprocal space arrays
    rho%of_g(:,:) = (0.0_dp, 0.0_dp)
    rho%of_r(:,:) = (0.0_dp)
    ! recall, if ovkan is false then this is strictly norm conserving and we don't need to deal with this
    IF(okvan .AND. is_allocated_bec_type(becp)) CALL deallocate_bec_type(becp)
    ! add up across bands
    CALL sum_band()

    ! check to see if LDA+U is on...it shouldn't be
    IF(lda_plus_U) CALL errore('set_tddft_hamiltonian','why are you running with LDA+U?',1)
   
    ! computes the HXC potential *and* energies...recall that pwcom is where MODULE ener lives
    CALL v_of_rho(rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v)

    ! compute the total local potential...this is where we will have to add perturbations...
    CALL setlocal()
    CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)

    ! now that we have the local potential, let's get the non-local potential (if USPP)
    ! TODO: we will certainly want to change the conditions for updating this...
    IF(okvan .AND. (ion_step_counter+electron_step_counter) == 0) CALL newd()

    ! compute the one-body contribution
    deband = get_one_electron_shift() 

    ! compute the Ewald energy
    ewld = ewald(alat, nat, nsp, ityp, zv, at, bg, tau, omega, g, gg, ngm, &
                 gcutm, gstart, gamma_only, strf)

    ! TODO: someday we will want to add other contributions, e.g., vdW and stuff...

    ! sum up the total energy
    etot = eband + deband + (etxc-etxcc) + ewld + ehart

    IF( (ion_step_counter+electron_step_counter) == 0) THEN
      WRITE(io_unit,'(5X, "Initial call to set_hamiltonian...building things before entering main TDDFT loop")')
      WRITE(io_unit,'(5X, "Note: This should match up with the output of the SCF calculation")')
      WRITE(io_unit,'(5X, 16X, 11X, "Total", 8X, "One-body", 9X, "Hartree", 14X, "XC", 11X, "Ewald")')
      WRITE(io_unit,'(5X,"Initial energy",2X,5F16.8)') etot, eband, ehart, etxc+etxcc, ewld
      WRITE(io_unit,*)
    ENDIF
    ! stop the clock on setting the Hamiltonian
    CALL stop_clock('set_hamiltonian')

  END SUBROUTINE set_tddft_hamiltonian
  
  SUBROUTINE allocate_tddft_preloop(this)
  !
  ! ...Allocates space for storing orbitals in the main TDDFT loop
  !
  USE constants,    ONLY : au_sec 
  USE wvfct,        ONLY : nbnd, npwx
  ! 
  IMPLICIT NONE
  ! input variable
  CLASS(tddft_type), INTENT(INOUT) :: this

  ! space for the KS orbitals across all stages of the propagator
  ALLOCATE(this%psi(npwx, nbnd, this%propagator%nstages))
  ! space for the Hamiltonian acting on the KS orbitals
  ALLOCATE(this%Hpsi(npwx, this%nbands_occupied_max))
  ! space for the overlap matrix acting on the KS orbitals
  ALLOCATE(this%Spsi(npwx, this%nbands_occupied_max)) 

  ! initialize to zero...
  this%psi(:,:,:) = (0.0_dp, 0.0_dp)
  this%Hpsi(:,:) = (0.0_dp, 0.0_dp)
  this%Spsi(:,:) = (0.0_dp, 0.0_dp)

  ! if this is an implicit propagator, we need to allocate the RHS
  ! and everything tied to GMRES workspace
  IF(this%propagator%limplicit)THEN
    ALLOCATE(this%rhs(npwx, this%nbands_occupied_max))
    this%rhs(:,:) = (0.0_dp, 0.0_dp)
    CALL this%propagator%implicit_solver%gmres_begin(npwx)
    ! IMPORTANT: cn_factor is set in Rydberg atomic units (thus the extra extra factor of 2...)
    this%propagator%cn_factor = (0.d0,1d-18)*this%dt_el/(4.0_dp*au_sec)

  ENDIF

  END SUBROUTINE allocate_tddft_preloop

  SUBROUTINE deallocate_tddft_postloop(this)
  ! 
  ! ...Deallocates things that were used in the main TDDFT loop
  ! 
  IMPLICIT NONE
  ! input variable
  CLASS(tddft_type), INTENT(INOUT) :: this

  DEALLOCATE(this%psi, this%Hpsi, this%Spsi)
  ! if we are using an implicit propagator, we need to deallocate the RHS
  ! and everything tied to GMRES
  IF(this%propagator%limplicit)THEN
    DEALLOCATE(this%rhs)
    CALL this%propagator%implicit_solver%gmres_end()
  ENDIF

  END SUBROUTINE deallocate_tddft_postloop

  SUBROUTINE propagate_tddft_all_orbitals(this, io_unit, ion_step_counter, electron_step_counter)
    ! 
    ! ... Moves all of the orbitals in the calculation forward by one time step
    ! 
    USE becmod,		ONLY : becp, &
			       allocate_bec_type, calbec, is_allocated_bec_type
    USE buffers,	ONLY : get_buffer, save_buffer
    USE io_files,	ONLY : iunwfc, nwordwfc
    USE klist,		ONLY : igk_k, ngk, nks, xk
    USE lsda_mod,	ONLY : current_spin, isk
    USE uspp,		ONLY : nkb, vkb
    USE wavefunctions,	ONLY : evc
    USE wvfct,		ONLY : current_k, nbnd, npwx
    IMPLICIT NONE
    ! input variable
    CLASS(tddft_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit
    INTEGER, INTENT(IN) :: ion_step_counter, electron_step_counter
    ! internal
    INTEGER :: ierr
    INTEGER :: band_counter, kpt_counter
    INTEGER :: npw
    COMPLEX(dp) :: cn_factor

    IF(this%propagator%limplicit)THEN
      ! grab the Crank-Nicolson factor...
      cn_factor = this%propagator%cn_factor 

      ! loop over k points
      DO kpt_counter = 1, nks

        ! if you don't set these effectively global variables, bad things happen...
        current_k = kpt_counter
        current_spin = isk(kpt_counter)

        ! get the number of plane waves at the current k point
        npw = ngk(kpt_counter)
        ! calculates the kinetic energy array at the current k point
        CALL g2_kin(kpt_counter) 
        ! calculates the Kleinman-Bylander projectors + structure factor in reciprocal space... 
        CALL init_us_2(npw, igk_k(1, kpt_counter), xk(1, kpt_counter), vkb)

        ! read orbitals from the file and calculate becp
        evc = (0.d0, 0.d0)
        ! if this is the very first step, then get the current orbitals from unit iunwfc, else iuntdorbs
        ! presently, we are saving at each time step
        IF(ion_step_counter*electron_step_counter==1)THEN
          CALL get_buffer(evc, nwordwfc, iunwfc, kpt_counter)
        ELSE
          CALL get_buffer(evc, nwordwfc, this%iuntdorbs, kpt_counter)
        ENDIF
        IF(.NOT. is_allocated_bec_type(becp))THEN
          CALL allocate_bec_type(nkb, nbnd, becp)
        ENDIF
        CALL calbec(npw, vkb, evc, becp) 
        
        ! compute matrix-vector products
        CALL tddft_h_psi(npwx, npw, this%nbands_occupied(kpt_counter), evc, this%Hpsi)
        CALL tddft_s_psi(npwx, npw, this%nbands_occupied(kpt_counter), evc, this%Spsi)

        ! load the right hand side
        this%rhs(:,:) = this%Spsi(:,:)-cn_factor*this%Hpsi(:,:)
        ! call GMRES
        CALL this%propagator%implicit_solver%gmres_solve(crank_nicolson_matvec, this%rhs(:,:), this%psi(:,:,1), npwx, npw, nbnd, cn_factor)
        ! for Crank-Nicolson, after the step, move the new wave functions to the position of the 'old' wave functions
        this%psi(:,:,2) = this%psi(:,:,1)

        ! update the wave functions for computation of energies, etc.
        evc(:,:) = this%psi(:,:,1)

        ! save the orbitals at the end of the step...
        CALL save_buffer(this%psi(:,:,1), nwordwfc, this%iuntdorbs, kpt_counter)

      ENDDO
    ELSE
      CALL errore('propagate_tddft_all_orbitals', 'no non-implicit propagators...yet', ierr)
    ENDIF

    RETURN

  END SUBROUTINE propagate_tddft_all_orbitals

  SUBROUTINE crank_nicolson_matvec(npwx, npw, cn_factor, vec_in, vec_out)
    !
    !...Applies the Crank-Nicolson propagator to a vector
    !
    USE kinds,		ONLY : dp
    ! 
    IMPLICIT NONE
    ! input variables
    INTEGER :: npwx, npw
    COMPLEX(dp) :: cn_factor
    COMPLEX(dp), INTENT(IN) :: vec_in(npwx)
    COMPLEX(dp), INTENT(INOUT) :: vec_out(npwx)
    ! internal variable
    COMPLEX(dp) :: vec_tmp(npwx)

    ! store S|vec_in> in |vec_out>
    CALL s_psi(npwx, npw, 1, vec_in, vec_out)
    ! store H|vec_in> in |vec_tmp>
    CALL h_psi(npwx, npw, 1, vec_in, vec_tmp)
    ! compute (S+idt/2 H)|vec_in> and store in |vec_out>
    vec_out(:) = vec_out(:) + cn_factor*vec_tmp(:)

    RETURN

  END SUBROUTINE crank_nicolson_matvec

END MODULE tddft_mod
