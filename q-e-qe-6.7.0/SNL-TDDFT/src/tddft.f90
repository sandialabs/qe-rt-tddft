! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

PROGRAM tddft

  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode, meta_ionode, meta_ionode_id
  USE mp,              ONLY : mp_bcast
  !  USE check_stop,      ONLY : check_stop_init
  USE becmod,          ONLY : becp, allocate_bec_type, is_allocated_bec_type, deallocate_bec_type
  USE constants,       ONLY : amu_ry, rytoev
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE dynamics_module, ONLY : allocate_dyn_vars, deallocate_dyn_vars, vel, verlet
  USE ions_base,       ONLY : amass, if_pos, tau
  USE mp_global,       ONLY : mp_startup
  USE mp_bands,        ONLY : nbgrp
  USE mp_world,        ONLY : world_comm
  USE environment,     ONLY : environment_start, environment_end
  USE wvfct,           ONLY : nbnd
  USE io_global,       ONLY : stdout
  USE noncollin_module,ONLY : noncolin
  ! for pluginization
  USE input_parameters, ONLY : nat_ => nat, ntyp_ => ntyp
  USE input_parameters, ONLY : assume_isolated_ => assume_isolated, &
                               ibrav_ => ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE cell_base,        ONLY : ibrav
  USE pwcom
  USE tddft_mod
  USE tddft_version
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   		:: code = 'TDDFT'
  LOGICAL, EXTERNAL   		:: check_para_diag
  CLASS(tddft_type), POINTER    :: this_calculation
  !------------------------------------------------------------------------
  INTEGER :: electron_step_counter, ion_step_counter
  REAL(dp) :: electron_time, ion_time

  INTEGER :: scratch_index
  REAL(dp) :: scratch_real

  ! initialization
#ifdef __MPI
  CALL mp_startup(start_images=.true.)
#else
  CALL mp_startup(start_images=.false.)
#endif
  CALL environment_start(code)

  WRITE(stdout,*)
  WRITE(stdout,'(5X,''***** SNL-TDDFT git revision '',A,'' *****'')') tddft_git_revision
  WRITE(stdout,*)

  ! create an instance of the tddft_type class that will contain our calculation
  ALLOCATE(this_calculation)

  ! read calculation settings from an input file
  CALL this_calculation%read_settings_file()

  ! set the IO level
  ! open_buffer will connect specified units to files for direct I/O access
  io_level = 1

  ! read the initial Kohn-Sham orbitals from a file, tmp_dir/prefix+postfix
  CALL read_file()

  ! set use_para_diag in parallel runs
  ! this probably isn't strictly necessary because TDDFT doesn't involve diagonalization...
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  ! open the files used in a TDDFT calculation
  CALL this_calculation%open_files()

  ! check to see whether gamma_only is erroneously being flagged...
  IF(gamma_only) CALL errore('tddft',, 'cannot run TDDFT with gamma_only == .TRUE.',1)

  ! check to see whether wf_collect was used in the SCF
#ifdef __BANDS
  IF( nbgrp>1 .AND. (twfcollect==.FALSE.) )&
  CALL errore('tddft', 'cannot use band-parallelization without wf_collect in SCF', 1)
#endif

  ! check to see whether noncolin is erroneously being flagged...
  IF( noncolin ) CALL errore('tddft', 'non-collinear spin is not yet supported', 1)

  ! do things to keep the code from crashing
  CALL this_calculation%perfunctory_business()

  ! set everything up before heading into the main TDDFT loop
  CALL this_calculation%initialize_calculation(stdout)
 
  ! print summary
  CALL this_calculation%print_summary(stdout)

  ! initialize parallelization over bands
#ifdef __BANDS
  CALL start_clock('init_par_band')
  CALL init_parallel_over_band(inter_bgrp_comm, nbnd)
  CALL stop_clock('init_par_band')
#endif

  ! allocate before the loop
  WRITE(stdout, '(5X, "Pre-loop allocation...")')
  CALL this_calculation%allocate_preloop()
  WRITE(stdout, *)

  ! allocate / deal with variables for ionic motion
  WRITE(stdout, '(5X, "Initializing quantities for molecular dynamics...")')
  CALL allocate_dyn_vars()
  ! note: if_pos is a multiplier for the ith component of the jth atom
  !       it can be set to zero to constrain certain atoms
  !       we should allow for the use of if_pos, at some point, but for now it
  !       is always 1 (i.e., no atoms are constrained) unless a projectile
  !       perturbation is on...
  ALLOCATE(if_pos(3, nat))
  IF(this_calculation%lprojectile_perturbation)THEN
    ! only move the projectile atom
    ! first set a scratch index variable for convenience (projectile index)
    scratch_index = this_calculation%projectile_perturbation%projectile_index
    ! and then a scratch real...(projectile velocity in Rydberg units...gross)
    scratch_real = SQRT(2.d0*(this_calculation%projectile_perturbation%projectile_kinetic_energy)/(amass(scratch_index)*amu_ry*rytoev))
    ! report velocity in Hartree atomic units (because it is what humans use for this
    ! sort of thing) and work with Rydberg atomic units (because it is apparently what computers use...)
    WRITE(stdout, '(5X, "Projectile velocity ",F12.4," (a.u.) ")') scratch_real/2.d0
    if_pos(:, :) = 0
    if_pos(:, scratch_index) = 1
    vel(:, :) = 0.d0
    vel(1, scratch_index) = scratch_real*this_calculation%projectile_perturbation%projectile_velocity(1)
    vel(2, scratch_index) = scratch_real*this_calculation%projectile_perturbation%projectile_velocity(2)
    vel(3, scratch_index) = scratch_real*this_calculation%projectile_perturbation%projectile_velocity(3)
  ELSE
    ! let all of the atoms move
    if_pos(:, :) = 1
    vel(:, :) = 0.d0
  ENDIF
  WRITE(stdout, *)
  
  ! main TDDFT loop
  electron_time = 0.0_dp
  ion_time = 0.0_dp
  WRITE(stdout, '(5X, "Entering main loop")')
  WRITE(stdout, *)
  DO ion_step_counter = 1, this_calculation%nsteps_ion

    ! for 'simple' stopping calculations... just update the position at the beginning of each ionic time step
    IF(this_calculation%lprojectile_perturbation)THEN
      tau(:, scratch_index) = tau(:, scratch_index) + vel(:, scratch_index)*this_calculation%dt_ion
    ENDIF

    DO electron_step_counter = 1, this_calculation%nsteps_el_per_nsteps_ion

      ! move all of the orbitals forward by one time step
      CALL this_calculation%propagate_all_orbitals(stdout, ion_step_counter, electron_step_counter)

      ! update the electron time after the step has been taken
      electron_time = electron_time + this_calculation%dt_el
 
      ! update the Hamiltonian for the next time step
      CALL this_calculation%set_hamiltonian(stdout, electron_step_counter, ion_step_counter)
      
      ! now that the step is complete, we report quantities of interest...
      IF(ionode)THEN
        WRITE(stdout,'(5X, " Time, energy:",2X, 6F16.8)') electron_time, etot, eband, ehart, etxc+etxcc, ewld 
        WRITE(stdout,*)
      ENDIF

    ENDDO

    ! compute forces 
    IF(is_allocated_bec_type(becp)) CALL deallocate_bec_type(becp)
    CALL forces()

    IF(this_calculation%lprojectile_perturbation)THEN
      WRITE(stdout,'(5X, " Projectile position:",2X, 3F16.8)') tau(:, scratch_index)
    ENDIF

    ! update the ion time after the step has been taken
    ion_time = ion_time + this_calculation%dt_ion

  ENDDO
  WRITE(stdout, '(5X, "Exiting main loop")')
  WRITE(stdout, *)

  ! deallocate after the loop
  WRITE(stdout, '(5X, "Post-loop deallocation...")')
  CALL this_calculation%deallocate_postloop()
  ! for Ehrenfest...
  CALL deallocate_dyn_vars()
  WRITE(stdout, *)
 
  ! for Ehrenfest 
  DEALLOCATE(if_pos)

  ! close the files that were opened at the beginning of the TDDFT calculation
  CALL this_calculation%close_files()

  ! print timings
  CALL this_calculation%print_summary_clock(stdout)

  CALL environment_end(code)

  ! synchronize before stopping...
  CALL this_calculation%stop_calculation( .TRUE. )

  STOP

END PROGRAM tddft
