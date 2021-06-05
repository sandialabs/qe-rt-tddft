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
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
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
  CALL init_parallel_over_band(inter_bgrp_comm, nbnd)
#endif

  ! allocate before the loop
  WRITE(stdout, '(5X, "Pre-loop allocation...")')
  CALL this_calculation%allocate_preloop()
  WRITE(stdout, *)

  ! main TDDFT loop
  electron_time = 0.0_dp
  ion_time = 0.0_dp
  WRITE(stdout, '(5X, "Entering main loop")')
  WRITE(stdout, *)
  DO ion_step_counter = 1, this_calculation%nsteps_ion

    DO electron_step_counter = 1, this_calculation%nsteps_el_per_nsteps_ion

      ! move all of the orbitals forward by one time step
      CALL this_calculation%propagate_all_orbitals(stdout, ion_step_counter, electron_step_counter)

      ! update the electron time after the step has been taken
      electron_time = electron_time + this_calculation%dt_el
 
      ! update the Hamiltonian for the next time step
      CALL this_calculation%set_hamiltonian(stdout, electron_step_counter, ion_step_counter)
      
      ! now that the step is complete, we report quantities of interest...
      IF(ionode)THEN
        WRITE(stdout,'(5X, " Energy:",2X, 6F16.8)') electron_time, etot, eband, ehart, etxc+etxcc, ewld 
        WRITE(stdout,*)
      ENDIF

    ENDDO

    ! update the ion time after the step has been taken
    ion_time = ion_time + this_calculation%dt_ion

  ENDDO
  WRITE(stdout, '(5X, "Exiting main loop")')
  WRITE(stdout, *)

  ! deallocate after the loop
  WRITE(stdout, '(5X, "Post-loop deallocation...")')
  CALL this_calculation%deallocate_postloop()
  WRITE(stdout, *)
  
  ! close the files that were opened at the beginning of the TDDFT calculation
  CALL this_calculation%close_files()

  ! print timings
  CALL this_calculation%print_summary_clock(stdout)

  CALL environment_end(code)

  ! synchronize before stopping...
  CALL this_calculation%stop_calculation( .TRUE. )

  STOP

END PROGRAM tddft
