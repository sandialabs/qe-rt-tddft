PROGRAM tddft

  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, meta_ionode, meta_ionode_id
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
  USE tddft_mod
  USE tddft_version
  USE iotk_module  
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   		:: code = 'TDDFT'
  LOGICAL, EXTERNAL   		:: check_para_diag
  CLASS(tddft_type), POINTER    :: this_calculation
  !------------------------------------------------------------------------

  ! initialization
#ifdef __MPI
  CALL mp_startup(start_images=.true.)
#else
  CALL mp_startup(start_images=.false.)
#endif
  CALL environment_start(code)

#ifndef __BANDS
  IF(nbgrp > 1) &
    call errore('tddft', 'configure and recompile TDDFT with --enable-band-parallel', 1)
#endif

  write(stdout,*)
  write(stdout,'(5X,''***** SNL-TDDFT git revision '',A,'' *****'')') tddft_git_revision
  write(stdout,*)

  ! create an instance of the tddft_type class that will contain our calculation
  allocate(this_calculation)

  ! read calculation settings from an input file
  CALL this_calculation%read_settings_file()

  ! set the IO level
  io_level = 1

  ! read the initial Kohn-Sham orbitals from a file, tmp_dir/prefix+postfix
  CALL read_file()
  
  ! set use_para_diag in parallel runs
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  ! open the files used in a TDDFT calculation

  ! check to see whether gamma_only is erroneously being flagged...
  IF(gamma_only) CALL errore('tddft',, 'Cannot run TDDFT with gamma_only == .TRUE.',1)

  ! close the files that were opened at the beginning of the TDDFT calculation





!  call tddft_openfil
!
!  if (gamma_only) call errore ('tdddft_main', 'Cannot run TDFFT with gamma_only == .true. ', 1)
!#ifdef __BANDS
!  if (nbgrp > 1 .and. (twfcollect .eqv. .false.)) &
!    call errore('tddft_main', 'Cannot use band-parallelization without wf_collect in SCF', 1)
!#endif
!  if (noncolin) call errore('tdddft_main', 'non-collinear not supported yet', 1)
!
!  nat_ = nat
!  ntyp_ = ntyp
!  ibrav_ = ibrav
!  assume_isolated_ = 'none'
!  call tddft_allocate()
!  call tddft_setup()
!  call tddft_summary()
!
!#ifdef __BANDS
!  call init_parallel_over_band(inter_bgrp_comm, nbnd)
!#endif
!
!  ! calculation
!  select case (trim(job))
!  case ('optical')
!     if (molecule) then
!        call molecule_optical_absorption
!     else
!        call errore('tddft_main', 'solids are not yet implemented', 1)
!     endif
!
!  case default
!     call errore('tddft_main', 'wrong or undefined job in input', 1)
!
!  end select
!  
!  ! print timings and stop the code
!  call tddft_closefil
!  call print_clock_tddft
!  call environment_end(code)
!  call stop_code( .true. )
  
STOP

END PROGRAM tddft
