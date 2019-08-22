PROGRAM tddft

  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, meta_ionode, meta_ionode_id
  USE mp,              ONLY : mp_bcast
!  USE check_stop,      ONLY : check_stop_init
!  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
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
  CHARACTER (LEN=9)   :: code = 'TDDFT'
  LOGICAL, EXTERNAL   :: check_para_diag
  TYPE(tddft_type)    :: this_calculation
  !------------------------------------------------------------------------

  ! initialize
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
 
  this_calculation%read_settings_file()

!  call tddft_readin()
!  call check_stop_init( max_seconds )
!
!  io_level = 1
! 
!  ! read ground state wavefunctions
!  call read_file
!#ifdef __MPI
!  use_para_diag = check_para_diag(nbnd)
!#else
!  use_para_diag = .false.
!#endif
!
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
