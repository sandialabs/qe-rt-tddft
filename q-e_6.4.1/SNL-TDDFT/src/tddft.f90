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
  INTEGER :: electron_step_counter, ion_step_counter
  REAL(dp) :: electron_time, ion_time

  ! initialization
#ifdef __MPI
  CALL mp_startup(start_images=.true.)
#else
  CALL mp_startup(start_images=.false.)
#endif
  CALL environment_start(code)

#ifndef __BANDS
  IF(nbgrp > 1) &
    CALL errore('tddft', 'configure and recompile TDDFT with --enable-band-parallel', 1)
#endif

  WRITE(stdout,*)
  WRITE(stdout,'(5X,''***** SNL-TDDFT git revision '',A,'' *****'')') tddft_git_revision
  WRITE(stdout,*)

  ! create an instance of the tddft_type class that will contain our calculation
  ALLOCATE(this_calculation)

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
  !

  ! check to see whether gamma_only is erroneously being flagged...
  IF(gamma_only) CALL errore('tddft',, 'cannot run TDDFT with gamma_only == .TRUE.',1)

  ! check to see whether wf_collect was used in the SCF
#ifdef __BANDS
  IF( nbgrp>1 .AND. (twfcollect==.FALSE.) )&
    CALL errore('tddft', 'cannot use band-parallelization without wf_collect in SCF', 1)     
#endif  
    
  ! check to see whether noncolin is erroneously being flagged...  
  IF( noncolin ) CALL errore('tddft', 'non-collinear spin is not yet supported', 1)

  ! pluginization stuff from Davide's code
  nat_ = nat
  ntyp_ = ntyp
  ibrav_ = ibrav
  assume_isolated_ = 'none'  

  ! allocate TDDFT stuff
  !

  ! TDDFT setup
  !

  ! TDDFT summary
  !

  ! initialize parallelization over bands
#ifdef __BANDS
  CALL init_parallel_over_band(inter_bgrp_comm, nbnd)
#endif  

  ! main TDDFT loop
  ! 
  electron_time = 0.0_dp
  ion_time = 0.0_dp
  DO ion_step_counter = 1, this_calculation%nsteps_ion
 
      DO electron_step_counter = 1, this_calculation%nsteps_el_per_nsteps_ion

          electron_time = electron_time + this_calculation%dt_el
	  WRITE(stdout,*) 'electron time: ', electron_time,  this_calculation%scalar_perturbation%scalar_envelope%evaluate(electron_time)

      ENDDO

      ion_time = ion_time + this_calculation%dt_ion

  ENDDO

  ! close the files that were opened at the beginning of the TDDFT calculation
  !

  ! print timings
  ! 

  CALL environment_end(code)

  ! 
  CALL this_calculation%stop_calculation( .TRUE. )
  
STOP

END PROGRAM tddft
