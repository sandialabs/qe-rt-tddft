MODULE tddft_mod
  !
  ! ... This module contains the variables used for TDDFT calculations
  !
  USE kinds,                   ONLY : dp
  USE tddft_perturbations_mod 
  IMPLICIT NONE
  ! 
  SAVE
  !
  TYPE, PUBLIC :: tddft_type

      INTEGER :: &
          iverbosity,		    &! integer indicating the level of verbosity
          nsteps_el,          	    &! total number of electronic steps
          nsteps_el_per_nsteps_ion, &! number of electronic steps per ionic steps
          nsteps_ion		     ! total number of ionic steps
      LOGICAL :: &
          lcorrect_ehrenfest_forces,&! flag = .TRUE. => compute the correct Ehrenfest forces (USPP/PAW feature)
          lcorrect_moving_ions,     &! flag = .TRUE. => compute the "gauge" correction for moving ions (USPP/PAW feature)
          lscalar_perturbation,     &! flag = .TRUE. => applies a perturbation through a homogeneous scalar potential
          lstopping_perturbation,   &! flag = .TRUE. => applies a perturbation by moving an ion at a fixed velocity
          lvector_perturbation,     &! flag = .TRUE. => applies a perturbation through a homogeneous vector potential
          lxray_perturbation         ! flag = .TRUE. => applies a perturbation through an inhomogeneous scalar potential
      REAL(dp) :: &
          dt_el,		    &! electronic time step
          dt_ion,                   &! ionic time step
          duration                   ! total duration in attoseconds
      CLASS(scalar_perturbation_type), POINTER :: scalar_perturbation 
      CLASS(stopping_perturbation_type), POINTER :: stopping_perturbation
      CLASS(vector_perturbation_type), POINTER :: vector_perturbation
      CLASS(xray_perturbation_type), POINTER :: xray_perturbation

      CONTAINS
        PROCEDURE :: read_settings_file => read_tddft_settings   ! reads the settings file for a TDDFT calculation in from a file
#ifdef __MPI	
	PROCEDURE :: broadcast_inputs => broadcast_tddft_inputs  ! broadcasts inputs to all tasks after reading in settings on the IO node
#endif

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
            lscalar_perturbation,     &
            lstopping_perturbation,   &
            lvector_perturbation,     &
            lxray_perturbation       
        REAL(dp) :: dt_el, dt_ion, duration
    
        NAMELIST /tddft/ prefix, tmp_dir, verbosity, &
                         nsteps_el, nsteps_ion, &
                         dt_el, dt_ion, duration, &
          	       lcorrect_ehrenfest_forces, lcorrect_moving_ions, &
          	       lscalar_perturbation, lstopping_perturbation, &
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
            lscalar_perturbation = .FALSE.
            lstopping_perturbation = .FALSE.
            lvector_perturbation = .FALSE.
            lxray_perturbation = .FALSE.
    
            ! default propagation is one step that is 1 as long
            dt_el = 1.0_dp
            dt_ion = 1.0_dp
            duration = 1.0_dp
            nsteps_el = 1
            nsteps_ion = 1
            nsteps_el_per_nsteps_ion = 1
     
            ! read tddft namelist from the input file    
            READ( 5, tddft, err = 200, iostat = ierr )
200 CALL errore('read_settings_file', 'reading tddft namelist', ierr)
    
            this%lcorrect_ehrenfest_forces = lcorrect_ehrenfest_forces
            this%lcorrect_moving_ions = lcorrect_moving_ions
            this%lscalar_perturbation = lscalar_perturbation
            this%lstopping_perturbation = lstopping_perturbation
            this%lvector_perturbation = lvector_perturbation
            this%lxray_perturbation = lxray_perturbation
    
            ! TODO: sanitize this, someday
            this%nsteps_el = CEILING(duration/dt_el)
            this%dt_el = duration/nsteps_el
            this%nsteps_ion = CEILING(duration/dt_ion)
            this%dt_ion = duration/nsteps_ion
            this%nsteps_el_per_nsteps_ion = nsteps_el/nsteps_ion
    
            ! set the integer verbosity flag
            SELECT CASE (verbosity)
                CASE('low')
                    this%iverbosity = 1
                CASE('medium')
                    this%iverbosity = 11
                CASE('high')
                    this%iverbosity = 21
                CASE DEFAULT
                    CALL errore('read_settings_file', 'verbosity can be ''low'', ''medium'' or ''high''', 1)
            END SELECT
    
            IF(this%lscalar_perturbation)THEN
                ALLOCATE(this%scalar_perturbation)
                CALL this%scalar_perturbation%read_settings_file()
            ENDIF
    
            IF(this%lstopping_perturbation)THEN
                ALLOCATE(this%stopping_perturbation)
                CALL this%stopping_perturbation%read_settings_file()
            ENDIF
    
            IF(this%lvector_perturbation)THEN
                ALLOCATE(this%vector_perturbation)
                CALL this%vector_perturbation%read_settings_file()
            ENDIF
    
            IF(this%lxray_perturbation)THEN
                ALLOCATE(this%xray_perturbation)
                CALL this%xray_perturbation%read_settings_file()
            ENDIF
     
        ENDIF
    
#ifdef __MPI
        ! broadcast input variables  
        CALL this%broadcast_inputs
#endif

        RETURN     
    
    END SUBROUTINE read_tddft_settings

#ifdef __MPI
    SUBROUTINE broadcast_tddft_inputs(this)
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
	CALL mp_bcast(this%dt_el, root, world_comm)
	CALL mp_bcast(this%dt_ion, root, world_comm)
	CALL mp_bcast(this%iverbosity, root, world_comm)
	CALL mp_bcast(this%lcorrect_ehrenfest_forces, root, world_comm)
	CALL mp_bcast(this%lcorrect_moving_ions, root, world_comm)
	CALL mp_bcast(this%lscalar_perturbation, root, world_comm)
	CALL mp_bcast(this%lstopping_perturbation, root, world_comm)
	CALL mp_bcast(this%lvector_perturbation, root, world_comm)
	CALL mp_bcast(this%lxray_perturbation, root, world_comm)
	CALL mp_bcast(this%nsteps_el, root, world_comm)
	CALL mp_bcast(this%nsteps_el_per_nsteps_ion, root, world_comm)
	CALL mp_bcast(this%nsteps_ion, root, world_comm)
	CALL mp_bcast(this%dt_el, root, world_comm)
	CALL mp_bcast(this%dt_ion, root, world_comm)
        IF(ASSOCIATED(this%scalar_perturbation)) CALL this%scalar_perturbation%broadcast_inputs
        IF(ASSOCIATED(this%stopping_perturbation)) CALL this%stopping_perturbation%broadcast_inputs
        IF(ASSOCIATED(this%vector_perturbation)) CALL this%vector_perturbation%broadcast_inputs
        IF(ASSOCIATED(this%xray_perturbation)) CALL this%xray_perturbation%broadcast_inputs

        RETURN

    END SUBROUTINE broadcast_tddft_inputs
#endif

END MODULE tddft_mod
