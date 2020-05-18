! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

MODULE it_solver_mod
  !
  ! ... Class that provides the iterative solution to Ax = b
  !
  USE kinds,	ONLY : dp
  !
  SAVE
  !
  TYPE it_solver_type

    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: &
    work		     ! working space
    INTEGER, ALLOCATABLE, DIMENSION(:) :: &
    inner_counts,	    &! counter for the inner loop on each RHS
    outer_counts             ! counter for the outer loop on each RHS
    INTEGER :: &
    number_rhs,             &! number of right hand sides (RHS), almost certainly the number of orbitals per k-point
    max_kry_dim = 15,	    &! maximum dimension of the Krylov subspace before a restart
    max_restarts = 200	     ! maximum number of restarts before calling it quits
    REAL(dp) :: &
    tol = 1.E-12_dp          ! tolerance for the iterative solve

  CONTAINS
    PROCEDURE :: set_max_kry_dim	! sets the maximum dimension of the Krylov subspace before a restart
    PROCEDURE :: set_max_restarts   ! sets the maximum number of restarts before calling it quits
    PROCEDURE :: set_tol		! sets the tolerance for *any* iterative solver
    PROCEDURE :: gmres_begin	! allocates the working space for GMRES
    PROCEDURE :: gmres_end		! deallocates the working space for GMRES
    PROCEDURE :: gmres_solve	! solves A*x = b

  END TYPE it_solver_type

CONTAINS

  SUBROUTINE set_max_kry_dim(this, max_kry_dim)
    !
    ! ... interface for setting the maximum dimension of the Krylov subspace for the iterative solver
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: max_kry_dim

    this%max_kry_dim = max_kry_dim

    RETURN

  END SUBROUTINE set_max_kry_dim

  SUBROUTINE set_max_restarts(this, max_restarts)
    !
    ! ... interface for setting the maximum number of restarts before calling it quits
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: max_restarts

    this%max_restarts = max_restarts

    RETURN

  END SUBROUTINE set_max_restarts

  SUBROUTINE set_tol(this, tol)
    !
    ! ... interface for setting the tolerance of the iterative solver
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this
    REAL(dp), INTENT(IN) :: tol

    this%tol = tol

    RETURN

  END SUBROUTINE set_tol

  SUBROUTINE gmres_begin(this, ndimx)
    !
    ! ... allocates the working space for GMRES
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: ndimx
    ! internal variable
    INTEGER :: ierr

    ! allocate the working space
    ALLOCATE(this%work(ndimx,3+this%max_kry_dim), stat=ierr)
    ! for GMRES, the working space is laid out as follows
    !    1st column = residual
    !    2nd column = w
    !    remaining columns = Krylov iterates

    ! check that the allocation was OK
    IF(ierr/=0) CALL errore('gmres_begin','error allocating',ierr)

    RETURN

  END SUBROUTINE gmres_begin

  SUBROUTINE gmres_end(this)
    !
    ! ... deallocates the working space
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this

    DEALLOCATE(this%work)

    RETURN

  END SUBROUTINE gmres_end

  SUBROUTINE gmres_solve(this, A, b, x, ndimx, ndim, nbnd)
    !
    ! ... solves A*x=b
    !
    USE mp_global,	ONLY : intra_pool_comm
    USE mp,		ONLY : mp_sum
    !
    IMPLICIT NONE
    ! input variables
    CLASS(it_solver_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: ndimx		! the maximum dimension of the vectors
    INTEGER, INTENT(IN) :: ndim			! the actual dimension of the vectors
    INTEGER, INTENT(IN) :: nbnd			! the number of bands
    EXTERNAL :: A 				! the subroutine that computes A*x
    COMPLEX(dp), INTENT(IN) :: b(ndimx,nbnd)	! the matrix of right hand sides (each column is an orbital/band)
    COMPLEX(dp), INTENT(INOUT) :: x(ndimx,nbnd) ! the matrix of solution vectors (each column is an orbital/band)
    ! internal variables
    COMPLEX(dp), EXTERNAL :: zdotc
    COMPLEX(dp) :: tmp
    INTEGER :: ierr				! error flag
    INTEGER :: irhs, inner, outer, k		! counters for loops
    REAL(dp) :: error, normb, normr, tolb	! variables germane to convergence
    LOGICAL :: converged_flag			! true = converged, false = not converged

    ! quantities projected into the subspace
    COMPLEX(dp) :: cs(this%max_kry_dim), sn(this%max_kry_dim)
    COMPLEX(dp) :: s(this%max_kry_dim+1), y(this%max_kry_dim+1)
    COMPLEX(dp) :: h(this%max_kry_dim+1,this%max_kry_dim)


    converged_flag = .FALSE.	! of course, we haven't converged yet
    this%number_rhs = nbnd	! set the number of RHS to the number of bands
    ALLOCATE(this%inner_counts(this%number_rhs), stat=ierr)
    IF(ierr/=0) CALL errore('gmres_solve','error allocating inner_counts',ierr)
    ALLOCATE(this%outer_counts(this%number_rhs), stat=ierr)
    IF(ierr/=0) CALL errore('gmres_solve','error allocating outer_counts',ierr)

    IF(.NOT.ALLOCATED(this%work)) CALL errore('gmres_solve','work not allocated', 1)

    DO irhs = 1, this%number_rhs ! loop over orbitals/bands

      this%inner_counts(irhs) = 0
      this%outer_counts(irhs) = 0

      ! compute the norm of the RHS for the current orbital/band
      normb = DBLE(zdotc(ndim, b(1,irhs), 1, b(1,irhs), 1))
#ifdef __MPI
      CALL mp_sum(normb, intra_pool_comm)
#endif
      normb = DSQRT(normb)
      tolb = normb*this%tol

      ! interface with Michele / compute A*x0, store in the second column of this%work
      !CALL A( ) ! x(1,irhs), this%work(:,2)
      ! copy A*x0 into the residual/first column of the working space
      CALL zcopy(ndim, this%work(:,2), 1, this%work(:,1), 1)
      ! subtract off b (i.e., the residual is A*x0-b)
      CALL zaxpy(ndim, -1.0_dp, b(1,irhs), 1, this%work(:,1), 1)
      ! change the sign
      CALL zscal(ndim, -1.0_dp, this%work(:,1), 1)

      ! compute the norm of the residual
      normr = DBLE(zdotc(ndim, this%work(:,1), 1, this%work(:,1), 1))
#ifdef __MPI
      CALL mp_sum(normr, intra_pool_comm)
#endif
      normr = DSQRT(normr)

      IF(normr<tolb)THEN
        converged_flag = .TRUE.  ! what to do with this...
        CYCLE
      ENDIF

      ! compute the Givens rotation coefficients
      cs(:) = CMPLX(0.0_dp,0.0_dp)
      sn(:) = CMPLX(0.0_dp,0.0_dp)
      ! Hessenberg matrix
      h(:,:) = CMPLX(0.0_dp,0.0_dp)
      s(1) = normr
      this%work(:,3) = this%work(:,1)/normr ! normalized b-A*x0

      DO outer = 1, this%max_restarts ! loop over restarts

        this%outer_counts(irhs) = outer

        IF(outer/=1)THEN ! we can save a matvec on the first iteration, because we know b-A*x0 was just computed
          ! compute A*x, store it in the space for w for efficiency
          !CALL A( ) ! x(1,irhs), this%work(:,2)
          ! copy A*x into the residual
          CALL zcopy(ndim, this%work(:,2), 1, this%work(:,1), 1)
          ! subtract off b
          CALL zaxpy(ndim, -1.0_dp, b(1,irhs), 1, this%work(:,1), 1)
          ! change the sign
          CALL zscal(ndim, -1.0_dp, this%work(:,1), 1)
          ! compute the norm of the residual
          normr = DBLE(zdotc(ndim, this%work(:,1), 1, this%work(:,1), 1))
#ifdef __MPI
          CALL mp_sum(normr, intra_pool_comm)
#endif
          normr = DSQRT(normr)
          s(1) = normr
          ! for starting the next inner loop...
          this%work(:,3) = this%work(:,1)/normr
        ENDIF

        DO inner = 1, this%max_kry_dim ! loop over the Krylov subspace

          this%inner_counts(irhs) = inner

          ! compute w
          !CALL A( ) ! this%work(:,2+inner), this%work(:,2)
          DO k = 1, inner

            tmp = zdotc(ndim, this%work(:,2+k), 1, this%work(:,2), 1)
#ifdef __MPI
            CALL mp_sum(tmp, intra_pool_comm)
#endif
            h(k,inner) = tmp
            ! w = w-h*v
            CALL zaxpy(ndim, -h(k,inner), this%work(:,2+k), 1, this%work(:,2), 1)

          ENDDO

          ! compute the last Hessenberg entry
          tmp = zdotc(ndim, this%work(:,2), 1, this%work(:,2), 1)
#ifdef __MPI
          CALL mp_sum(tmp, intra_pool_comm)
#endif
          h(inner+1,inner) = SQRT(tmp)

          this%work(:,3+inner) = this%work(:,2)/h(inner+1,inner)

          DO k = 1, inner-1
            tmp = CONJG(cs(k))*h(k,inner)+CONJG(sn(k))*h(k+1,inner)
            h(k+1,inner) = -sn(k)*h(k,inner) + cs(k)*h(k+1,inner)
            h(k,inner) = tmp
          ENDDO

          CALL givens_rotation(h(inner,inner), h(inner+1,inner), cs(inner), sn(inner))
          tmp = CONJG(cs(inner))*s(inner)
          s(inner+1) = -sn(inner)*s(inner)
          s(inner) = tmp
          h(inner,inner) = CONJG(cs(inner))*h(inner,inner) + CONJG(sn(inner))*h(inner+1,inner)

          error = abs(s(inner+1))/normb
          IF(error<=this%tol)THEN
            y(1:inner) = s(1:inner)
            CALL gaussian_elimination(this%max_kry_dim+1, inner, h, y)
            DO k = 1, inner
              x(:,irhs) = x(:,irhs) + this%work(:,2+k)*y(k)
            ENDDO
            RETURN
          ENDIF

        ENDDO ! end loop over the Krylov subspace

        y(1:this%max_kry_dim) = s(this%max_kry_dim)
        CALL gaussian_elimination(this%max_kry_dim+1, this%max_kry_dim, h, y)
        DO k = 1, this%max_kry_dim
          x(:,irhs) = x(:,irhs) + this%work(:,2+k)*y(k)
        ENDDO

        ! compute the residual, just like the other two times
        !CALL A( ) ! x(1,irhs), this%work(:,2)
        CALL zcopy(ndim, this%work(:,2), 1, this%work(:,1), 1)
        CALL zaxpy(ndim, -1.0_dp, b(1,irhs), 1, this%work(:,1), 1)
        CALL zscal(ndim, -1.0_dp, this%work(:,1), 1)
        normr = DBLE(zdotc(ndim, this%work(:,1), 1, this%work(:,1), 1))
#ifdef __MPI
        CALL mp_sum(normr, intra_pool_comm)
#endif
        normr = DSQRT(normr)
        s(this%max_kry_dim+1) = normr

        error = abs(s(this%max_kry_dim+1))/normb
        IF(error<=this%tol)THEN
          RETURN
        ENDIF

      ENDDO ! end loop over restarts

    ENDDO ! end loop over orbitals/bands

    ! check flag
    IF(.NOT. converged_flag)THEN
      CALL errore('gmres_solve', 'gmres did not achieve convergence, increase maximum inner/outer iterations', 1)
      STOP
    ENDIF

    RETURN

  END SUBROUTINE gmres_solve

END MODULE it_solver_mod

SUBROUTINE gaussian_elimination(m, n, A, b)
  !
  ! ... simple implementation of Gaussian elimination, as required by GMRES
  !
  USE kinds,	ONLY : dp
  IMPLICIT NONE
  ! input variables
  INTEGER, INTENT(IN) :: m, n
  COMPLEX(dp), INTENT(INOUT) :: A(m,m-1), b(m)
  ! internal variables
  INTEGER :: i, j, k
  INTEGER :: ksave(1)
  COMPLEX(dp) :: tmp_array(size(A,1))
  COMPLEX(dp) :: summ, ratio, tmp

  DO i = 1, n-1 ! loop over the rows of A

    ! find the entry in this row that is largest (for stability)
    ksave = maxloc(abs(A(i:n,i)))

    ! shift the index to compensate for the slicing in the call above (i.e., i::n)
    k = ksave(1) + i - 1

    IF(k /= i)THEN
      tmp_array(1:n) = A(i,1:n)
      A(i,1:n) = A(k,1:n)
      A(k,1:n) = tmp_array(1:n)
      tmp = b(i)
      b(i) = b(k)
      b(k) = tmp
    ENDIF

    DO j = i+1, n ! loop over all rows of A below the current one

      ratio = A(j,i)/A(i,i)
      A(j,1:n) = A(j,1:n) - ratio * A(i,1:n)
      b(j) = b(j) - ratio * b(i)

    ENDDO ! end loop over remaining rows of A

  ENDDO ! end loop over rows of A

  DO i = n, 1, -1 ! loop backwards over the elements of b
    summ = b(i)
    DO j = i+1, n
      summ = summ - A(i,j) * b(j)
    ENDDO
    b(i) = summ/A(i,i)
  ENDDO ! end reverse loop over the elements of b

  RETURN

END SUBROUTINE gaussian_elimination

SUBROUTINE givens_rotation(a, b, c, s)
  !
  ! ... simple implementation of a Givens rotation, as required by GMRES
  !
  USE kinds,	ONLY : dp
  IMPLICIT NONE
  ! input variables
  COMPLEX(dp) :: a, b, c, s
  ! internal variable
  COMPLEX(dp) :: tmp

  IF(b == cmplx(0.0_dp,0.0_dp) )THEN
    c = cmplx(1.0_dp,0.0_dp)
    s = cmplx(0.0_dp,0.0_dp)
  ELSEIF(abs(b)>abs(a))THEN
    tmp = abs(a)/b
    s = 1.0_dp/sqrt(1.0_dp+tmp*tmp)
    c = a*s/b
  ELSE
    tmp = b/abs(a)
    c = a/abs(a)/sqrt(1.0_dp+tmp*tmp)
    s = tmp/sqrt(1.0_dp+tmp*tmp)
  ENDIF

END SUBROUTINE givens_rotation
