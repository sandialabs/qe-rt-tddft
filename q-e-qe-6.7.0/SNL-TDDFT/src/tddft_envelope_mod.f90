! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

MODULE tddft_envelope_mod
  !
  ! ... This module contains the variables used to define a time-dependent envelope for real-time TDDFT calculations
  !
  USE kinds,		ONLY : dp
  USE constants, 	ONLY : pi, tpi
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE tddft_envelope_type

    INTEGER :: envelope_index		! an integer index that uniquely describes the functional form of the envelope
    ! 1 = linear
    ! 2 = Gaussian
    ! 3 = integral of a Gaussian
    ! 4 = derivative of a Gaussian
    ! 5 = pure cosine
    ! 6 = cosine-squared envelope with carrier
    REAL(dp) :: linear_slope
    REAL(dp) :: amplitude
    REAL(dp) :: carrier_frequency
    REAL(dp) :: delay
    REAL(dp) :: width

  CONTAINS
    PROCEDURE :: print_summary => print_envelope_summary
    PROCEDURE :: evaluate => evaluate_envelope
    PROCEDURE :: linear => linear_envelope
    PROCEDURE :: gaussian => gaussian_envelope
    PROCEDURE :: integral_gaussian => integral_gaussian_envelope
    PROCEDURE :: derivative_gaussian => derivative_gaussian_envelope
    PROCEDURE :: pure_cosine => pure_cosine_envelope
    PROCEDURE :: cosine_squared_times_carrier => cosine_squared_times_carrier_envelope
  END TYPE tddft_envelope_type

CONTAINS

  SUBROUTINE print_envelope_summary(this, io_unit)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: io_unit
    ! internal variables
    INTEGER :: ierr

    SELECT CASE (this%envelope_index)
    CASE(1)
      WRITE(io_unit,'(5x," Linear ")')
      WRITE(io_unit,'(5x," Slope                     =",F12.4," (as)^-1")') this%linear_slope
    CASE(2)
      WRITE(io_unit,'(5x," Gaussian ")')
      WRITE(io_unit,'(5x," Amplitude                 =",F12.4," (as) ")') this%amplitude
      WRITE(io_unit,'(5x," Delay                     =",F12.4," (as) ")') this%delay
      WRITE(io_unit,'(5x," Width                     =",F12.4," (as) ")') this%width
    CASE(3)
      WRITE(io_unit,'(5x," Integral of a Gaussian ")')
      WRITE(io_unit,'(5x," Amplitude                 =",F12.4," (as) ")') this%amplitude
      WRITE(io_unit,'(5x," Delay                     =",F12.4," (as) ")') this%delay
      WRITE(io_unit,'(5x," Width                     =",F12.4," (as) ")') this%width
    CASE(4)
      WRITE(io_unit,'(5x," Derivative of a Gaussian ")')
      WRITE(io_unit,'(5x," Amplitude                 =",F12.4," (as) ")') this%amplitude
      WRITE(io_unit,'(5x," Delay                     =",F12.4," (as) ")') this%delay
      WRITE(io_unit,'(5x," Width                     =",F12.4," (as) ")') this%width
    CASE(5)
      WRITE(io_unit,'(5x," Pure cosine ")')
      WRITE(io_unit,'(5x," Amplitude                 =",F12.4," (unitless) ")') this%amplitude
      WRITE(io_unit,'(5x," Carrier frequency         =",F12.4," (eV) ")') this%carrier_frequency
    CASE(6)
      WRITE(io_unit,'(5x," Cosine squared envelope plus carrier ")')
      WRITE(io_unit,'(5x," Amplitude                 =",F12.4," (unitless) ")') this%amplitude
      WRITE(io_unit,'(5x," Carrier frequency         =",F12.4," (eV) ")') this%carrier_frequency
      WRITE(io_unit,'(5x," Width                     =",F12.4," (as) ")') this%width
    CASE DEFAULT
      CALL errore('print_envelope_summary', 'envelope index without implementation', ierr)
    END SELECT
    WRITE(io_unit,'(5x," ")')

    RETURN

  END SUBROUTINE print_envelope_summary

  FUNCTION evaluate_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value
    ! internal variables
    INTEGER :: ierr

    SELECT CASE (this%envelope_index)
    CASE(1)
      envelope_value = this%linear(time)
    CASE(2)
      envelope_value = this%gaussian(time)
    CASE(3)
      envelope_value = this%integral_gaussian(time)
    CASE(4)
      envelope_value = this%derivative_gaussian(time)
    CASE(5)
      envelope_value = this%pure_cosine(time)
    CASE(6)
      envelope_value = this%cosine_squared_times_carrier(time)
    CASE DEFAULT
      CALL errore('evaluate_envelope', 'envelope index without implementation', ierr)
    END SELECT

  END FUNCTION evaluate_envelope

  FUNCTION linear_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value

    envelope_value = this%linear_slope*time

  END FUNCTION linear_envelope

  FUNCTION gaussian_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value
    ! internal variable
    REAL(dp) :: gaussian_argument

    gaussian_argument = (time-this%delay)/(SQRT(2.0_dp)*this%width)
    envelope_value = this%amplitude*EXP(-gaussian_argument*gaussian_argument)/(SQRT(tpi)*this%width)

  END FUNCTION gaussian_envelope

  FUNCTION integral_gaussian_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value
    ! internal variable
    REAL(dp) :: erf_argument
    REAL(dp), EXTERNAL :: qe_erf

    erf_argument = (time-this%delay)/(SQRT(2.0_dp)*this%width)
    envelope_value = 0.5_dp*(1.0_dp + qe_erf(erf_argument))

  END FUNCTION integral_gaussian_envelope

  FUNCTION derivative_gaussian_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value
    ! internal variable
    REAL(dp) :: gaussian_argument

    gaussian_argument = (time-this%delay)/(SQRT(2.0_dp)*this%width)
    envelope_value = -(time-this%delay)*this%amplitude*EXP(-gaussian_argument*gaussian_argument)/(SQRT(tpi)*this%width**3)

  END FUNCTION derivative_gaussian_envelope

  FUNCTION pure_cosine_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value

    envelope_value = this%amplitude * COS(this%carrier_frequency*time)

  END FUNCTION pure_cosine_envelope

  FUNCTION cosine_squared_times_carrier_envelope(this, time) RESULT(envelope_value)

    IMPLICIT NONE
    ! input variables
    CLASS(tddft_envelope_type), INTENT(INOUT) :: this
    REAL(dp) :: time
    ! output variable
    REAL(dp) :: envelope_value

    IF(time <= this%width)THEN
      envelope_value = this%amplitude * ((COS(pi*time/this%width))**2) * COS(this%carrier_frequency*time)
    ELSE
      envelope_value = 0.0_dp
    ENDIF

  END FUNCTION cosine_squared_times_carrier_envelope

END MODULE tddft_envelope_mod
