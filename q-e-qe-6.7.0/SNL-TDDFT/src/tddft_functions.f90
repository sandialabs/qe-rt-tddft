! Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

FUNCTION get_one_electron_shift() RESULT(one_electron_shift)
  ! 
  ! ...Another routine more from CE-TDDFT with light relabeling
  !
  USE kinds,            ONLY : dp
  USE cell_base,        ONLY : omega
  USE fft_base,         ONLY : dfftp
  USE funct,            ONLY : dft_is_meta
  USE ldaU,             ONLY : lda_plus_U
  USE lsda_mod,         ONLY : nspin
  USE mp,               ONLY : mp_barrier, mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : noncolin
  USE paw_variables,    ONLY : okpaw, ddd_paw
  USE scf,              ONLY : rho, scf_type, v
  
  IMPLICIT NONE
  
  INTEGER :: grid_counter
  REAL(dp) :: one_electron_shift

  one_electron_shift = 0.0_dp

  ! add on the contribution from the local potential
  IF(nspin==2)THEN
    DO grid_counter = 1, dfftp%nnr
       one_electron_shift = one_electron_shift &
                            - (rho%of_r(grid_counter,1) + rho%of_r(grid_counter,2))*v%of_r(grid_counter,1) &
                            - (rho%of_r(grid_counter,1) - rho%of_r(grid_counter,2))*v%of_r(grid_counter,2)
    ENDDO
  ELSE
    one_electron_shift = -SUM(rho%of_r(:,:)*v%of_r(:,:))
  ENDIF
  
  ! add on the contribution for meta-GGAs, if relevant
  IF(dft_is_meta()) one_electron_shift = one_electron_shift - SUM(rho%kin_r(:,:)*v%of_r(:,:))

  ! apply the appropriate scale factor for the real-space sum
  one_electron_shift = one_electron_shift*(omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3))

  ! accumulate sums across the group
  CALL mp_sum(one_electron_shift, intra_bgrp_comm)

  IF(lda_plus_U .OR. noncolin) CALL errore('get_one_electron_shift','why are you running w/LDA+U or noncollinear on?',1)

  ! add the PAW contribution...
  IF(okpaw) one_electron_shift = one_electron_shift - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))

  RETURN

END FUNCTION get_one_electron_shift
