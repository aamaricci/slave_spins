MODULE SLAVE_SPINS
  USE SS_INPUT_VARS, only: &
       ss_read_input,                                      &
       Nlat,Norb,Nspin,filling,Uloc,Ust,Jh,xmu,beta,eps,   &
       wini,wfin,Nloop,Nsuccess,Lmats,Lreal,Pfile,LOGfile

  USE SS_IO,    only: &
       ss_get_dens,&
       ss_get_zeta,&
       ss_get_lambda,&
       ss_get_lambda0,&
       ss_get_self,&
       ss_get_sz

  USE SS_MAIN, only: ss_solve
  implicit none
END MODULE SLAVE_SPINS
