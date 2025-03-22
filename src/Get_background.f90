subroutine get_background(DATA1)

  USE module_global
  IMPLICIT NONE
  include 'netcdf.inc'

  CHARACTER(len=100)                                 ::      background_ensrf,filename1
  REAL                                               ::      DATA1(bnx,bny,bnvar,ensemble)
  CHARACTER(len=10)                                  ::      string_tmp,res
  INTEGER                                            ::      numofen,k
  character(len=8) sday


  write(sday,'(i8)') icdate
  filename1=trim(mozartdir)//'/data/emis/emissions.intermediate.'//sday//'.nc'
!  print*, trim(filename1)

    !****************************************************
    !           open ncfiles
    !****************************************************

    retval = nf_open(trim(filename1), NF_NOWRITE, ncid)
    if (retval .ne. nf_noerr) call handle_err_bkg(retval)


    retval = nf_inq_varid    ( ncid, "bgens", varid )
    if (retval .ne. nf_noerr) call handle_err_bkg(retval)
    retval = nf_get_var_real ( ncid, varid, DATA1 )
    if (retval .ne. nf_noerr) call handle_err_bkg(retval)

    !close ncfile
    retval = nf_close(ncid)



END subroutine

!*********************************************8


subroutine handle_err_bkg(errcode)
  implicit none
  include 'netcdf.inc'
  integer errcode

  print *, 'Get Ensrf Background Error: ', nf_strerror(errcode)
  stop 2
end
!***************************
