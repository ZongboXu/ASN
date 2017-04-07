PROGRAM test
IMPLICIT NONE
INTEGER(KIND=4)::NPTS=0,num=0
CHARACTER(len=100),allocatable,DIMENSION(:)::arrsname
CHARACTER(len=100),allocatable,DIMENSION(:)::arrtname
CHARACTER(128)::fname
INTEGER::snum,tnum,fid,stsensor,endsensor
INTEGER::i=0,sn1,sn2
INTEGER::status
!using arrt1/2 to repesent data in time domain ,al
!so frequence domain
OPEN(12,FILE='snamenum',STATUS='OLD')
  READ(12,*) snum
CLOSE(12)
OPEN(12,FILE='tnamenum',STATUS='OLD')
  READ(12,*) tnum
CLOSE(12)
allocate(arrsname(0:snum-1),STAT=status)
allocate(arrtname(0:tnum-1),STAT=status)
!loop
OPEN(12,FILE='loop',STATUS='OLD')
  READ(12,*) stsensor,endsensor
CLOSE(12)
!
fname="snamelist"
OPEN(12,FILE=fname,STATUS='OLD')
  READ(12,*) (arrsname(i),i=0,snum-1)
CLOSE(12)
fname="tnamelist"
OPEN(12,FILE=fname,STATUS='OLD')
  READ(12,*) (arrtname(i),i=0,tnum-1)
CLOSE(12)

DO sn1=stsensor-1,endsensor-1
!$OMP PARALLEL DO shared(arrsname,arrtname,sn1,snum,tnum)
DO sn2=0,snum-1
    CALL stcoherence(arrsname,arrtname,sn1,sn2,snum,tnum)
END DO
!$OMP END PARALLEL DO

END DO
  DEALLOCATE(arrsname,STAT=status)
  DEALLOCATE(arrtname,STAT=status)
END PROGRAM test

SUBROUTINE stcoherence(arrsname,arrtname,sn1,sn2,snum,tnum)
IMPLICIT NONE
INTEGER::i=0,sn1,sn2,tn,n=0
INTEGER::status,fid=12,tnum,snum
CHARACTER(len=100),DIMENSION(0:snum-1)::arrsname
CHARACTER(len=100),DIMENSION(0:tnum-1)::arrtname
CHARACTER(128)::fname
REAL(KIND=8)::mm=0.d0,c=0.001
INTEGER(KIND=8)::p1,p2,p3
REAL(KIND=4)::delta=0.,maxarr=0.
INTEGER(KIND=4)::NPTS=0,num=0
REAL,allocatable,DIMENSION(:)::head,coharr
DOUBLE COMPLEX,allocatable,DIMENSION(:)::arrt1,arrt2,arrf1,arrf2,farr,fa1,fa2
INTEGER,PARAMETER::fftw_measure=1
INTEGER,PARAMETER::fftw_backward=1
INTEGER,PARAMETER::fftw_forward=-1
DO tn=0,tnum-1 
  fname=trim(adjustl(arrsname(sn1)))//'_'//trim(adjustl(arrtname(tn)))
  status=access(fname," ")
  IF(status.eq.0) EXIT
END DO
  CALL sacinfo(fid,fname,NPTS,delta)
  num=2*NPTS-1
  allocate(head(0:NPTS-1),STAT=status)
  allocate(arrt1(0:num-1),STAT=status)
  allocate(arrt2(0:num-1),STAT=status)
  allocate(arrf1(0:num-1),STAT=status)
  allocate(arrf2(0:num-1),STAT=status)
  allocate(coharr(0:num-1),STAT=status)
  allocate(farr(0:num-1),STAT=status)
  allocate(fa1(0:num-1),STAT=status)
  allocate(fa2(0:num-1),STAT=status)
!*****************************************
!prepare to do fft to two arraries
  call dfftw_plan_dft_1d(p1,num,arrt1(0:num-1),arrf1(0:num-1),fftw_forward,fftw_measure)
  call dfftw_plan_dft_1d(p2,num,arrt2(0:num-1),arrf2(0:num-1),fftw_forward,fftw_measure)
  call dfftw_plan_dft_1d(p3,num,arrt2(0:num-1),arrt1(0:num-1),fftw_backward,fftw_measure)
!**********************************
  coharr=0.0
  farr=CMPLX(0.d0,0.d0)
  fa1=CMPLX(0.d0,0.d0)
  fa2=CMPLX(0.d0,0.d0)
  fid=sn1*10+sn2+10
  n=0
DO tn=0,tnum-1
  arrt1=CMPLX(0.d0,0.d0)
  arrt2=CMPLX(0.d0,0.d0)
  arrf1=CMPLX(0.d0,0.d0)
  arrf2=CMPLX(0.d0,0.d0)
!************examination of two files' existence
  fname=trim(adjustl(arrsname(sn1)))//'_'//trim(adjustl(arrtname(tn)))
  status=access(fname," ")
  IF(status.ne.0)CYCLE
  fname=trim(adjustl(arrsname(sn2)))//'_'//trim(adjustl(arrtname(tn)))
  status=access(fname," ")
  IF(status.ne.0)CYCLE 
 n=n+1  
!***********************************************
  fname=trim(adjustl(arrsname(sn1)))//'_'//trim(adjustl(arrtname(tn)))
  CALL rdsac(fid,fname,head,NPTS)
  DO i=0,NPTS-1
     arrt1(i)=CMPLX(head(i),0)
  END DO
!another arrary for file
  fname=trim(adjustl(arrsname(sn2)))//'_'//trim(adjustl(arrtname(tn)))
  CALL rdsac(fid,fname,head,NPTS)
!**********************************
  DO i=0,NPTS-1
    arrt2(i)=CMPLX(head(i),0)
  END DO
!************fft for data*************************
  call dfftw_execute_dft(p1,arrt1(0:num-1),arrf1(0:num-1))
  call dfftw_execute_dft(p2,arrt2(0:num-1),arrf2(0:num-1))
!************************************************* 
!make arrt2 complex conjura
  arrf2(0:num-1)=CONJG(arrf2(0:num-1))
  mm= MAXVAL(abs(arrf1(1:num-1))*abs(arrf2(1:num-1)))
!cross-coherence 
!******************1 
  DO i=0,num-1
   arrt2(i)=(arrf1(i)*arrf2(i))/MAX(ABS(arrf1(i))*ABS(arrf2(i)),c*mm) ! (ABS(arrf1(i))*ABS(arrf2(i))) 
  END DO
!*********************
    farr(0:num-1)=farr(0:num-1)+arrt2(0:num-1)
!****************ifft******************************
!  call dfftw_execute_dft(p3,arrt2,arrt1)
!    maxarr=maxval(REAL(arrt1(0:num-1)))
!  coharr(0:num-1)=coharr+REAL(arrt1(0:num-1))/maxarr
END DO
!****************time domain************************
  call dfftw_execute_dft(p3,farr(0:num-1),arrt1(0:num-1))
    maxarr=maxval(REAL(arrt1(0:num-1)))
  coharr(0:num-1)=REAL(arrt1(0:num-1))/maxarr
  fname=trim(arrsname(sn1))//trim(arrsname(sn2))//'.coh'
  WRITE(*,*) fname
  CALL wsac(fid,fname,coharr(0:num-1),num,delta)
!*************frequency domain*****************************
  coharr(0:num-1)=REAL(farr(0:num-1))/n
  fname=trim(arrsname(sn1))//trim(arrsname(sn2))//'.fcoh'
  WRITE(*,*) fname
  CALL wsac(fid,fname,coharr(0:num-1),num,delta)  
!**************destroy p1/2/3*******************
  call dfftw_destroy_plan(p1)
  call dfftw_destroy_plan(p2)
  call dfftw_destroy_plan(p3)
!***************************************************
!***************************************************
  DEALLOCATE(head,STAT=status)
  DEALLOCATE(arrt1,STAT=status)
  DEALLOCATE(arrt2,STAT=status)
  DEALLOCATE(arrf1,STAT=status)
  DEALLOCATE(arrf2,STAT=status)
  DEALLOCATE(coharr,STAT=status)
  DEALLOCATE(farr,STAT=status)
END SUBROUTINE stcoherence
