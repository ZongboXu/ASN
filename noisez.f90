PROGRAM znoiserayleigh
!This program is aimed at simulating vertical/z component data of ambient seismic noise.
!Assuming there are only fundamental Rayleigh waves in noise. 
!Using plane wave assumption to simulate noise in frequency domain(frecord) then transform to time domain(trecord). 
!Loading surface-wave phase velocities by reading 'rayleighpv.dat'.
!Saving simulated noise in SAC files, named as H1,2...
!Reference Lawrence 2012 JGR.
!Author: Zongbo XU
!$use omp_lib
IMPLICIT NONE
CHARACTER(len=128) :: fname,addnumber
INTEGER::nsource,nreceiver=60,fre=100,fm=10,t=3600,twindow=3600
INTEGER::i,nr,ns,nf,nt,npv,npvf,status
REAL,ALLOCATABLE,DIMENSION(:)::sacdata
DOUBLE COMPLEX,ALLOCATABLE,DIMENSION(:)::arrt,arrf
REAL,ALLOCATABLE,DIMENSION(:,:)::sourcexy,receiverxy,cv
DOUBLE COMPLEX,ALLOCATABLE,DIMENSION(:,:)::frecord,trecord
REAL::ang,pi=3.1415926,c,a=0,f,r,k,dt,distance,srt,dpv,dpf,x,timeb,timee!a means attenuate 
INTEGER(KIND=8)::p1
INTEGER,PARAMETER::fftw_measure=1
INTEGER,PARAMETER::fftw_backward=1,fftw_forward=-1
!calculate number of sources
i=0
CALL CPU_TIME(timeb)
!nr here mean r

!nsource if noise sources are evenly distributed, nsource better > 3600
!if only inline or offline sources, nsource can equal 500.
nsource=1500

ALLOCATE(sourcexy(0:nsource-1,1:2),STAT=status)
ALLOCATE(receiverxy(0:nreceiver-1,1:2),STAT=status)
ALLOCATE(frecord(0:nreceiver-1,0:t*fre-1),STAT=status)
ALLOCATE(trecord(0:nreceiver-1,0:t*fre-1),STAT=status)
ALLOCATE(arrt(0:t*fre-1),STAT=status)
ALLOCATE(arrf(0:t*fre-1),STAT=status)
ALLOCATE(sacdata(0:t*fre-1),STAT=status)
call dfftw_plan_dft_1d(p1,t*fre,arrt(0:t*fre-1),arrf(0:t*fre-1),fftw_forward,fftw_measure)
!wavelet
dt=1/real(fre)
DO i=0,t*fre-1
  srt=i*dt-1
  r=(1-2*(pi*fm*srt)*(pi*fm*srt))*EXP(-1*(pi*fm*srt)*(pi*fm*srt))
  arrt(i)=CMPLX(r,0)
END DO

!read phase velocity
fname='rayleighpv.dat'
OPEN(12,FILE=fname,STATUS='old')
READ(12,*) npv
ALLOCATE(cv(0:1,0:npv-1),STAT=status)
DO i=0,npv-1
  READ(12,*) cv(0,i),cv(1,i)
END DO
CLOSE(12)

!wavelet in frequency domain
CALL dfftw_execute_dft(p1,arrt(0:t*fre-1),arrf(0:t*fre-1))
CALL dfftw_destroy_plan(p1)
call dfftw_plan_dft_1d(p1,t*fre,frecord(0,0:t*fre-1),trecord(0,0:t*fre-1),fftw_backward,fftw_measure)
!call dfftw_plan_dft_1d(p1,t*fre,frecord(0,0:t*fre-1),trecord(0,0:t*fre-1),fftw_backward,fftw_measure)

!position of sources and receivers
!source in even distribution 
i=0
!ang=2*pi
!DO nr=10,20,1
!  ang=ang-2*pi
!  DO WHILE(ang<=2*pi)
!    sourcexy(i,1)=nr
!    sourcexy(i,2)=ang
!    ang=ang+0.1
!    i=i+1
!  END DO 
!END DO

!source in random distribution
CALL RANDOM_SEED()
DO i=0,nsource-1
   CALL RANDOM_NUMBER(r)
   sourcexy(i,1)=1+4*r
END DO
!DO i=0,nsource-1
!   CALL RANDOM_NUMBER(ang)
!   sourcexy(i,2)=2*pi*ang
!END DO
!source in directional random distribution
!
!Inline
DO i=0,nsource/3-1
   CALL RANDOM_NUMBER(ang)
   sourcexy(i,2)=pi/6*ang-pi/12;
END DO
!offline
DO i=nsource/3,nsource-1
!DO i=0,nsource-1
   CALL RANDOM_NUMBER(ang)
!   sourcexy(i,2)=pi/6*ang+pi/4;
    sourcexy(i,2)=pi/6*ang+pi/12;
END DO
ang=0
!!receivers' xy
DO nr=0,2*nreceiver/5-1
   receiverxy(nr,1)=-0.12+0.005*nr
   receiverxy(nr,2)=0
END DO
DO nr=2*nreceiver/5,nreceiver-1
   receiverxy(nr,1)=-0.06+0.005*(nr-2*nreceiver/5)
   receiverxy(nr,2)=pi/2
END DO
!initialise whole simulated spectral record
DO nf=0,t*fre-1
  DO nr=0,nreceiver-1
    frecord(nr,nf)=CMPLX(0.0,0.0)
  END DO
END DO

!
!source loop
!$OMP PARALLEL DO PRIVATE(ns,nr,srt,r,f,nf,nt,npvf,dpf,dpv,c,k)
DO ns=0,nsource-1
  CALL RANDOM_NUMBER(x)
  srt=0.95*real(twindow)*x
!srt=0.95*real(twindow)*ns/nsource
!receiver loop
   DO nr=0,nreceiver-1
!receivers' xy
!     receiverxy(nr,1)=-0.063+0.005*nr
!     receiverxy(nr,2)=pi
!distance between receiver and source
     r=SQRT(sourcexy(ns,1)*sourcexy(ns,1)+receiverxy(nr,1)*receiverxy(nr,1)&
-2*sourcexy(ns,1)*receiverxy(nr,1)*COS(sourcexy(ns,2)-receiverxy(nr,2)))
      DO nf=0,t*fre-1
         f=real(nf)/real(t)
!impuse loop/ every timewindow /totally 60 shots
        DO nt=0,t/twindow-1
          IF(f<=45.AND.f>=0.1)THEN
!achieve phase velocity at f        
            DO npvf=0,npv-1
             IF(f<cv(0,npvf)) exit
            END DO
!interpert the phase velocity at frequency f which is between two certain frequencies
            dpf=cv(0,npvf)-cv(0,npvf-1)
            dpv=cv(1,npvf)-cv(1,npvf-1)
            c=(cv(1,npvf-1)+dpv/dpf*(f-cv(0,npvf-1)))/1000
            k=2*pi*f/c !
!            ang=-1*(receiverxy(nr,1)*receiverxy(nr,1)+r*r-sourcexy(ns,1)*sourcexy(ns,1))&
!/(2*receiverxy(nr,1)*r)
!            IF(abs(ang)>1) ang=ang/abs(ang)   
!            ang=SQRT(1-ang*ang)
!            IF(sourcexy(ns,2)>pi.AND.sourcexy(ns,2)<2*pi) ang=-1*ang
            frecord(nr,nf)=frecord(nr,nf)+arrf(nf)*EXP(-1*a*r)/SQRT(r)*CMPLX(COS(-srt*2*pi*f-k*r),&
SIN(-srt*2*pi*f-k*r))
          END IF
        END DO
      END DO
  END DO
END DO
!$OMP END PARALLEL DO

!stock data in sac files
DO nr=0,nreceiver-1
  CALL dfftw_execute_dft(p1,frecord(nr,0:t*fre-1),trecord(nr,0:t*fre-1))
!  WRITE(11,*) (real(trecord(nr,nt))/(t*fre),nt=0,t*fre-1)
  DO nt=0,t*fre-1
    sacdata(nt)=real(trecord(nr,nt))/(t*fre)
  END DO
!name sac file using order number in array.
  WRITE(addnumber,'(i2)') nr
  fname='H'//trim(adjustl(addnumber))//'.sac'
!write distance into the sac file
  distance=0
  CALL wsac(11,fname,sacdata,t*fre,dt,distance)
END DO
!
CALL dfftw_destroy_plan(p1)
DEALLOCATE(sourcexy,STAT=status)
DEALLOCATE(receiverxy,STAT=status)
DEALLOCATE(frecord,STAT=status)
DEALLOCATE(trecord,STAT=status)
DEALLOCATE(arrt,STAT=status)
DEALLOCATE(arrf,STAT=status)
DEALLOCATE(sacdata,STAT=status)
DEALLOCATE(cv,STAT=status)
CALL CPU_TIME(timee)
WRITE(*,*) 'TIME',timee-timeb
END PROGRAM znoiserayleigh

  SUBROUTINE wsac(fid,fname,head,NPTS,delta,distance)
    IMPLICIT NONE
    REAL,DIMENSION(0:NPTS-1)::head
    CHARACTER(len=128) :: fname
    REAL(KIND=4),DIMENSION(70):: head1
    INTEGER(KIND=4),DIMENSION(40) :: head2
    CHARACTER(len=8),DIMENSION(24) :: head3
    INTEGER :: row,col,fid,i
    INTEGER(KIND=4)::NPTS
    REAL(KIND=4)::delta,distance
    row=0
    col=0
    DO i=1,70
      head1(i)=-12345.
    END DO
    DO i=1,40
      head2(i)=-12345
    END DO
    DO i=1,24
      head3='-12345'
    END DO
!delta, b, e ,leven, npts, iftpe, nvhdr are necessary for sac file
    head1(1)=delta
    head1(6)=0!begin time
    head1(7)=(NPTS-1)*delta!end time
    head1(51)=distance
    head2(10)=NPTS
    head2(7)=6!nvhdr
    head2(16)=1!iftpe, data sections 
    head2(36)=1!leven, same time interval
    OPEN(fid, FILE=fname ,  ACCESS = 'Direct' , FORM = 'Unformatted' ,STATUS='REPLACE',recl=4*5)
    DO row=1,14
       WRITE( fid,rec=row) (head1((row-1)*5+col),col=1,5)
    END DO
    CLOSE(fid)
    OPEN(fid, FILE=fname ,  ACCESS = 'Direct' , FORM = 'Unformatted' ,STATUS='OLD',recl=4*5)
    DO row=15,22
       WRITE(fid,rec=row)  (head2((row-15)*5+col),col=1,5)
    END DO
    CLOSE(fid)
    OPEN(fid, FILE=fname ,  ACCESS = 'Direct' , FORM = 'Unformatted' ,STATUS='OLD',recl=8)
    DO i=1,24
       WRITE(fid,rec=i+55) head3(i)
    END DO
    CLOSE(fid)
    OPEN(fid, FILE=fname ,  ACCESS = 'Direct' , FORM = 'Unformatted' ,STATUS='OLD',recl=4)
    DO i=0,NPTS-1
       WRITE(fid,REC=i+158+1) head(i)
    END DO
    CLOSE(fid)
    END SUBROUTINE wsac
