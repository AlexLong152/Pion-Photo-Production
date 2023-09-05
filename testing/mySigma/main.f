      program TEST
c     hold(Sp,Msp,S,Ms)
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      real*8 vec(3)
      integer i, Sp, S
      character(len=:), allocatable :: a
c     real*8 m, mp
      complex*16 holdOut, output
    
c     variable for line break
      a = "############################################################"
      hold=c0
      write(*,'(A)') "Evaluating (\vec{sigma}_1+\vec{sigma}_2).vec"
      vec=(/1,1,1/)
      write(*,'(A,I3,I3,I3)') "vec=",int(vec)

      do s=0,1
      do sp=0,1
          call singlesigma(hold,vec,sp,s)
          call holdOutput(hold,sp,s)
      end do!sp
      end do!s

      write(*,'(A)') a,a
      write(*,*) ""
      write(*,*) ""

      write(*,'(A)') "Evaluating (\vec{sigma}_1-\vec{sigma}_2).vec"
      write(*,*) ""
      vec=(/1,1,1/)
      write(*,'(A,I3,I3,I3)') "vec=",int(vec)

      do s=0,1
      do sp=0,1
          call singlesigmaasy(hold,vec,sp,s)
          call holdOutput(hold,sp,s)
      end do!sp
      end do!s

        write(*,'(A)') a,a
        write(*,*) ""
        write(*,'(A)') "Order of indicies in matrix listed here as m,mp"
        do mp=-1,1
        do m=-1,1
            write(*,'(I2,A1,I2)',advance="no") m,",",mp
            write(*,'(A3)',advance="no") " "
        end do
        write(*,*) ""
        end do
      end program TEST


      subroutine holdOutput(hold,sp,s)
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer sp,s
      logical printflag
      character(19) formt
      printflag=.FALSE.

      do m=-1,1
      do mp=-1,1
        if(hold(sp,mp,s,m).ne.cmplx(0.d0,0.d0)) then
            printflag=.TRUE.
        end if
      end do!mp
      end do!m

      formt= '(F0.3,SP,F0.3,"i")'

      if (printflag) then
        write(*,*) ""
        write(*,'(A5,I2,I2)') "s,sp=",s,sp
        do mp=-1,1,1
        do m=-1,1,1
            write(*,formt,advance="no") hold(sp,mp,s,m)
            write(*,'(A3)',advance="no") " "
        end do
        write(*,*) ""
        end do 
        write(*,*) ""
        write(*,*) ""
      end if


      end subroutine


      subroutine singlesigma(hold,vec,Sp,S)
c     calculates 2*S.A, where S=(sigma1+sigma2)/2
c     
c     OUTPUT VARIABLE:
c     
      complex*16 hold(0:1,-1:1,0:1,-1:1)
c     arguments for hold are (Sp,Msp,S,Ms)
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 Ax,Ay,Az
      real*8 vec(3)
      integer Sp,S
      complex*16 c0,ci
c     
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus
c     
c     
c********************************************************************
c     
      c0=cmplx(0.d0,0.d0)
      ci=cmplx(0.d0,1.d0)
      hold=c0
    
      Ax=vec(1)
      Ay=vec(2)
      Az=vec(3)
      if ((Sp .eq. 1) .and. (S .eq. 1)) then
         
         Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
         Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
         
         hold(1,1,1,1)=2.d0*Az
c     hold(1,0,1,0)=0
         hold(1,-1,1,-1)=-2.d0*Az
         hold(1,0,1,1)=-2.d0*Aplus
c     hold(1,-1,1,1)=0
         hold(1,1,1,0)=2.d0*Aminus   
         hold(1,-1,1,0)=-2.d0*Aplus !check       
c     hold(1,1,1,-1)=0       
         hold(1,0,1,-1)=2.d0*Aminus
      end if
      end subroutine

      subroutine singlesigmaasy(hold,vec,Sp,S)
c     
c     calculates (sigma1-sigma2).A
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
c********************************************************************
c     
c     OUTPUT VARIABLE:
c     
      complex*16 hold(0:1,-1:1,0:1,-1:1)
c     arguments for hold are (Sp,Msp,S,Ms)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 Ax,Ay,Az,factor, vec(3)
      integer verbosity
      integer Sp,S
c     
      complex*16 c0,ci
c********************************************************************
c     
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus
c     
c********************************************************************
c     
      factor=1.d0
      c0=cmplx(0.d0,0.d0)
      ci=cmplx(0.d0,1.d0)
      hold=c0
    
      Ax=vec(1)
      Ay=vec(2)
      Az=vec(3)

      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
      hold=c0
      
      if ((Sp .eq. 0) .and. (S .eq. 1)) then
         hold(0,0,1,1)=-2.d0*Aplus
         hold(0,0,1,0)=2.d0*Az
         hold(0,0,1,-1)=-2.d0*Aminus
      else if ((Sp .eq. 1) .and. (S .eq. 0)) then
         hold(1,1,0,0)=-2.d0*Aminus
         hold(1,0,0,0)=2.d0*Az
         hold(1,-1,0,0)=-2.d0*Aplus
      end if
      if (verbosity.eq.1000) continue
      end
