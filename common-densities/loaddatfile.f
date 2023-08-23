      subroutine getfilelength(multipolePath, nlines)
      IMPLICIT NONE
     
      character*500, intent(in) :: multipolePath
      integer, intent(out) :: nlines
      integer tmpUnit
      open(newunit=tmpUnit, file=multipolePath, status='OLD')!, iostat=test, action='read')
      nlines=0
      do
        read(tmpUnit,*,end=99)
        nlines = nlines + 1
      end do

 99    close(unit=tmpUnit)
      end subroutine

      subroutine loaddatfile(multipolePath, nlines, pole)
      IMPLICIT NONE
      character*500, intent(in) :: multipolePath
      integer, intent(in) ::  nlines
      real*8, intent(out) :: pole(2, nlines)
      integer iRead,unitNum
c     logical existCheck

c     inquire(file=trim(multipolePath), exist=existCheck )
c     write(*,*) "multipolePath=",trim(multipolePath), "  file_exists=", existCheck

      open(newunit=unitNum, file=trim(multipolePath), status='OLD')
c     write(*,*) "opened it"

      do iRead=1,nlines
        read (unitNum,*) pole(1,iRead), pole(2,iRead)
      end do

      close(unit=unitNum)
      end subroutine

      subroutine loadatSValue(multipolePath, sValue, poleValue)
      IMPLICIT NONE
     
      character*500, intent(in) :: multipolePath
      real*8 currentS, lastS
      real*8 currentPole, lastPole
c     note S is actualy sqrt(s) where s is the mandalstam variable
      real*8, intent(in):: sValue
      complex*16, intent(out):: poleValue
      real*8 interpolate
      integer tmpUnit

c     i=0
      currentS=0.d0
      lastS=0.d0
      currentS=0.d0
      lastS=0.d0
      open(newunit=tmpUnit, file=multipolePath, status='OLD')!, iostat=test, action='read')
      do
c       read (unitPole,*, end=98) pole(1,iRead), pole(2,iRead)
        lastPole=currentPole
        lastS=currentS
        read (tmpUnit,*, end=98) currentS, currentPole
        if ((lastS.le.sValue).and.(sValue.le.currentS)) then
            poleValue = interpolate(sValue, lastS, currentS, lastPole, currentPole)
            exit 
        end if
      end do
 98   close(unit=tmpUnit)
      end subroutine

      Function  interpolate(sValue, lastS, currentS, lastPole, currentPole)
c       Linear interpolation
c       y0,y1=last pole, current pole
c       x0,x1=last sqrt(mandalstam S), current sqrt(mandalstamS)
          real*8 :: interpolate 
          real*8, intent(in) :: sValue, lastS, currentS, lastPole, currentPole
          real*8 m, x, b
          m=(currentPole-lastPole)/(currentS-lastS)
          b=lastPole
          x=sValue-lastS
c         write(*,*) "sValue=", sValue 
c         write(*,*) "lastS=", lastS 
c         write(*,*) "currentS=", currentS 
c         write(*,*) "lastPole=", lastPole 
c         write(*,*) "currentPole=", currentPole 
c         write(*,*) "m=", m 
c         write(*,*) "b=", b 
c         write(*,*) "x=", x
c         write(*,*) ""
          interpolate=m*x+b
          return 
      End Function interpolate 
      
