c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains: CalcPionPhoto2BAx()
c               Calcholdasy
c               singlesigmaasy
c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug-Oct 2016/hgrie Feb 2017: Arman added OQ4 diagrams, corrected
c     isospin factor:
c     in subroutine Calcholdasy, all hold(..) must be multiplied by 2.d0
c     This is presently taken care of through the "factor" factors
c====================================================================
c     
      subroutine CalcPionPhoto2BAxasy(PiPhoto2Bx,factor,
     &     eps,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram A for x->x.
c     
c     Indices in Comp2Bab are that first index gives NN spin state:
c     S=0 or S=1, second index gives spin projection. This is for final
c     state. Third and fourth indices give same for initial state. 
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 PiPhoto2Bx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 eps(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     singlesigma contains (sig_1+sig_2).A structure, antisymmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair: 0 or 1
c     
c********************************************************************
c     
      call singlesigmaasy(hold,eps(1),eps(2),eps(3),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            PiPhoto2Bx(Sp,Msp,S,Ms)=PiPhoto2Bx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcPionPhoto2BAyasy(PiPhoto2Bx,factor,
     &     eps,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram A for x->x.
c     
c     Indices in Comp2Bab are that first index gives NN spin state:
c     S=0 or S=1, second index gives spin projection. This is for final
c     state. Third and fourth indices give same for initial state. 
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 PiPhoto2Bx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 eps(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     singlesigma contains (sig_1+sig_2).A structure, antisymmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair: 0 or 1
c     
c********************************************************************
c     
      call singlesigmaasy(hold,eps(1),eps(2),eps(3),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            PiPhoto2Bx(Sp,Msp,S,Ms)=PiPhoto2Bx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcPionPhoto2BBasy(Pion2Bout,factor,
     &     q1,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for x->y.
c     
c     Indices in Pion2Bab are that first index gives NN spin state:
c     S=0 or S=1, second index gives spin projection. This is for final
c     state. Third and fourth indices give same for initial state. 
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 Pion2Bout(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 q1(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagator
c     Sp,S-final- and initial-state total spin of pair      
c     
c********************************************************************
c     
      call singlesigmaasy(hold,q1(1),q1(2),q1(3),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Pion2Bout(Sp,Msp,S,Ms)=Pion2Bout(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
      subroutine Calcholdasy(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)   
c     
c     Calculates anti-symmetric part of spin structure sig1.A sig2.B.
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     OUTPUT VARIABLE:
      complex*16 hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 Ax,Ay,Az,Bx,By,Bz,factor
      integer Sp,S
      integer verbosity
c
c     factor-overall factor
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     
c     LOCAL VARIABLES:
      complex*16 Aplus,Aminus,Bplus,Bminus
c     
c********************************************************************
c     
      hold=c0
      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
      Bplus=-(Bx+ci*By)/(dsqrt(2.d0))
      Bminus=(Bx-ci*By)/(dsqrt(2.d0))
      if ((Sp .eq. 0) .and. (S .eq. 1)) then
         hold(0,0,1,1)=factor*2.d0*(-Aplus*Bz+Az*Bplus)
         hold(0,0,1,0)=-factor*2.d0*(Aplus*Bminus-Aminus*Bplus)
         hold(0,0,1,-1)=factor*2.d0*(Aminus*Bz-Az*Bminus)
      else if ((Sp .eq. 1) .and. (S .eq. 0)) then
         hold(1,1,0,0)=factor*2.d0*(Aminus*Bz-Az*Bminus)
         hold(1,0,0,0)=factor*2.d0*(Aplus*Bminus-Aminus*Bplus)
         hold(1,-1,0,0)=factor*2.d0*(-Aplus*Bz+Az*Bplus)
      end if
      
      if (verbosity.eq.1000) continue
      end
c
c********************************************************************
c********************************************************************
      subroutine singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     
c     calculates (sigma1-sigma2).A
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
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
      real*8 Ax,Ay,Az,factor
      integer verbosity
      integer Sp,S
c     
c********************************************************************
c     
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus
c     
c********************************************************************
c     
      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
      hold=c0
      
      if ((Sp .eq. 0) .and. (S .eq. 1)) then
         hold(0,0,1,1)=-factor*2.d0*Aplus
         hold(0,0,1,0)=factor*2.d0*Az
         hold(0,0,1,-1)=-factor*2.d0*Aminus
      else if ((Sp .eq. 1) .and. (S .eq. 0)) then
         hold(1,1,0,0)=-factor*2.d0*Aminus
         hold(1,0,0,0)=factor*2.d0*Az
         hold(1,-1,0,0)=-factor*2.d0*Aplus
      end if
      if (verbosity.eq.1000) continue
      end
c
