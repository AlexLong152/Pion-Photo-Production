c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   TODO: rewrite in term of S_Dot_A instead of singlesigma
c     contains:
c              CalcPionPhoto2BAx
c              CalcPionPhoto2BAy
c              Calchold
c              singlesigma
c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug-Oct 2016/hgrie Feb 2017: Arman added OQ4 diagrams
c====================================================================

      subroutine CalcPionPhoto2BAx(PiPhoto2Bx,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for x->x.
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
      complex*16 PiPhoto2Bx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
      complex*16 PiPhoto2By(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 Ax,Ay,Az
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     singlesigma contains (sig_1+sig_2).A structure, symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            PiPhoto2Bx(Sp,Msp,S,Ms)=PiPhoto2Bx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
      
      subroutine CalcPionPhoto2BAy(PiPhoto2By,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for 
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
      complex*16 PiPhoto2Bx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
      complex*16 PiPhoto2By(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 Ax,Ay,Az
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     singlesigma contains (sig_1+sig_2).A structure, symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            PiPhoto2Bx(Sp,Msp,S,Ms)=PiPhoto2Bx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine CalcPionPhoto2BB(Pion2Bout,factor,
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
      real*8 eps(3), q1(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagator
c     Sp,S-final- and initial-state total spin of pair      
c     
c********************************************************************
c     
      call singlesigma(hold,q1(1),q2(2),q3(3),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Pion2By(Sp,Msp,S,Ms)=Pion2By(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcPionPhoto2BBx(Pion2Bxy,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
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
      complex*16 Pion2Bxy(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 Ax,Ay,Az
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagator
c     Sp,S-final- and initial-state total spin of pair      
c     
c********************************************************************
c     
      call singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Pion2Bxy(Sp,Msp,S,Ms)=Pion2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine Calchold(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)
c     
c     Calculates symmetric part of spin structure sig.A sig.B.
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
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 Ax,Ay,Az,Bx,By,Bz,factor,AdotB
      integer verbosity
      integer Sp,S
c     
c     factor-overall factor
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus,Bplus,Bminus
c     
c     
c********************************************************************
c     
      hold=c0
      AdotB=Ax*Bx+Ay*By+Az*Bz
      
      if ((Sp .eq. 0) .and. (S .eq. 0)) then
         hold(0,0,0,0)=-factor*2.d0*AdotB
      else if ((Sp .eq. 1) .and. (S .eq. 1)) then
         
         Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
         Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
         Bplus=-(Bx+ci*By)/(dsqrt(2.d0))
         Bminus=(Bx-ci*By)/(dsqrt(2.d0))
         
         hold(1,1,1,1)=factor*(2.d0*Az*Bz)
         hold(1,0,1,1)=-factor*2.d0*(Az*Bplus+Aplus*Bz)
         hold(1,-1,1,1)=factor*4.d0*Aplus*Bplus
         hold(1,1,1,0)=factor*2.d0*(Aminus*Bz+Az*Bminus)
         hold(1,0,1,0)=-2.d0*factor*(Aplus*Bminus+Aminus*Bplus+Az*Bz)
         hold(1,-1,1,0)=-hold(1,0,1,1)
         hold(1,1,1,-1)=factor*4.d0*Aminus*Bminus
         hold(1,0,1,-1)=-hold(1,1,1,0)
         hold(1,-1,1,-1)=hold(1,1,1,1)
      end if
      if (verbosity.eq.1000) continue
      end
c
      subroutine S_Dot_A(output,A,Sp,S,verbosity)
c     This is an interface to singlesigma thats more useful for my case
c********************************************************************
      complex*16 output(0:1,-1:1,0:1,-1:1)!same as hold
      real*8 factor
      real*8 A(3)
      integer verbosity
      integer Sp,S
      factor=1.d0
      subroutine(output,A(1),A(2),A(3),factor,Sp,S,verbosity)

      end subroutine

cccc  new subroutine to calculate the matrix elements of factor*A.S ccc
cccc same should be created in spintrickasy.f for factor*A.(sigma1-sigma2)      
      subroutine singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     TODO: why is "factor" involved here? Im pretty sure this actually calcualtes
c     TODO: factor*A.(sigma_1+sigma_2)
c     calculates 2*S.A, where S=(sigma1+sigma2)/2
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
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 Ax,Ay,Az,factor
      integer verbosity
      integer Sp,S
c     
c     factor-overall factor
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
      hold=c0
      
      if ((Sp .eq. 1) .and. (S .eq. 1)) then
         
         Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
         Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
         
         hold(1,1,1,1)=factor*2.d0*Az
c     hold(1,0,1,0)=0
         hold(1,-1,1,-1)=hold(1,1,1,1)
         hold(1,0,1,1)=-factor*2.d0*Aplus
c     hold(1,-1,1,1)=0
         hold(1,1,1,0)=factor*2.d0*Aminus   
         hold(1,-1,1,0)=-factor*2.d0*Aplus !check       
c     hold(1,1,1,-1)=0       
         hold(1,0,1,-1)=factor*2.d0*Aminus
      end if
      
      if (verbosity.eq.1000) continue
      end
c     
