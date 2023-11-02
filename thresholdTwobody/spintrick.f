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

      subroutine CalcPionPhoto2BA(PiPhoto2B,factor,
     &     eps,Sp,S,verbosity)
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
      complex*16 PiPhoto2B(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     complex*16 PiPhoto2By(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
c     real*8 Ax,Ay,Az
      real*8 eps(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     singlesigma contains (sig_1+sig_2).A structure, symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     

      call S_Dot_A(hold,factor,eps,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            PiPhoto2B(Sp,Msp,S,Ms)=PiPhoto2B(Sp,Msp,S,Ms)+factor*hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
      
c
      subroutine CalcPionPhoto2BB(Pion2Bout,factor,
     &     vec,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram B.
c     "vec" gets dotted into sigma_1+sigma_2
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
      real*8 vec(3)
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     Sp,S-final- and initial-state total spin of pair      
c     
c********************************************************************
c     
      call singlesigma(hold,vec(1),vec(2),vec(3),factor,Sp,S,verbosity)
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
      
      subroutine S_Dot_A(output,factor,A,Sp,S,verbosity)
c     This is an interface to singlesigma thats more useful for my case
c********************************************************************
      complex*16 output(0:1,-1:1,0:1,-1:1)!same as hold
      real*8 factor
      real*8 A(3)
      integer verbosity
      integer Sp,S
c     factor=1.d0
      call singlesigma(output,A(1),A(2),A(3),factor,Sp,S,verbosity)
      end subroutine


cccc  new subroutine to calculate the matrix elements of factor*A.S ccc
cccc same should be created in spintrickasy.f for factor*A.(sigma1-sigma2)      

      subroutine singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
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
c     arguments for hold are (Sp,Msp,S,Ms)
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
c     why only Sp=1=S case ??? 

      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))

      if ((Sp .eq. 1) .and. (S .eq. 1)) then
          hold(1,1,1,1)=factor*Az
          hold(1,-1,1,-1)=-1.d0*factor*Az
          hold(1,0,1,1)=-1.d0*Aplus
          hold(1,1,1,0)=factor*Aminus   
          hold(1,-1,1,0)=-factor*Aplus
          hold(1,0,1,-1)=factor*Aminus
      end if

      if ((Sp .eq. 1) .and. (S .eq. 0)) then
          hold(0,0,1,-1)=factor*Aplus
          hold(0,0,1,0)=factor*Az
          hold(0,0,1,1)=factor*Aminus
      end if

      if ((Sp .eq. 0) .and. (S .eq. 1)) then
          hold(1,-1,0,0)=factor*Aminus
          hold(1,0,0,0)=factor*Az
          hold(1,1,0,0)=factor*Aplus
      end if
      if (verbosity.eq.1000) continue
      end
