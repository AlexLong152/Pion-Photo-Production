c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains:
c              
c              CalcCompton2BAxx
c              CalcCompton2BAxy
c              CalcCompton2BAyx
c              CalcCompton2BAyy
c              CalcCompton2BB
c              CalcCompton2BCx
c              CalcCompton2BCy
c              CalcCompton2BDx
c              CalcCompton2BDy
c              CalcCompton2BE
c              CalcCompton2Bfg
c              CalcCompton2Bfg2
c              CalcCompton2Bhi
c              CalcCompton2Bhi2
c              CalcCompton2Bjm
c              CalcCompton2Bno
c              CalcCompton2Bfgfg2ni
c              CalcCompton2Bhini
c              CalcCompton2Bhini2
c              CalcCompton2Baa2hihi2ni
c              CalcCompton2Bd2jmni
c              CalcCompton2Bdnoni
c              Calchold
c              singlesigma
c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug-Oct 2016/hgrie Feb 2017: Arman added OQ4 diagrams
c====================================================================

      subroutine CalcPionPhoto2BAx(Pion2Bxx,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for x->x.
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
      complex*16 Pion2Bxx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
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
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Pion2Bxx(Sp,Msp,S,Ms)=Pion2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine CalcPionPhoto2BAy(Comp2Bxy,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for x->y.
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
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
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
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcPionPhoto2BAyx(Comp2Byx,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for y->x.
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
      complex*16 Comp2Byx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 Ax,Ay,Az,Bx,By,Bz
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     eps.eps' structure
c     
c********************************************************************
c     
      call Calchold(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcPionPhoto2BAyy(Comp2Byy,factor,
     &     Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram A for x->y.
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
      complex*16 Comp2Byy(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor
      real*8 Ax,Ay,Az,Bx,By,Bz
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call Calchold(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine CalcCompton2BB(Comp2Bxx,Comp2Byy,factor,
     &     thetacm,Ax,Ay,Az,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram B. Note that this diagram
c     only contributes to two of the Compton polarization transitions.
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,polnfacx,polnfacy
      real*8 Ax,Ay,Az,Bx,By,Bz
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     thetacm-cm scattering angle, needed for poln vectors
c     
c********************************************************************
c     
      polnfacx=dcos(thetacm)
      polnfacy=1.d0
      call Calchold(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
     &           *polnfacx
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcCompton2BCx(Comp2Bxx,Comp2Bxy,factor,factor12,
     &     thetacm,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram C for x->x and x->y.
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Bxy(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x,
     &     polnfac12y
      real*8 qx,qy,qz,qppx,qppy,qppz
      real*8 q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for outgoing polarizations
c     x and y respectively
c     thetacm-cm scattering angle, needed for poln vectors
c     
c********************************************************************
c     
      polnfacx=((qx+qppx)*dcos(thetacm)-(qz+qppz)*dsin(thetacm))
      polnfacy=(qy+qppy)
      polnfac12x=((q12x+qpp12x)*dcos(thetacm)-(q12z+qpp12z)*
     &     dsin(thetacm))
      polnfac12y=(q12y+qpp12y)
      call Calchold(hold,1.d0,0.d0,0.d0,qppx,qppy,qppz,factor,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,1.d0,0.d0,0.d0,qpp12x,qpp12y,qpp12z,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcCompton2BCy(Comp2Byx,Comp2Byy,factor,factor12,
     &     thetacm,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram C for y->x and y->y.
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
      complex*16 Comp2Byx(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x,
     &     polnfac12y
      real*8 qx,qy,qz,qppx,qppy,qppz
      real*8 q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for outgoing polarizations
c     x and y respectively
c     thetacm-cm scattering angle, needed for poln vectors
c     
c********************************************************************
c     
      polnfacx=((qx+qppx)*dcos(thetacm)-(qz+qppz)*dsin(thetacm))
      polnfacy=(qy+qppy)
      polnfac12x=((q12x+qpp12x)*dcos(thetacm)-(q12z+qpp12z)*
     &     dsin(thetacm))
      polnfac12y=(q12y+qpp12y)
      call Calchold(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,
     &     factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,qpp12x,qpp12y,qpp12z,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine CalcCompton2BDx(Comp2Bxx,Comp2Byx,factor,factor12,
     &     thetacm,qx,qy,qpx,qpy,qpz,
     &     q12x,q12y,qp12x,qp12y,qp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram D for x->x and y->x.
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x,
     &     polnfac12y
      real*8 qx,qy,qpx,qpy,qpz
      real*8 q12x,q12y,qp12x,qp12y,qp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for incoming
c     polarizations x and y respectively
c     thetacm-cm scattering angle, needed for poln vectors
c     
c********************************************************************
c     
      polnfacx=qx+qpx
      polnfacy=qy+qpy
      polnfac12x=q12x+qp12x
      polnfac12y=q12y+qp12y
      call Calchold(hold,qpx,qpy,qpz,dcos(thetacm),0.d0,-dsin(thetacm),
     &     factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qp12x,qp12y,qp12z,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcCompton2BDy(Comp2Bxy,Comp2Byy,factor,factor12,
     &     qx,qy,qpx,qpy,qpz,
     &     q12x,q12y,qp12x,qp12y,qp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram D for x->y and y->y.
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
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 factor,factor12,polnfacx,polnfacy,polnfac12x,
     &     polnfac12y
      real*8 qx,qy,qpx,qpy,qpz
      real*8 q12x,q12y,qp12x,qp12y,qp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factor-contains meson propagators, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for incoming polarizations
c     x and y respectively
c     thetacm-cm scattering angle, needed for poln vectors      
c     
c********************************************************************
c     
      polnfacx=qx+qpx
      polnfacy=qy+qpy
      polnfac12x=q12x+qp12x
      polnfac12y=q12y+qp12y
      call Calchold(hold,qpx,qpy,qpz,0.d0,1.d0,0.d0,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qp12x,qp12y,qp12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c
c====================================================================
c
      subroutine CalcCompton2BE(Comp2Bxx,Comp2Bxy,Comp2Byx,Comp2Byy,
     &     factor,factE,factE12,thetacm,qx,qy,qz,qpx,qpy,qpz,qppx,qppy,
     &     qppz,q12x,q12y,q12z,qp12x,qp12y,qpp12x,qpp12y,qpp12z
     &     ,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates symmetric part of diagram E for all four polarization
c     transitions.
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
c     INPUT/OUTPUT VARIABLES:
c     
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Bxy(0:1,-1:1,0:1,-1:1)
     &     ,Comp2Byx(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
     &     ,hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factE,factE12
      real*8 qx,qy,qz,qpx,qpy,qpz,qppx,qppy,qppz
      real*8 q12x,q12y,q12z,qp12x,qp12y,qpp12x,qpp12y,qpp12z
      integer verbosity
c     
c     qx,qy,qz-components of vector q=p - p' + (k + k')/2
c     qpx,qpy,qpz-components of vector q'=p - p' + (k - k')/2
c     qppx,qppy,qppz-components of vector q''=p - p' + (k' - k)/2
c     q12x,q12y,q12z-components of vector q12=p' - p + (k + k')/2
c     hold-contains sig1.A sig.2B structure, symmetric part
c     factorE-contains meson propagators, overall factor
c     factorE12-same after 1<->2 interchange
c     thetacm-cm scattering angle, needed for poln vectors
c     
c********************************************************************
c     
c     LOCAL VARIABLES:
c     
      real*8 polnfacx,polnfacy,polnfacxp,polnfacyp
      real*8 polnfac12x,polnfac12y,polnfac12xp,polnfac12yp
      integer Ms,Msp,Sp,S
c     
c     polnfacx,polnfacy-initial-state polarization factors
c     polnfacxp,polnfacyp-final-state polarization factors
c     Sp,S-final- and initial-state total spin of pair      
c     
c**********************************************************************
c     
      polnfacxp=dcos(thetacm)*(qx + qppx) - dsin(thetacm)*(qz + qppz)
      polnfacx=qx + qpx
      polnfacyp=qy + qppy
      polnfacy=qy + qpy
      polnfac12xp=dcos(thetacm)*(q12x + qpp12x) - dsin(thetacm)*
     &     (q12z+qpp12z)
      polnfac12x=q12x + qp12x
      polnfac12yp=q12y + qpp12y
      polnfac12y=q12y + qp12y
c     if (firsttime) then
      call Calchold(hold,qpx,qpy,qpz,qppx,qppy,qppz,factor,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*
     &           (polnfacx*polnfacxp*factE+polnfac12x*polnfac12xp*factE12)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*
     &           (polnfacx*polnfacyp*factE+polnfac12x*polnfac12yp*factE12)  
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*
     &           (polnfacy*polnfacxp*factE+polnfac12y*polnfac12xp*factE12)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*
     &           (polnfacy*polnfacyp*factE+polnfac12y*polnfac12yp*factE12)
         end do
      end do
c     end if
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
ccccc Arman: Order-4 diagrams, spi-symmetric parts followccc
      subroutine CalcCompton2Bfg(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     Comp2Bpx,Comp2Bpy,factor,factor12,
     &     thetacm,k,qpppx,qpppy,qppx,qppy,qppz,
     &     qppp12x,qppp12y,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16  Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1),Comp2Bpy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qpppx,qpppy,qppx,qppy,qppz
      real*8 qppp12x,qppp12y,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qpppx
      polnfacy=qpppy
      polnfac12x=qppp12x
      polnfac12y=qppp12y
      
      
      call Calchold(hold,qppx,qppy,qppz,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      call Calchold(hold,qppx,qppy,qppz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     c these take care of singlesigma part of fg2 diagrams too
      call singlesigma(hold,qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
         end do
      end do
      call singlesigma(hold,qpp12x,qpp12y,qpp12z,factor12,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bfg2(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     Comp2Bx,Comp2By,factor,factor12,
     &     thetacm,qpppx,qpppy,qpppz,qppx,qppy,qppz,
     &     qppp12x,qppp12y,qppp12z,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16  Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bx(0:1,-1:1,0:1,-1:1),Comp2By(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8     polnfac12y
      real*8 qpppx,qpppy,qpppz,qppx,qppy,qppz
      real*8 qppp12x,qppp12y,qppp12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qpppx*dcos(thetacm)-qpppz*dsin(thetacm)
      polnfacy=qpppy
      polnfac12x=qppp12x*dcos(thetacm)-qppp12z*dsin(thetacm)
      polnfac12y=qppp12y
      call Calchold(hold,qppx,qppy,qppz,1.d0,0.d0,
     &     0.d0,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,1.d0,0.d0,
     &     0.d0,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      call Calchold(hold,qppx,qppy,qppz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bhi(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     Comp2Bpx,Comp2Bpy,factor,factor12,
     &     thetacm,k,qpppx,qpppy,qx,qy,qz,
     &     qppp12x,qppp12y,q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1),Comp2Bpy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qpppx,qpppy,qx,qy,qz
      real*8 qppp12x,qppp12y,q12x,q12y,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qpppx
      polnfacy=qpppy
      polnfac12x=qppp12x
      polnfac12y=qppp12y
      call Calchold(hold,qx,qy,qz,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      call Calchold(hold,qx,qy,qz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
ccccc factor=1 in the subroutine, then I put it here by hand, becauce (1)<->(2) does not change the spin structure.
      
      
      call singlesigma(hold,dcos(thetacm),0.d0,-dsin(thetacm),
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (factor*qy+factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-
     &           (factor*qx+factor12*q12x)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
cccccccc
      call singlesigma(hold,0.d0,1.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (factor*qy+factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           (factor*qx+factor12*q12x)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)  
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bhi2(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,Comp2Bx,Comp2By,factor,factor12,
     &     thetacm,k,qpppx,qpppy,qpppz,qx,qy,qz,
     &     qppp12x,qppp12y,qppp12z,q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bx(0:1,-1:1,0:1,-1:1),Comp2By(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qpppx,qpppy,qpppz,qx,qy,qz
      real*8 qppp12x,qppp12y,qppp12z,q12x,q12y,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qppp12x*dcos(thetacm)-qppp12z*dsin(thetacm)
      polnfacy=qppp12y
      polnfac12x=qpppx*dcos(thetacm)-qpppz*dsin(thetacm)
      polnfac12y=qpppy
      
      call Calchold(hold,qx,qy,qz,1.d0,0.d0,
     &     0.d0,factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,1.d0,0.d0,
     &     0.d0,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      call Calchold(hold,qx,qy,qz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call singlesigma(hold,1.d0,0.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (factor*qy+factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (factor*(qx*dcos(thetacm)-qz*dsin(thetacm))+
     &           factor12*(q12x*dcos(thetacm)-q12z*dsin(thetacm)))*
     &           k*ci*(1.d0+kappanu)*hold(Sp,Msp,S,Ms)
         end do
      end do
      
      call singlesigma(hold,0.d0,1.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           (factor*qy+factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
     &           (factor*(qx*dcos(thetacm)-qz*dsin(thetacm))+
     &           factor12*(q12x*dcos(thetacm)-q12z*dsin(thetacm)))*
     &           k*ci*(1.d0+kappanu)*hold(Sp,Msp,S,Ms)
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bjm(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,Comp2Bpx,Comp2Bpy,factor,factor12,
     &     thetacm,k,qx,qy,qz,qppx,qppy,qppz,qpppx,qpppy,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,qppp12x,qppp12y,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1),Comp2Bpy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacxx,polnfacxy,polnfacyx,polnfacyy
      real*8 polnfac12xx,polnfac12xy,polnfac12yx,polnfac12yy
      real*8 polnfacpx, polnfac12px,polnfacpy,polnfac12py
      real*8 qpppx,qpppy,qx,qy,qz,qppx,qppy,qppz,k
      real*8 qppp12x,qppp12y,q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      
      polnfacpx=qx*dcos(thetacm)-qz*dsin(thetacm)
      polnfacpy=qy
      
      polnfac12px=q12x*dcos(thetacm)-q12z*dsin(thetacm)
      polnfac12py=q12y
      
      polnfacxx=qpppx*polnfacpx
      polnfacxy=qpppx*qy
      polnfacyx=qpppy*polnfacpx
      polnfacyy=qpppy*qy
      polnfac12xx=qppp12x*polnfac12px
      polnfac12xy=qppp12x*q12y
      polnfac12yx=qppp12y*polnfac12px
      polnfac12yy=qppp12y*q12y
      
      call Calchold(hold,qx,qy,qz,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacxx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacxy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacyx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacyy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-polnfacpx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-polnfacpy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12xx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12xy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12yx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12yy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpx(Sp,Msp,S,Ms)=Comp2Bpx(Sp,Msp,S,Ms)-polnfac12px*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bpy(Sp,Msp,S,Ms)=Comp2Bpy(Sp,Msp,S,Ms)-polnfac12py*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call singlesigma(hold,qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qy*hold(Sp,Msp,S,Ms)*
     &           (qx*dcos(thetacm)-qz*dsin(thetacm))
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qy*hold(Sp,Msp,S,Ms)*
     &           qy
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qx*hold(Sp,Msp,S,Ms)*
     &           (qx*dcos(thetacm)-qz*dsin(thetacm))
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qx*hold(Sp,Msp,S,Ms)*
     &           qy
         end do
      end do
      call singlesigma(hold,qpp12x,qpp12y,qpp12z,factor12,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*q12y*hold(Sp,Msp,S,Ms)*
     &           (q12x*dcos(thetacm)-q12z*dsin(thetacm))
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*q12y*hold(Sp,Msp,S,Ms)*
     &           q12y
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*q12x*hold(Sp,Msp,S,Ms)*
     &           (q12x*dcos(thetacm)-q12z*dsin(thetacm))
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*q12x*hold(Sp,Msp,S,Ms)*
     &           qy
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bno(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,Comp2Bx,Comp2By,factor,factor12,
     &     thetacm,k,qx,qy,qz,qppx,qppy,qppz,qpppx,qpppy,qpppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,qppp12x,qppp12y,qppp12z,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bx(0:1,-1:1,0:1,-1:1),Comp2By(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacxx,polnfacxy,polnfacyx,polnfacyy
      real*8 polnfac12xx,polnfac12xy,polnfac12yx,polnfac12yy
      real*8 polnfacx,polnfacy,polnfac12x,polnfac12y
      real*8 qpppx,qpppy,qpppz,qx,qy,qz,qppx,qppy,qppz,k
      real*8 qppp12x,qppp12y,qppp12z,q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qppx
      polnfacy=qppy
      polnfac12x=qpp12x
      polnfac12y=qpp12y
      
      polnfacxx=qppx*(qpppx*dcos(thetacm)-qpppz*dsin(thetacm))
      polnfacxy=qppx*qpppy
      polnfacyx=qppy*polnfacxx/qppx
      polnfacyy=qppy*qpppy
      polnfac12xx=qpp12x*(qppp12x*dcos(thetacm)-qppp12z*dsin(thetacm))
      polnfac12xy=qpp12x*qppp12y
      polnfac12yx=qpp12y*polnfac12xx/qpp12x
      polnfac12yy=qpp12y*qppp12y
      
      
      call Calchold(hold,-q12x,-q12y,-q12z,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacxx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacxy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacyx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacyy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,-qx,-qy,-qz,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12xx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12xy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12yx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12yy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2By(Sp,Msp,S,Ms)=Comp2By(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call singlesigma(hold,qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qppx*hold(Sp,Msp,S,Ms)*
     &           q12y
            
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qppx*hold(Sp,Msp,S,Ms)*
     &           (q12x*dcos(thetacm)-q12z*dsin(thetacm))
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qppy*hold(Sp,Msp,S,Ms)*
     &           q12y
            
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qppy*hold(Sp,Msp,S,Ms)*
     &           (q12x*dcos(thetacm)-q12z*dsin(thetacm))
         end do
      end do
      call singlesigma(hold,qpp12x,qpp12y,qpp12z,factor12,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qpp12x*hold(Sp,Msp,S,Ms)*
     &           qy
            
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qpp12x*hold(Sp,Msp,S,Ms)*
     &           (qx*dcos(thetacm)-qz*dsin(thetacm))
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*qpp12y*hold(Sp,Msp,S,Ms)*
     &           qy
            
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           ci*(1.d0+kappanu)*k*qpp12y*hold(Sp,Msp,S,Ms)*
     &           (qx*dcos(thetacm)-qz*dsin(thetacm))
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
      
ccccccccccc new isospin structure contributions follow   ccccccccc
      subroutine CalcCompton2Bfgfg2ni(Comp2Bxx,Comp2Byx,
     &     Comp2Bxy,Comp2Byy,
     &     factor,factor12,
     &     thetacm,k,qppx,qppy,qppz,qx,qy,qz,
     &     qpp12x,qpp12y,qpp12z,q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     ``ni" stand for new isospin structure
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfac1x,polnfac1y,polnfac112x
      real*8 polnfac112y,k,polnfac2x,polnfac2y,polnfac212x
      real*8 polnfac212y,qppx,qppy,qppz,qx,qy,qz
      real*8 qpp12x,qpp12y,qpp12z,q12x,q12y,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      
cccccccccccccccccccccccccccccccccccccc
      call Calchold(hold,qppx,qppy,qppz,k*dsin(thetacm),0.d0,
     &     k*(1.d0+dcos(thetacm)),factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-dcos(thetacm)*
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,k*dsin(thetacm),0.d0,
     &     k*(1.d0+dcos(thetacm)),factor12,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+dcos(thetacm)*
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            
         end do
      end do
cccccccccccccccccccccccccccccccccccccc
      
      polnfac1x=qx+(1.d0+kappas)*k*dsin(thetacm)
      polnfac1y=qy
      polnfac112x=q12x+(1.d0+kappas)*k*dsin(thetacm)
      polnfac112y=q12y
      
      call Calchold(hold,qppx,qppy,qppz,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac1x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac1y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac112x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac112y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
c     
      call Calchold(hold,qppx,qppy,qppz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac1x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac1y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac112x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac112y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
ccccccccccccccccccccccccccccc
      polnfac2x=q12x*dcos(thetacm)-(q12z+(1.d0+kappas)*k)*dsin(thetacm)
      polnfac2y=q12y
      polnfac212x=qx*dcos(thetacm)-(qz+(1.d0+kappas)*k)*dsin(thetacm)
      polnfac212y=qy
      
      call Calchold(hold,qppx,qppy,qppz,1.d0,0.d0,
     &     0.d0,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac2x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac2y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,1.d0,0.d0,
     &     0.d0,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac212x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac212y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
c     
      call Calchold(hold,qppx,qppy,qppz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac2x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac2y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac212x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac212y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
ccccccc in (1)<->(2) isospin structure changes sign. It is put by hand in calculating the outputs.
      subroutine CalcCompton2Bhini(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     factor,factor12,
     &     thetacm,k,qx,qy,qz,
     &     q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qx,qy,qz
      real*8 q12x,q12y,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
c     c in this subroutine the sigma1.epsilonp sigma2.epsilon
c     c terms of diagrams hi and hi2 are included (in the new isospin parts).
c     c Simple modification of the variables polnfac and polnfac12 takes care of this
cccc  term in diagrams hi2.      
      polnfacx=-qx
      polnfacy=-qy
      polnfac12x=-q12x
      polnfac12y=-q12y
      
      
      call Calchold(hold,qx,qy,qz-k*(1.d0+kappas),dcos(thetacm),
     &     0.d0,-dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z-k*(1.d0+kappas),dcos(thetacm),
     &     0.d0,-dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
c     
      call Calchold(hold,qx,qy,qz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Bhini2(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qy,qz,
     &     q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qx,qy,qz
      real*8 q12x,q12y,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qx*dcos(thetacm)-qz*dsin(thetacm)
      polnfacy=qy
      polnfac12x=q12x*dcos(thetacm)-q12z*dsin(thetacm)
      polnfac12y=q12y
      
      call Calchold(hold,(1.d0+kappas)*k*dsin(thetacm)-qx,-qy,
     &     (1.d0+kappas)*k*dcos(thetacm)-qz,
     &     1.d0,0.d0,0.d0,factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,(1.d0+kappas)*k*dsin(thetacm)-q12x,-q12y,
     &     (1.d0+kappas)*k*dcos(thetacm)-q12z,1.d0,0.d0,0.d0,factor12,
     &     Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
c     
      call Calchold(hold,(1.d0+kappas)*k*dsin(thetacm)-qx,-qy,
     &     (1.d0+kappas)*k*dcos(thetacm)-qz,
     &     0.d0,1.d0,0.d0,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,(1.d0+kappas)*k*dsin(thetacm)-q12x,-q12y,
     &     (1.d0+kappas)*k*dcos(thetacm)-q12z,0.d0,1.d0,0.d0,
     &     factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
ccccccccccccccccccccccccccc
      subroutine CalcCompton2Baa2hihi2ni(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,
     &     factor,factor12,
     &     thetacm,k,qx,qz,
     &     q12x,q12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12
      real*8 polnfac,polnfac12,k
      real*8 qx,qz
      real*8 q12x,q12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
c     c in this subroutine the sigma1.epsilonp sigma2.epsilon
c     c terms of diagrams hi and hi2 are included (in the new isospin parts).
      
      
      
      polnfac=(1.d0+kappas)*k*((1.d0+dcos(thetacm))*qz+
     &     dsin(thetacm)*qx)-2.d0*k**2
      polnfac12=-(1.d0+kappas)*k*((1.d0+dcos(thetacm))*q12z+
     &     dsin(thetacm)*q12x)+2.d0*k**2
      
      
ccccccccccccccccccccccccccc
c     c what follows in this subroutine is the sigma1.epsilonp sigma2.epsilon
c     c terms of diagrams hi and hi2.
c     c  Also it includes diagrams a and a2. 
      
      
      
      call Calchold(hold,1.d0,0.d0,0.d0,dcos(thetacm),
     &     0.d0,-dsin(thetacm),1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (polnfac*factor+polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,dcos(thetacm),
     &     0.d0,-dsin(thetacm),1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           (polnfac*factor+polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,1.d0,0.d0,0.d0,0.d0,
     &     1.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (polnfac*factor+polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,0.d0,
     &     1.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
     &           (polnfac*factor+polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
ccccccccccccccccccccccccc      
      subroutine CalcCompton2Bd2jmni(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacxx,polnfacxy,polnfacyx,polnfacyy
      real*8 polnfac12xx,polnfac12xy,polnfac12yx,polnfac12yy
      real*8 polnfacpx, polnfacpy, polnfac12px, polnfac12py
      real*8 qx,qy,qz,qppx,qppy,qppz,k
      real*8 q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      
      polnfacpx=qx*dcos(thetacm)-qz*dsin(thetacm)
      polnfacpy=qy
      
      polnfac12px=q12x*dcos(thetacm)-q12z*dsin(thetacm)
      polnfac12py=q12y
      
      polnfacxx=qx*polnfacpx
      polnfacxy=qx*qy
      polnfacyx=qy*polnfacpx
      polnfacyy=qy*qy
      polnfac12xx=q12x*polnfac12px
      polnfac12xy=q12x*q12y
      polnfac12yx=q12y*polnfac12px
      polnfac12yy=q12y*q12y
      
      call Calchold(hold,qx,qy,qz+(1.d0+kappas)*k,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacxx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacxy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacyx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacyy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,q12x,q12y,q12z+(1.d0+kappas)*k,
     &     qpp12x,qpp12y,qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12xx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac12xy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12yx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac12yy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,1.d0,0.d0,0.d0,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacpx*
     &           (-(1.d0+kappas)*k*qz+k**2)*k*qz*hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacpy*
     &           (-(1.d0+kappas)*k*qz+k**2)*k*qz*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacpx*
     &           (-(1.d0+kappas)*k*qz+k**2)*hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacpy*
     &           (-(1.d0+kappas)*k*qz+k**2)*hold(Sp,Msp,S,Ms)
         end do
      end do
cccccccccccc
      call Calchold(hold,1.d0,0.d0,0.d0,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12px*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac12py*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12px*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac12py*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
cccccccccccccccccc
      subroutine CalcCompton2Bdnoni(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacxx,polnfacxy,polnfacyx,polnfacyy
      real*8 polnfac12xx,polnfac12xy,polnfac12yx,polnfac12yy
      real*8 polnfacx,polnfacy,polnfac12x,polnfac12y
      real*8 qx,qy,qz,qppx,qppy,qppz,k
      real*8 q12x,q12y,q12z,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig.q' sig.q'' structure
c     thetacm-c.m. scattering angle
c     factorB-contains meson propagators, overall factor
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=-qppx*(1.d0+kappas)*k*(q12x*dsin(thetacm)+
     &     q12z*dcos(thetacm))+qppx*k**2
      polnfacy=-qppy*(1.d0+kappas)*k*(q12x*dsin(thetacm)+
     &     q12z*dcos(thetacm))+qppy*k**2
      polnfac12x=-qpp12x*(1.d0+kappas)*k*(qx*dsin(thetacm)+
     &     qz*dcos(thetacm))+qpp12x*k**2
      polnfac12y=-qpp12y*(1.d0+kappas)*k*(qx*dsin(thetacm)+
     &     qz*dcos(thetacm))+qpp12y*k**2
      
      polnfacxx=qppx*(-q12x*dcos(thetacm)+q12z*dsin(thetacm))
      polnfacxy=-qppx*q12y
      polnfacyx=qppy*polnfacxx/qppx
      polnfacyy=-qppy*q12y
      polnfac12xx=qpp12x*(-qx*dcos(thetacm)+qz*dsin(thetacm))
      polnfac12xy=-qpp12x*qy
      polnfac12yx=qpp12y*polnfac12xx/qpp12x
      polnfac12yy=-qpp12y*qy
      
      
      
      call Calchold(hold,-(1.d0+kappas)*k*dsin(thetacm)-q12x,
     &     -q12y,-(1.d0+kappas)*k*dcos(thetacm)-q12z,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacxx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacxy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacyx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacyy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,-(1.d0+kappas)*k*dsin(thetacm)-qx,
     &     -qy,-(1.d0+kappas)*k*dcos(thetacm)-qz,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12xx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac12xy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12yx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-polnfac12yy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,dcos(thetacm),0.d0,-dsin(thetacm),qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
ccccccccccccccc
      call Calchold(hold,dcos(thetacm),0.d0,-dsin(thetacm),
     &     qpp12x,qpp12y,qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calchold(hold,0.d0,1.d0,0.d0,qpp12x,qpp12y,
     &     qpp12z,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
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
cccc  new subroutine to calculate the matrix elements of factor*A.S ccc
cccc same should be created in spintrickasy.f for factor*A.(sigma1-sigma2)      
      subroutine singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     
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
