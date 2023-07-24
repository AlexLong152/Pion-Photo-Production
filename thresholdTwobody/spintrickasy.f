c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains: CalcCompton2BAxxasy()
c               CalcCompton2BAxyasy()
c               CalcCompton2BAyxasy()
c               CalcCompton2BAyyasy()
c               CalcCompton2BBasy()
c               CalcCompton2BCxasy
c               CalcCompton2BCyasy
c               CalcCompton2BDxasy
c               CalcCompton2BDyasy
c               CalcCompton2BEasy
c               CalcCompton2BFgasy
c               CalcCompton2BFg2asy
c               CalcCompton2Bhiasy
c               CalcCompton2Bhi2asy
c               CalcCompton2Bjmasy
c               CalcCompton2Bnoasy
c               CalcCompton2Bfgfg2niasy
c               CalcCompton2Bhiniasy
c               CalcCompton2Bhini2asy
c               CalcCompton2Baa2hihi2niasy
c               CalcCompton2Bd2jmniasy
c               CalcCompton2Bdnoniasy
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
      subroutine CalcPionPhoto2BAxxasy(Comp2Bxx,factor,
     &     Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair: 0 or 1
c     
c********************************************************************
c     
      call singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do
c     
      if (verbosity.eq.1000) continue
      return
      end
c
c====================================================================
c
      subroutine CalcPionPhoto2BAxyasy(Comp2Bxy,factor,
     &     Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram A for x->y.
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
      real*8 Ax,Ay,Az,Bx,By,Bz
      integer Ms,Msp,Sp,S
      integer verbosity
c     
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
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
      subroutine CalcPionPhoto2BAyxasy(Comp2Byx,factor,
     &     Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram A for y->x.
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
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
      subroutine CalcPionPhoto2BAyyasy(Comp2Byy,factor,
     &     Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram A for y->y.
c     
c     Indices in Comp2Bab are that first index gives NN spin state:
c     S=0 or S=1, second index gives spin projection. This is for final
c     state. Third and fourth indices give same for initial state. c     
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     
c********************************************************************
c     
      call singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
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
      subroutine CalcCompton2BBasy(Comp2Bxx,Comp2Byy,factor,
     &     thetacm,Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram B. Note that this
c     diagram only contributesto two of the polarization transitions.
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
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),
     &     Comp2Byy(0:1,-1:1,0:1,-1:1),hold(0:1,-1:1,0:1,-1:1)
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     thetacm-c.m. scattering angle, needed for polarization factors
c     factor-contains meson propagator, overall factor
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for xx and yy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=dcos(thetacm)
      polnfacy=1.d0
      call Calcholdasy(hold,Ax,Ay,Az,Bx,By,Bz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           hold(Sp,Msp,S,Ms)*polnfacx
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
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
      subroutine CalcCompton2BCxasy(Comp2Bxx,Comp2Bxy,factor,factor12,
     &     thetacm,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram C for x->x and x->y.
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     thetacm-c.m. scattering angle, needed for polarization factors
c     factor-contains meson propagator, overall factor
c     factor12-same after 1<->2 interchange      
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for xx and xy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=((qx+qppx)*dcos(thetacm)-(qz+qppz)*dsin(thetacm))
      polnfacy=(qy+qppy)
      polnfac12x=((q12x+qpp12x)*dcos(thetacm)-(q12z+qpp12z)*
     &     dsin(thetacm))
      polnfac12y=(q12y+qpp12y)
      call Calcholdasy(hold,1.d0,0.d0,0.d0,qppx,qppy,qppz,factor,Sp,S,
     &     verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,1.d0,0.d0,0.d0,qpp12x,qpp12y,qpp12z,
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
      subroutine CalcCompton2BCyasy(Comp2Byx,Comp2Byy,factor,factor12,
     &     thetacm,qx,qy,qz,qppx,qppy,qppz,
     &     q12x,q12y,q12z,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram C for y->x and y->y.
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     thetacm-c.m. scattering angle, needed for polarization factors
c     factor-contains meson propagator, overall factor
c     factor12-same after 1<->2 interchange      
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for yx and yy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=((qx+qppx)*dcos(thetacm)-(qz+qppz)*dsin(thetacm))
      polnfacy=(qy+qppy)
      polnfac12x=((q12x+qpp12x)*dcos(thetacm)-(q12z+qpp12z)*
     &     dsin(thetacm))
      polnfac12y=(q12y+qpp12y)
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qpp12x,qpp12y,qpp12z,
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
      subroutine CalcCompton2BDxasy(Comp2Bxx,Comp2Byx,factor,factor12,
     &     thetacm,qx,qy,qpx,qpy,qpz,
     &     q12x,q12y,qp12x,qp12y,qp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram D for x->x and y->x.
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     thetacm-c.m. scattering angle, needed for polarization factors
c     factor-contains meson propagator, overall factor
c     factor12-same after 1<->2 interchange      
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for xx and yx from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qx+qpx
      polnfacy=qy+qpy
      polnfac12x=q12x+qp12x
      polnfac12y=q12y+qp12y
      call Calcholdasy(hold,qpx,qpy,qpz,dcos(thetacm),0.d0,
     &     -dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,qp12x,qp12y,qp12z,dcos(thetacm),0.d0,
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
      subroutine CalcCompton2BDyasy(Comp2Bxy,Comp2Byy,factor,factor12,
     &     qx,qy,qpx,qpy,qpz,
     &     q12x,q12y,qp12x,qp12y,qp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram D for x->y and y->y.
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
c     thetacm-c.m. scattering angle, needed for polarization factors
c     factor-contains meson propagator, overall factor
c     factor12-same after 1<->2 interchange
c     Sp,S-final- and initial-state total spin of pair
c     polnfacx,polnfacy-polarization factors for xy and yy from 
c     eps.eps' structure
c     
c********************************************************************
c     
      polnfacx=qx+qpx
      polnfacy=qy+qpy
      polnfac12x=q12x+qp12x
      polnfac12y=q12y+qp12y
      call Calcholdasy(hold,qpx,qpy,qpz,0.d0,1.d0,0.d0,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,qp12x,qp12y,qp12z,0.d0,1.d0,0.d0,
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
      subroutine CalcCompton2BEasy(Comp2Bxx,Comp2Bxy,Comp2Byx,Comp2Byy,
     &     factor,factE,factE12,thetacm,qx,qy,qz,qpx,qpy,qpz,qppx,qppy,
     &     qppz,q12x,q12y,q12z,qp12x,qp12y,qpp12x,qpp12y,qpp12z,Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates anti-symmetric part of diagram E for all four
c      polarization transitions.
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
     &     ,Comp2Byx(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1),
     &     hold(0:1,-1:1,0:1,-1:1)
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
c     thetacm-c.m. scattering angle
c     factor-contains two unchanged meson propagators, overall factor
c     factorE-contains third meson propagator
c     factorE12-same after 1<->2 interchange
c     qx,qy,qz-components of vector q=p - p' + (k + k')/2
c     qpx,qpy,qpz-components of vector q'=p - p' + (k - k')/2
c     qppx,qppy,qppz-components of vector q''=p - p' + (k' - k)/2
c     q12x,q12y,q12z-components of vector q12=p' - p + (k + k')/2
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
c     hold-contains sig1.A sig.2B structure, anti-symmetric part
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
      call Calcholdasy(hold,qpx,qpy,qpz,qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+ hold(Sp,Msp,S,Ms)*
     &           (polnfacx*polnfacxp*factE+polnfac12x*polnfac12xp*factE12)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+ hold(Sp,Msp,S,Ms)*
     &           (polnfacx*polnfacyp*factE+polnfac12x*polnfac12yp*factE12)  
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+ hold(Sp,Msp,S,Ms)*
     &           (polnfacy*polnfacxp*factE+polnfac12y*polnfac12xp*factE12)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+ hold(Sp,Msp,S,Ms)*
     &           (polnfacy*polnfacyp*factE+polnfac12y*polnfac12yp*factE12)
         end do
      end do
c     end if
      
      if (verbosity.eq.1000) continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Arman's Order-4 diagrams, spin-asymmetric parts follow 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine CalcCompton2Bfgasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
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
      
      
      call Calcholdasy(hold,dcos(thetacm),0.d0,
     &     -dsin(thetacm),qppx,qppy,qppz,factor,Sp,S,verbosity)
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,dcos(thetacm),0.d0,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,factor,Sp,
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
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
      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           ci*(1.d0+kappanu)*k*hold(Sp,Msp,S,Ms)*(1-dcos(thetacm))
         end do
      end do
      call singlesigmaasy(hold,qpp12x,qpp12y,qpp12z,factor12,
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
      subroutine CalcCompton2Bfg2asy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
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
      call Calcholdasy(hold,1.d0,0.d0,0.d0,
     &     qppx,qppy,qppz,factor,Sp,S,verbosity)
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,1.d0,0.d0,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,factor,Sp,
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
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
      subroutine CalcCompton2Bhiasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     Comp2Bpx,Comp2Bpy,factor,factor12,
     &     thetacm,k,qpppx,qpppy,qx,qy,qz,
     &     qppp12x,qppp12y,q12x,q12y,q12z,Sp,S,verbosity)
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
      call Calcholdasy(hold,qx,qy,qz,dcos(thetacm),0.d0,
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
      call Calcholdasy(hold,dcos(thetacm),0.d0,-dsin(thetacm),
     &     q12x,q12y,q12z,factor12,Sp,S,verbosity)
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
      call Calcholdasy(hold,qx,qy,qz,0.d0,1.d0,0.d0,factor,Sp,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,q12x,q12y,q12z,
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
ccccc factor=1 in the subrotine, then I put it here by hand, becauce (1)<->(2) does not change the spin structure.
      
      
      call singlesigmaasy(hold,-dcos(thetacm),0.d0,dsin(thetacm),
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (factor*qy-factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-
     &           (factor*qx-factor12*q12x)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
cccccccc
      call singlesigmaasy(hold,0.d0,-1.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (factor*qy-factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           (factor*qx-factor12*q12x)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)  
         end do
      end do
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bhi2asy(Comp2Bxx,Comp2Byx,Comp2Bxy,
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
      
      call Calcholdasy(hold,1.d0,0.d0,
     &     0.d0,qx,qy,qz,factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bx(Sp,Msp,S,Ms)=Comp2Bx(Sp,Msp,S,Ms)-hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,q12x,q12y,q12z,1.d0,0.d0,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qx,qy,qz,factor,Sp,
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
      call Calcholdasy(hold,q12x,q12y,q12z,0.d0,1.d0,0.d0,
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
ccccccccccccc
      call singlesigmaasy(hold,1.d0,0.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (factor*qy-factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (factor*(qx*dcos(thetacm)-qz*dsin(thetacm))-
     &           factor12*(q12x*dcos(thetacm)-q12z*dsin(thetacm)))*
     &           k*ci*(1.d0+kappanu)*hold(Sp,Msp,S,Ms)
         end do
      end do
      
      call singlesigmaasy(hold,0.d0,1.d0,0.d0,
     &     1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           (factor*qy-factor12*q12y)*k*ci*(1.d0+kappanu)*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
     &           (factor*(qx*dcos(thetacm)-qz*dsin(thetacm))-
     &           factor12*(q12x*dcos(thetacm)-q12z*dsin(thetacm)))*
     &           k*ci*(1.d0+kappanu)*hold(Sp,Msp,S,Ms)
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
      subroutine CalcCompton2Bjmasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
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
      
      call Calcholdasy(hold,qx,qy,qz,qppx,qppy,
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,
     &     q12x,q12y,q12z,factor12,Sp,S,verbosity)
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
      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,
     &     Sp,S,verbosity)
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
      call singlesigmaasy(hold,qpp12x,qpp12y,qpp12z,factor12,
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
      subroutine CalcCompton2Bnoasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
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
      
      
      call Calcholdasy(hold,-q12x,-q12y,-q12z,qppx,qppy,
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,
     &     -qx,-qy,-qz,factor12,Sp,S,verbosity)
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
      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)
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

      call singlesigmaasy(hold,qpp12x,qpp12y,qpp12z,factor12,
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
c********************************************************************
c********************************************************************
ccccccccccc new isospin structure contributions follow   ccccccccc
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Bfgfg2niasy(Comp2Bxx,Comp2Byx,
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
      call Calcholdasy(hold,k*dsin(thetacm),0.d0,k*(1.d0+dcos(thetacm)),
     &     qppx,qppy,qppz,factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-dcos(thetacm)*
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)-
     &           (1.d0+kappas)*hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,k*dsin(thetacm),
     &     0.d0,k*(1.d0+dcos(thetacm)),factor12,Sp,S,verbosity)
      
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
      
      call Calcholdasy(hold,dcos(thetacm),0.d0,-dsin(thetacm),
     &     qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac1x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac1y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,dcos(thetacm),0.d0,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac1x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac1y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
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
      
      call Calcholdasy(hold,1.d0,0.d0,0.d0,
     &     qppx,qppy,qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac2x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfac2y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,1.d0,0.d0,
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
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,qppz,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac2x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfac2y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,0.d0,1.d0,0.d0,
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
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Bhiniasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qy,qz,
     &     q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
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
      
      call Calcholdasy(hold,qx,qy,qz-k*(1.d0+kappas),dcos(thetacm),
     &     0.d0,-dsin(thetacm),factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,dcos(thetacm),0.d0,-dsin(thetacm),
     &     q12x,q12y,q12z-k*(1.d0+kappas),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfac12y*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
c     
      call Calcholdasy(hold,qx,qy,qz,0.d0,1.d0,0.d0,factor,Sp,
     &     S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,0.d0,1.d0,0.d0,q12x,q12y,q12z,
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
      subroutine CalcCompton2Bhini2asy(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qy,qz,
     &     q12x,q12y,q12z,Sp,S,verbosity)
c     
c********************************************************************
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
      
      call Calcholdasy(hold,1.d0,0.d0,0.d0,
     &     (1.d0+kappas)*k*dsin(thetacm)-qx,-qy,
     &     (1.d0+kappas)*k*dcos(thetacm)-qz,
     &     factor,Sp,S,verbosity)
      
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,(1.d0+kappas)*k*dsin(thetacm)-q12x,-q12y,
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
      call Calcholdasy(hold, 0.d0,1.d0,0.d0,
     &     (1.d0+kappas)*k*dsin(thetacm)-qx,-qy,
     &     (1.d0+kappas)*k*dcos(thetacm)-qz,
     &     factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+polnfacx*
     &           hold(Sp,Msp,S,Ms)
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+polnfacy*
     &           hold(Sp,Msp,S,Ms)
            
         end do
      end do
      call Calcholdasy(hold,(1.d0+kappas)*k*dsin(thetacm)-q12x,-q12y,
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
ccccccccccccccccccccccccc
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Baa2hihi2niasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
     &     Comp2Byy,factor,factor12,
     &     thetacm,k,qx,qz,
     &     q12x,q12z,Sp,S,verbosity)
c     
c********************************************************************
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
c     c Simple modification of the variables polnfac and polnfac12 takes care of this
cccc  term in diagrams hi2.      
      
      polnfac=(1.d0+kappas)*k*((1.d0+dcos(thetacm))*qz+
     &     dsin(thetacm)*qx)-2.d0*k**2
      polnfac12=-(1.d0+kappas)*k*((1.d0+dcos(thetacm))*q12z+
     &     dsin(thetacm)*q12x)+2.d0*k**2
      
ccccc polnfac and polnfac12 have a sign change of hi2 contributions vs the symmetric part, this is to compensate for sigma1.epsilon sigma2.epsilonp in hi vs sigma1.epsilonp sigma2.epsilon in hi2.     
      
ccccccccccccccccccccccccccc
c     c what follows in this subroutine is the sigma1.epsilon sigma2.epsilonp
c     c terms of diagrams hi and hi2.
c     c  Also it includes diagrams a and a2. 
      
      
      
      call Calcholdasy(hold,1.d0,0.d0,0.d0,dcos(thetacm),
     &     0.d0,-dsin(thetacm),1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+
     &           (polnfac*factor-polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,0.d0,1.d0,0.d0,dcos(thetacm),
     &     0.d0,-dsin(thetacm),1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+
     &           (polnfac*factor-polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,1.d0,0.d0,0.d0,0.d0,
     &     1.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+
     &           (polnfac*factor-polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,0.d0,1.d0,0.d0,0.d0,
     &     1.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byy(Sp,Msp,S,Ms)=Comp2Byy(Sp,Msp,S,Ms)+
     &           (polnfac*factor-polnfac12*factor12)*hold(Sp,Msp,S,Ms)
         end do
      end do
      
      if (verbosity.eq.1000) continue
      return
      end
ccccccccccccccccccccccccc
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Bd2jmniasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
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
      
      call Calcholdasy(hold,qx,qy,qz+(1.d0+kappas)*k,qppx,qppy,
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
      
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,
     &     q12x,q12y,q12z+(1.d0+kappas)*k,factor12,Sp,S,verbosity)
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
      call Calcholdasy(hold,1.d0,0.d0,0.d0,qppx,qppy,
     &     qppz,factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+polnfacpx*
     &           (-(1.d0+kappas)*k*qz+k**2)*hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+polnfacpy*
     &           (-(1.d0+kappas)*k*qz+k**2)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,0.d0,1.d0,0.d0,qppx,qppy,
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
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,
     &     1.d0,0.d0,0.d0,factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12px*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)-polnfac12py*
     &           (-(1.d0+kappas)*k*q12z+k**2)*hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calcholdasy(hold,qpp12x,qpp12y,qpp12z,
     &     0.d0,1.d0,0.d0,factor12,Sp,S,verbosity)
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
c********************************************************************
c********************************************************************
      subroutine CalcCompton2Bdnoniasy(Comp2Bxx,Comp2Byx,Comp2Bxy,
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
      call Calchold(hold,qpp12x,qpp12y,qpp12z,
     &     -(1.d0+kappas)*k*dsin(thetacm)-qx,
     &     -qy,-(1.d0+kappas)*k*dcos(thetacm)-qz,
     &     factor12,Sp,S,verbosity)
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
      call Calchold(hold,qpp12x,qpp12y,qpp12z,
     &     dcos(thetacm),0.d0,-dsin(thetacm),factor12,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)-polnfac12x*
     &           hold(Sp,Msp,S,Ms)
            
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)-polnfac12y*
     &           hold(Sp,Msp,S,Ms)
         end do
      end do
      call Calchold(hold,qpp12x,qpp12y,qpp12z,
     &     0.d0,1.d0,0.d0,factor12,Sp,S,verbosity)
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
c====================================================================
c====================================================================
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
