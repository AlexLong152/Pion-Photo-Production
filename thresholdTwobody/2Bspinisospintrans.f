c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine Calc2Bspinisospintrans(PiPhoto2Bx,PiPhoto2By,PiPhoto2Bz,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,k,px,py,pz,ppx,ppy,ppz,calctype,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Original version by Deepshikha, Spring 2006
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug 2016 modifications by Arman
c     -added OQ4 diagrams
c     -corrected overall T=0 vs. T=1 factor in symmetric part
c     -found missing factor of -2 in anti-symmetric part

c     List of changes to OQ3 calculation cf. Deepshikha's original:
c     (1) OQ3 diagrams must carry (-1)**(t21) in factor[ABCC12DD12] and factorAasy
c     (2) factors in asymmetric diagrams:
c     (a) additional (-2) in all asymmetric diagrams, i.e. including (1) above, it should be (-1)**(t21)
c     in diagram A
c     (b) factors of other diagrams more tricky since 1<->2 symmetry is more complicated there
c     -- correct version now specified with explicit "+"    and "-" signs; Deepshikha's original 
c     had the signs reversed (see point (a) above), except in diagram E
c     (c) [DP] in diagram E Deepshikha appears to've also missed the fact that 1<->2 leaves spin
c     structure unchanged, and so minus signs are not needed there.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie/DP Feb 2017: Added switch for OQ4 diagrams, corrected diagram E, commenting
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Feb 2017: modified order call such that OQ4 automatically calculates OQ3+OQ4 
c     note: readinput.f sets twobody at
c     Odelta0=OQ2 => ZERO: code exited inside readinput.f
c     Odelta2=OQ3 => calculated below: Beane diagrams
c     Odelta0=OQ3 => ZERO: code exited inside readinput.f
c     Odelta4=OQ4 => calculated after "if.calctype.eq.OQ3 return" below: Arman's diagrams
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Calls:
c     calculateqs: to calculate momenta
c     CalcCompton...: to calculate Compton amplitudes. Although much of the work is done
c     in those routines via "Calchold"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1)!x
      complex*16 Comp2Bxy(0:1,-1:1,0:1,-1:1)!y
      complex*16 Comp2Byx(0:1,-1:1,0:1,-1:1)!z
      complex*16 Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2By(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpy(0:1,-1:1,0:1,-1:1)

      complex*16,intent(out) :: PiPhoto2Bx(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2By(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2Bz(0:1,-1:1,0:1,-1:1)

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Note that Comp2Bab computes the amplitude for polarization a->polarization b
c     Indices in Comp2Bab are that first index gives NN spin state: S=0 or S=1,
c     second index gives spin projection. This is for final state. Third and fourth
c     indices give same for initial state. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      integer,intent(in) :: calctype
      real*8,intent(in)  :: thetacm,k
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p
      real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz
               
      integer,intent(in) :: verbosity
      
      real*8 qpx,qpy,qpz,qppx,qppy,qppz,qx,qy,qz
      real*8 q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
      real*8 qsq,qpsq,qppsq,q12sq,qp12sq,qpp12sq
      real*8 qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z
      real*8 qpppsq,qppp12sq
      real*8 dl12by2
      real*8 factorA,factorB
      real*8 factorA2,factorB2
      real*8 factorAasy,factorBasy
      real*8 denomVec(3), p(3), pp(3),kVec(3)

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Definitions of momenta repeated here for convenience
c     (All quantities in this comment to be read as vectors)
c     
c     q=p - p' + (k+k')/2: momentum of meson after photon 1 strikes, but
c     before photon 2 leaves
c     q'=p - p' + (k'-k)/2
c     q''=p - p' + (k-k')/2: momentum of meson after photon 1 strikes, and
c     after photon 2 leaves
c     
c     q12, q12', q12''=these quantities with p->-p and p'->-p'=>
c     q12''=-q' and q12'=-q''
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Factors:
c     factorA->factorE: for symmetric (in spin) part of OQ3 diagrams
c     factorAasy: for anti-symmetric (in spin) part of OQ3 diagram A
c     factorC12,factorD12,factorE12: factors associated with corresponding
c     diagrams after 1<->2 exchange=>momenta are "12" momenta cf. original
c     factors
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*************************************************************************************
c     
c     First a little initialization:
c     
      PiPhoto2Bx=c0
      PiPhoto2By=c0
      PiPhoto2Bz=c0

      Comp2Bxx=c0
      Comp2Bxy=c0
      Comp2Byx=c0
c     Comp2Byy=c0
c     Comp2Bx=c0
c     Comp2By=c0
c     Comp2Bpx=c0
c     Comp2Bpy=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':
c     
      call calculateqs(qx,qy,qz,q12x,q12y,q12z,qpx,qpy,qpz,
     &     qp12x,qp12y,qp12z,qppx,qppy,qppz,qpp12x,qpp12y,qpp12z,
     &     qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z,
     &     qsq,qpsq,qppsq,qpppsq,q12sq,qp12sq,qpp12sq,qppp12sq,px,py,pz,
     &     ppx,ppy,ppz,
     &     k,thetacm,verbosity)
c     p=(/px,py,pz/)
c     pp=(/ppx,ppy,ppz/)
c     kVec=(/0,0,k/)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OQ3 MEC contributions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     symmetric part! so that ResultXYXY in fm¹: twobody kernel should have units MeV^-4=fm/MeV^3
      factorA=  -(-1)**(t12)*(1.d0/((px-ppx)**2+(py-ppy)**2+(pz-ppz+k/2)**2))/HC
      factorB=+2*(-1)**(t12)*(1.d0/((px-ppx)**2+(py-ppy)**2+(pz-ppz+k/2)**2))*
     &     (1.d0/((px-ppx)**2+(py-ppy)**2+(pz-ppz-k/2)**2+mpi2))/HC

c     factorA2=  -(-1)**(t12)*(1.d0/(DOT_PRODUCT(p-pp+kVec/2,p-pp+kVec/2)))/HC
c     factorB2=+2*(-1)**(t12)*(1.d0/(DOT_PRODUCT(p-pp+kVec/2,p-pp+kVec/2)))*
c    &     (1.d0/(DOT_PRODUCT(p-pp+kVec/2,p-pp+kVec/2)+mpi2)/HC
c     if(factorA.ne.factorA2) then
c         write(*,*) "factorA!=factorA2"
c     end if
c     if(factorB.ne.factorB2) then
c         write(*,*) "factorB!=factorB2"
c     end if
c     antisymmetric part: turns out to be the same, only the vaue of t12 will be different
c     factorAasy=factorA
c     factorBasy=factorB
      
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcCompton2BA(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &           factorA,
     &           s12p,s12,verbosity)
            call CalcCompton2BB(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &           factorB,
     &           px-ppx,py-ppy,pz-ppz-k/2, ! preceding is vector dotted with σ
     &           px-ppx,py-ppy,pz-ppz, ! preceding is vector dotted with ε
     &           s12p,s12,verbosity)
c     
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only
c     
            call CalcCompton2BAasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &           factorAasy,
     &           s12p,s12,verbosity)
            call CalcCompton2BBasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &           factorBasy,
     &           px-ppx,py-ppy,pz-ppz-k/2, ! preceding is vector dotted with σ
     &           px-ppx,py-ppy,pz-ppz, ! preceding is vector dotted with ε
     &           s12p,s12,verbosity)

         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      PiPhoto2Bx=Comp2Bxx
      PiPhoto2By=Comp2Bxy
      PiPhoto2Bz=Comp2Byx

      if (verbosity.eq.1000) continue
      end
