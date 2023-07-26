c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:

      subroutine Calc2Bspinisospintrans(PiPhoto2Bx,PiPhoto2By,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,k,px,py,pz,ppx,ppy,ppz,calctype,Mnucl,verbosity)
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
c     Only care about neutral pion photoproduction 
      complex*16,intent(out) :: PiPhoto2By(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2Bx(0:1,-1:1,0:1,-1:1)

c     complex*16,intent(out) :: PiPhoto2Bxx(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Bxy(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Byx(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Byy(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Bx(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2By(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Bpx(0:1,-1:1,0:1,-1:1)
c     complex*16,intent(out) :: PiPhoto2Bpy(0:1,-1:1,0:1,-1:1)
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
      real*8 Mnucl
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p
      real*8  :: px,py,pz,ppx,ppy,ppz
               
      integer,intent(in) :: verbosity
      real*8 qpx,qpy,qpz,qppx,qppy,qppz,qx,qy,qz
      real*8 q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
      real*8 qsq,qpsq,qppsq,q12sq,qp12sq,qpp12sq
      real*8 qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z
      real*8 qpppsq,qppp12sq
      real*8 dl12by2
      real*8 factorA!,factorB,factorC,factorD,factorE
      real*8 factorAasy,factorC12,factorD12,factorE12
c     real*8 factorfg, factorfg12, factorfg2
c     real*8 factorfg212, factorhi, factorhi12
c     real*8 factorhi2, factorhi212        
c     real*8 factorjm, factorjm12, factorno, factorno12  
      real*8 K2n
      real*8 mPion
      real*8 q(3), kp(3)
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
c     PiPhoto2Bpx=c0
c     PiPhoto2Bpy=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':
c     
c     call calculateqs(qx,qy,qz,q12x,q12y,q12z,qpx,qpy,qpz,
c    &     qp12x,qp12y,qp12z,qppx,qppy,qppz,qpp12x,qpp12y,qpp12z,
c    &     qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z,
c    &     qsq,qpsq,qppsq,qpppsq,q12sq,qp12sq,qpp12sq,qppp12sq,px,py,pz,
c    &     ppx,ppy,ppz,
c    &     k,thetacm,verbosity)

      mPion=134.97
      write(*,*) "######################################################################"
      write(*,*) ""
      call calculateqsmass(px,py,pz,ppx,ppy,ppz,q,k,kp,thetacm,mPion,mNucl,verbosity)
      write(*,*) "In 2Bspinisospintrans.f DOT_PRODUCT(q,q)=", DOT_PRODUCT(q,q)
      write(*,*) ""
      write(*,*) "######################################################################"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OQ3 MEC contributions
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.  (mt12p .eq. 0)) then
c     
c     Here we calculate the non-isospin-changing part of the matrix element.
c     This is the only piece at OQ3.
c     
c     Define overall factors for spin-symmetric parts of matrix element
c     
         K2n=0.135*0.001 !in pion mass units
c   denom = (p12 - p12' +k_gamma/2)^2 =denomVec^2
c   denomVec = (k1-k2 + k2p-k1p+k_gamma)/2

c   DiagramA, Lenkewitz
c   ---------------
c   \frac{3}{2} * epsilon \cdot(\vec{\sigma}_1  + \vec{sigma_2}) (\vec{\tau_1} \cdot \vec{\tau_2} - \tau_1^z \tau_2^z)/denomVec^22
c   In the documentation they dont distinguish between k1 and k2
c   
c   In my derivation I just got 1/q^2 at threshold
        if(DOT_PRODUCT(q,q).le.0.0001) then
            write(*,*) "q is really small"
        end if 
        factorA=((-1)**t12)*1.5* K2n/(DOT_PRODUCT(q,q))

c     
c     Note that factorE, factorE12 only work if used in concert with factor B
************************************************************************************************
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram A, symmetric part
c     
c----------------------------------------------------------------------
            call CalcPionPhoto2BAx(PiPhoto2Bx,factorA,
     &           1.d0,0.d0,0.d0,s12p,s12,verbosity)
            call CalcPionPhoto2BAy(PiPhoto2By,factorA,
     &           0.d0,1.d0,0.d0,s12p,s12,verbosity)
c           call CalcPionPhoto2BAxx(PiPhoto2Bxx,factorA,
c    &           1.d0,0.d0,0.d0,dcos(thetacm),0.d0
c    &           ,-dsin(thetacm),s12p,s12,verbosity)
c           call CalcPionPhoto2BAxy(PiPhoto2Bxy,factorA,
c    &           1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
c    &           s12p,s12,verbosity)
c           call CalcPionPhoto2BAyx(PiPhoto2Byx,factorA,
c    &           0.d0,1.d0,0.d0,dcos(thetacm),0.d0,
c    &           -dsin(thetacm),s12p,s12,verbosity)
c           call CalcPionPhoto2BAyy(PiPhoto2Byy,factorA,
c    &           0.d0,1.d0,0.d0,0.d0,1.d0,0.d0,
c    &           s12p,s12,verbosity)
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram B, symmetric part
c     
c----------------------------------------------------------------------
c           call CalcCompton2BB(Comp2Bxx,Comp2Byy,
c    &           factorB,thetacm,qpx,qpy,qpz,qppx,qppy,qppz,
c    &           s12p,s12,verbosity)     
         else                   !l12-l12p is odd;  s12-s12p=+/- 1 => spin asymmetric part of operator.               
************************************************************************************************
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram A, anti-symmetric part
c     
c----------------------------------------------------------------------
            factorAasy=factorA
c           
            call CalcPionPhoto2BAxasy(PiPhoto2Bx,factorAasy,
     &           1.d0,0.d0,0.d0,s12p,s12,verbosity)
c           call CalcPionPhoto2BAxyasy(PiPhoto2Bxy,factorAasy,
c    &           1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
c    &           s12p,s12,verbosity)
c           call CalcPionPhoto2BAyxasy(PiPhoto2Byx,factorAasy,
c    &           0.d0,1.d0,0.d0,dcos(thetacm),0.d0,
c    &           -dsin(thetacm),s12p,s12,verbosity)
c           call CalcPionPhoto2BAyyasy(PiPhoto2Byy,factorAasy,
c    &           0.d0,1.d0,0.d0,0.d0,1.d0,0.d0,
c    &           s12p,s12,verbosity)
         end if                 ! s12 question
      end if                    !t12 question
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return 
      if (verbosity.eq.1000) continue
      end
