c     Alex 2023: Pion Photoproduction 
c          Based on Lenkewitz derivation
c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine Calc2Bspinisospintrans(PiPhoto2Bx,PiPhoto2By,PiPhoto2Bz,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,k,px,py,pz,ppx,ppy,ppz,calctype,Mnucl,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Alex July 2023: 
c     Converted to pion photoproduction
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
      complex*16,intent(out) :: PiPhoto2Bx(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2By(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2Bz(0:1,-1:1,0:1,-1:1)
      complex*16 diff(0:1,-1:1,0:1,-1:1)

      complex*16 DiagramA(0:1,-1:1,0:1,-1:1)
      complex*16 DiagramB(0:1,-1:1,0:1,-1:1)
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
c     real*8 q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
c     real*8 qsq,qpsq,qppsq,q12sq,qp12sq,qpp12sq
c     real*8 qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z
c     real*8 qpppsq,qppp12sq
      real*8 dl12by2
      real*8 factorA,factorB!,factorC,factorD,factorE
      real*8 factorAasy, factorBasy
c     real*8 Bnumer
      real*8 K2n
      real*8 mPion
      real*8 q(3), q1(3)
      real*8 eps(3)
      real*8 k1(3),k1p(3),k2(3),k2p(3), p12(3), p12p(3), pp(3)
c     debugging variables
      real*8 p(3),tmp1(3), tmp2(3), kVec(3), kp(3)
      real*8 tmpVec(3)
      real*8 denomVec(3)
      real*8 mu
      logical :: allZerox, allZeroy
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
      mu=0.d0
      PiPhoto2Bx=c0
      PiPhoto2By=c0
      PiPhoto2Bz=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':

      mPion=134.97
      call calculateqsmass(px,py,pz,ppx,ppy,ppz,q,k,q1,kp,k1,k2,k1p,k2p,thetacm,mPion,mNucl,verbosity)
      kVec=(/0.d0,0.d0,k/)
c     p=(/px,py,pz/)
c     pp=(/ppx,ppy,ppz/)


      p12=(k1-k2)/2
      p12p=(k1p-k2p)/2
      K2n=(0.135*0.001/135)*197.3 !fm^2
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
c     139.6=pi^+ mass
c     1fm = 197.3 MeV^-1
c   
c     
************************************************************************************************
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only
c----------------------------------------------------------------------
c     Calculate two-body diagram A, symmetric part
c     In my derivation I just got 1/q^2 at threshold
c----------------------------------------------------------------------
            denomVec=p12-p12p+(kVec/2)!q
            eps=(/1.d0,0.d0,0.d0/)
            factorA=((-1)**t12)*K2n*0.5/(DOT_PRODUCT(denomVec,denomVec)+mu)
            call CalcPionPhoto2BA(PiPhoto2Bx,factorA,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,1.d0,0.d0/)
            call CalcPionPhoto2BA(PiPhoto2By,factorA,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,0.d0,1.d0/)
            call CalcPionPhoto2BA(PiPhoto2Bz,factorA,
     &           eps,s12p,s12,verbosity)

c----------------------------------------------------------------------
c     
c     Calculate two-body diagram B, symmetric part
c     subroutine CalcPionPhoto2BB(Pion2Bout,factor,
c    &     vec,Sp,S,verbosity)
c     Where vec is the vector being dotted into sigma_1+sigma_2
c----------------------------------------------------------------------
            eps=(/1.d0,0.d0,0.d0/)
            q1=p12-p12p-(kVec/2)
            q=p12-p12p+(kVec/2)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BB(PiPhoto2Bx,factorB,
     &          q1,s12p,s12,verbosity)

            eps=(/0.d0,1.d0,0.d0/)
            q1=p12-p12p-(kVec/2)
            q=p12-p12p+(kVec/2)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BB(PiPhoto2By,factorB,
     &          q1,s12p,s12,verbosity)

            eps=(/0.d0,0.d0,1.d0/)
            q1=p12-p12p-(kVec/2)
            q=p12-p12p+(kVec/2)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BB(PiPhoto2Bz,factorB,
     &          q1,s12p,s12,verbosity)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        BEGIN ASYMMETRIC PART
         else                   !l12-l12p is odd;  s12-s12p=+/- 1 => spin asymmetric part of operator.               
************************************************************************************************
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram A, anti-symmetric part
c     
c----------------------------------------------------------------------

            denomVec=p12-p12p+(kVec/2)!q
            eps=(/1.d0,0.d0,0.d0/)
            factorAasy=((-1)**t12)*K2n*0.5/(DOT_PRODUCT(denomVec,denomVec)+mu)
            call CalcPionPhoto2BAasy(PiPhoto2Bx,factorAasy,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,1.d0,0.d0/)
            call CalcPionPhoto2BAasy(PiPhoto2By,factorAasy,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,0.d0,1.d0/)
            call CalcPionPhoto2BAasy(PiPhoto2Bz,factorAasy,
     &           eps,s12p,s12,verbosity)
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram B, anti-symmetric part
c     recall q=q2
c----------------------------------------------------------------------

            q1=p12-p12p-(kVec/2)
            q=p12-p12p+(kVec/2)

            eps=(/1.d0,0.d0,0.d0/)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BBasy(PiPhoto2Bx,factorB,
     &          q1,s12p,s12,verbosity)


            eps=(/0.d0,1.d0,0.d0/)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BBasy(PiPhoto2By,factorB,
     &          q1,s12p,s12,verbosity)

            eps=(/0.d0,0.d0,1.d0/)
            factorB=K2n*((-1)**t12)*
     &          (DOT_PRODUCT(eps,p12-p12p))/(
     &          (DOT_PRODUCT(q1,q1)+mPion**2)*
     &          (DOT_PRODUCT(q,q))
     &             +mu)
            call CalcPionPhoto2BBasy(PiPhoto2Bz,factorB,
     &          q1,s12p,s12,verbosity)
         end if                 ! s12 question
      end if                    !t12 question
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      if (verbosity.eq.1000) continue
      end

      subroutine printDiff(diff)
c     function for debugging
      complex*16 diff(0:1,-1:1,0:1,-1:1)
      character(19) formt
      integer s,sp,m,mp
      formt = '(F0.4,SP,F0.4,"i")'
c     hold(Sp,Msp,S,Ms)
      do s=0,1
      do sp=0,1
      do m=-1,1,1
      do mp=-1,1,1
        if (diff(sp,mp,s,m).ne.c0) then
            write(*,'(A6)',advance='no') "diff="
            write(*,formt,advance='no') diff(sp,mp,s,m)
c           write(*,*) diff(sp,mp,s,m)
            write(*,'(A19,I3,I3,I3,I3)') " for (sp,mp,s,m) ",sp,mp,s,m
        end if
      end do
      end do
      end do
      end do
      end subroutine

      subroutine check_all_zero(allZero,diffs)
      logical :: allZero 
      complex*16 diffs(0:1,-1:1,0:1,-1:1)
      integer s,sp,m,mp
      allZero=.true.
      do s=0,1
      do sp=0,1
      do m=-1,1,1
      do mp=-1,1,1
      if (diffs(sp,mp,s,m).ne.c0) then
        allZero=.false.
      end if
      end do
      end do
      end do
      end do
      end subroutine




