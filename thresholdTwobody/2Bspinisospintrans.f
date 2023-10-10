c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     TODO: Get prefactors right
c     TODO: Get Diagram B asymmetric part working
      subroutine Calc2Bspinisospintrans(PiPhoto2Bx,PiPhoto2By,
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
      complex*16,intent(out) :: PiPhoto2By(0:1,-1:1,0:1,-1:1)
      complex*16 diff(0:1,-1:1,0:1,-1:1)
      complex*16,intent(out) :: PiPhoto2Bx(0:1,-1:1,0:1,-1:1)
c     

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
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':

      mPion=134.97
      call calculateqsmass(px,py,pz,ppx,ppy,ppz,q,k,q1,kp,thetacm,mPion,Mnucl,verbosity)

      kVec=(/0.d0,0.d0,k/)
      p=(/px,py,pz/)
      pp=(/ppx,ppy,ppz/)

c     TODO: get these from calculateqsmass instead of recalculating them
      k1=p-(kVec/2)
      k2=(-1*p)-(kVec/2)
      k1p=pp-kp/2
      k2p=(-1*pp)-kp/2

      p12=(k1-k2)/2
      p12p=(k1p-k2p)/2
c     write(*,*) "In 2Bspinisospintrans.f: (/px,py,pz)=",(/px,py,pz/) 
c     write(*,*) "In 2Bspinisospintrans.f: (/ppx,ppy,ppz)=",(/ppx,ppy,ppz/) 
c     p=(/px,py,pz/)
c     write(*,*) "In 2Bspinisospintrans.f: p-kVec/2=",p-kVec/2 
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
         K2n=(0.135*0.001/139.6)*197.3 !fm^2


c   
c     
************************************************************************************************
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram A, symmetric part
c     In my derivation I just got 1/q^2 at threshold
c     
c     DiagramA, Lenkewitz
c     ---------------
c     &F_{T/L}^{(a)} \vec{\epsilon}_{\lambda,\text{T/L}}\cdot\vec{S}=
c     &= \frac{3}{2}  \langle M_J'|
c     \frac{
c     \vec{\epsilon}_{\lambda,\text{T/L}}\cdot (\vec \sigma_1 + \vec \sigma_2)}{
c     \left(\vec{p}_{12}-\vec{p}_{12}^{\;\prime}+\vec{k}_\gamma/2\right)^2}}
c     | M_J\rangle
c     In the documentation they dont distinguish between k1 and k2
c     In my derivation (Alex Long) though I just got 1/q^2
c----------------------------------------------------------------------
c           factorA=((-1)**t12)*1.5* K2n/(DOT_PRODUCT(q,q))
            denomVec=p12-p12p+(kVec/2)

            factorA=((-1)**t12)*0.5/(DOT_PRODUCT(denomVec,denomVec))

            eps=(/1.d0,0.d0,0.d0/)
            call CalcPionPhoto2BA(PiPhoto2Bx,factorA,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,1.d0,0.d0/)
            call CalcPionPhoto2BA(PiPhoto2By,factorA,
     &           eps,s12p,s12,verbosity)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           The stuff below here is for debugging
c           diff=c0
c           eps=(/1.d0,0.d0,0.d0/)
c           call CalcPionPhoto2BA(diff,factorA,
c    &           eps,s12p,s12,verbosity)
c           call check_all_zero(allZerox,diff)

c           eps=(/0.d0,1.d0,0.d0/)
c           call CalcPionPhoto2BA(diff,factorA,
c    &           eps,s12p,s12,verbosity)
c           call check_all_zero(allZeroy,diff)

c           if ((.not.(allZerox)).or.(.not.(allZeroy))) then
c               write(*,'(A,I3,I3,I3,I3,I3,I3)') "mt12,mt12p,s12,s12p,t12,t12p=",
c    &                  mt12,mt12p,s12,s12p,t12,t12p
c               write(*,*) "x polarization"

c               eps=(/1.d0,0.d0,0.d0/)
c               call CalcPionPhoto2BA(diff,factorA,
c    &                eps,s12p,s12,verbosity)
c               call printDiff(diff/factorA)

c               eps=(/0.d0,1.d0,0.d0/)
c               call CalcPionPhoto2BA(diff,factorA,eps,s12p,s12,verbosity)
c               write(*,*) "y polarization"
c               call printDiff(diff/factorA)

c               write(*,*) ""
c               write(*,*) ""
c           end if
c           diff=c0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram B, symmetric part

c     Diagram B Lenkewtiz
c     -------------------------------------------------
c     F_{T/L}^{(b)} \vec{\epsilon}_{\lambda,\text{T/L}}\cdot\vec{S}_{M_J^\prime  M_J} =3 \langle M_J' |
c     \frac{
c     (\vec{p}_{12}-\vec{p}_{12}'-\vec{k}_\gamma/2)\cdot (\vec\sigma_1 + \vec\sigma_2)\vec{\epsilon}\cdot(\vec{p}_{12}-\vec{p}_{12}')
c     }
c     { \\denom of fraction
c     [(\vec{p}_{12}-\vec{p}_{12}^{\;\prime}- \vec{k}_\gamma/2 )^2+M_\pi^2 ]
c     \big[\vec{p}_{12}-\vec{p}_{12}^{\;\prime}+\vec{k}_\gamma/2 \big]^2
c     } 
c     |M_J\rangle_\psi,
c----------------------------------------------------------------------
c           eps=(/1.d0,0.d0,0.d0/)
c           factorB=K2n*((-1)**t12)*1.5*DOT_PRODUCT(eps,q1+q)/(
c    &          (DOT_PRODUCT(q1,q1)+mPion**2)
c    &          *(DOT_PRODUCT(q,q))
c    &          +mu)

c           call CalcPionPhoto2BB(PiPhoto2Bx,factorB,
c    &          q1,s12p,s12,verbosity)

c           eps=(/0.d0,1.d0,0.d0/)
c           factorB=K2n*((-1)**t12)*1.5*DOT_PRODUCT(eps,q1+q)/(
c    &          (DOT_PRODUCT(q1,q1)+mPion**2)
c    &          *(DOT_PRODUCT(q,q))
c    &          +mu)

c           call CalcPionPhoto2BB(PiPhoto2By,factorB,
c    &          q1,s12p,s12,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        BEGIN ASYMMETRIC PART
         else                   !l12-l12p is odd;  s12-s12p=+/- 1 => spin asymmetric part of operator.               
************************************************************************************************
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram A, anti-symmetric part
c     
c----------------------------------------------------------------------

            denomVec=p12-p12p+(kVec/2)
            factorAasy=((-1)**t12)*0.5/(DOT_PRODUCT(denomVec,denomVec))

            eps=(/1.d0,0.d0,0.d0/)
            call CalcPionPhoto2BAasy(PiPhoto2Bx,factorAasy,
     &           eps,s12p,s12,verbosity)

            eps=(/0.d0,1.d0,0.d0/)
            call CalcPionPhoto2BAasy(PiPhoto2By,factorAasy,
     &           eps,s12p,s12,verbosity)

c           diff=c0
c           call CalcPionPhoto2BAasy(diff,factorAasy,
c    &           eps,s12p,s12,verbosity)
c           call printDiff(diff)
c           diff=c0
c           call CalcPionPhoto2BAasy(diff,factorAasy,
c    &           eps,s12p,s12,verbosity)
c           call printDiff(diff)
c           diff=c0
c----------------------------------------------------------------------
c     
c     Calculate two-body diagram B, anti-symmetric part
c     
c----------------------------------------------------------------------

c           eps=(/1.d0,0.d0,0.d0/)
c           factorB=K2n*((-1)**t12)*1.5*DOT_PRODUCT(eps,q1+q)/(
c    &          (DOT_PRODUCT(q1,q1)+mPion**2)
c    &          *(DOT_PRODUCT(q,q))
c    &          +mu)

c           call CalcPionPhoto2BBasy(PiPhoto2Bx,factorB,
c    &           q1,s12p,s12,verbosity)

c           eps=(/0.d0,1.d0,0.d0/)
c           factorB=K2n*((-1)**t12)*1.5*DOT_PRODUCT(eps,q1+q)/(
c    &          (DOT_PRODUCT(q1,q1)+mPion**2)
c    &          *(DOT_PRODUCT(q,q))
c    &          +mu)

c           call CalcPionPhoto2BBasy(PiPhoto2By,factorB,
c    &           q1,s12p,s12,verbosity)
         end if                 ! s12 question
      end if                    !t12 question
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      if (verbosity.eq.1000) continue
      end

      subroutine printDiff(diff)
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




