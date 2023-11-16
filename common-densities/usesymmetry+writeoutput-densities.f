c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c
c     hgrie May 2018:
c     merges old common/useparity.f and parity/output routines at the end of main.*.f
c     contains usesymmetryroutine() [was called useparityroutine()]
c              outputroutine()
c
c     twoSmax/twoMz dependence: via array sizes & do-loops

c     After rewrite for general-spin target: 
c     CAREFUL: CARTESIAN SYMMETRY ROUTINE PRETTY CERTAINLY WRONG!!!!!!!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Oct 2022: in usesymmetryroutine(): routine is applied now to only the MEs which are still undefined after 
c     transcarttosphere:
c     for Mzp=0, run only over Mz>0 
c     since for twoMi=twoMf=0, all 4 spherical amplitudes have already been constructed in transcarttospher().
c     
c     hgrie Aug 2020: added call to outputtomath() which produces
c     mathematica-friendly output of the first 2(2S+1)² MEs.
c      These are not related by time-reversal symmetry.
c     hgrie June 2018:
c     renamed from useparity* => usesymmetry*: Daniel showed that it's not parity but time-reversal invariance
c     see his notes pasted into "d Compton above Threshold I (2013-2016)", p. 46 Inv
c     This also fixed the factor as (-)^(Mf-Mi) for arbitrary target spin!
c     The routine remains the same, but all mention of "parity" is replaced by "symmetry".
c     
c     hgrie May 2018:
c     merges old common/useparity.f and parity/output routines at the end of main.*.f
c     contains useparityroutine()
c              outputroutine()
c
c     twoSmax/twoMz dependence: via array sizes & do-loops

c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c      
c     hgrie May 2018: more extensive description of changes in main.*.f
c     rewritten such that magnetic quantum numbers are "2*Mz" etc
c     changed phases of parity operation in useparityroutine() to cover 3He and deuteron,
c      
c     Implemented symmetry for arbitrary nucleon spin:
c     Use Mzp>=0, and for Mzp=0, run only over Mz>=0
c     -- that's still 2 more than necessary since ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time.
c     see manuscript "Compton Densities Approach" pp11-12
c*************************************************************************************
c*************************************************************************************
c*************************************************************************************
      subroutine usesymmetryroutine(amp11,amp1m1,ampm11,ampm1m1,twoSnucl,verbosity)
c
c     hgrie Oct 2022: see note above!
c      
c     Version 7 hgrie Feb 2013 -- NEW FILE
c     adapted from deuteron to 3He by hgrie Oct 2014  
c         
c     Use Symmetry to generate Compton matrix elements
c     see hgrie manuscript Deuteron Tensor pp68-74 and
c     notes "Implementing Parity and Time-Reversal Invariance", Feb 2013
c
c     UPDAE JUNE 2018: Symmetry is actually T-reversal -- see above for location of Daniel's notes!:
c     
c     implement A[-Mf,-lf;-Mi,-li] = (-)^(Mf-Mi) A[Mf,lf;Mi,li]
c      
c     see manuscript "Compton Densities Approach" pp. 8-10.
c     BUT SUPERSEDED BY notes pasted into "d Compton above Threshold I (2013-2016)", p. 46 Inv
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities lf, li (final, initial):
c     name of variable is ampλiλf, i.e. "IN" first, but the array is labelled (twoMzp,twoMz), i.e. "OUT" first!!! 
c     
      complex*16,intent(inout) :: amp11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)     ! 11:  (lf;li)=(1,1)
      complex*16,intent(inout) :: amp1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)    ! 1m1: (lf;li)=(-1,1)
      complex*16,intent(inout) :: ampm11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)    ! m11: (lf;li)=(1,-1)
      complex*16,intent(inout) :: ampm1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)   ! m1m1:(lf;li)=(-1,-1)
      integer,intent(in)       :: twoSnucl
      integer,intent(in)       :: verbosity

      integer                  :: twoMf,twoMi,twoMilimit
c     
c**********************************************************************
c     
      if (verbosity.eq.1000) continue
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      do twoMf=twoSnucl,0,-2
c     for Mzp=0, run only over Mz>0; for twoMi=twoMf=0, all 4 sherical amplitudes have already been constructed in transcarttospher()
         if (twoMf.eq.0) then
            twoMilimit = 2 ! !!!! see above
         else
            twoMilimit = -twoSnucl
         end if
         do twoMi=twoSnucl,twoMilimit,-2
            if (verbosity.eq.5) then 
               write(*,*) twoMf,"   ",twoMi
               write(*,*) amp11(twoMf,twoMi)
               write(*,*) ampm1m1(-twoMf,-twoMi)
               write(*,*) amp1m1(twoMf,twoMi)
               write(*,*) ampm11(-twoMf,-twoMi)
               write(*,*) ampm11(twoMf,twoMi)
               write(*,*) amp1m1(-twoMf,-twoMi)
               write(*,*) ampm1m1(twoMf,twoMi)
               write(*,*) amp11(-twoMf,-twoMi)
            end if              ! verbosity
            
            ampm1m1(-twoMf,-twoMi) = (-1)**((twoMf-twoMi)/2) * amp11(twoMf,twoMi)
            ampm11(-twoMf,-twoMi)  = (-1)**((twoMf-twoMi)/2) * amp1m1(twoMf,twoMi)
            amp1m1(-twoMf,-twoMi)  = (-1)**((twoMf-twoMi)/2) * ampm11(twoMf,twoMi)
            amp11(-twoMf,-twoMi)   = (-1)**((twoMf-twoMi)/2) * ampm1m1(twoMf,twoMi)
            
            if (verbosity.eq.5) then 
               write(*,*) "****************************************"
               write(*,*) amp11(twoMf,twoMi)
               write(*,*) ampm1m1(-twoMf,-twoMi)
               write(*,*) amp1m1(twoMf,twoMi)
               write(*,*) ampm11(-twoMf,-twoMi)
               write(*,*) ampm11(twoMf,twoMi)
               write(*,*) amp1m1(-twoMf,-twoMi)
               write(*,*) ampm1m1(twoMf,twoMi)
               write(*,*) amp11(-twoMf,-twoMi)
               write(*,*) "****************************************"
               write(*,*) "****************************************"
            end if              ! verbosity
         end do   
      end do   
      return
      end                       ! usesymmetryroutine
c
c*************************************************************************************
c*************************************************************************************
c*************************************************************************************

      subroutine outputroutine(outUnitno,cartesian,twoSnucl,twoMzplimit,
     &     Resultxx,Resultxy,Resultyx,Resultyy,verbosity)
c     hgrie May 2018: new routines, outsourced from main.*.f
c
c     construct symmetry-partners if necessary and outpout in cartesian or spherical coordinates
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities 
c     
      complex*16,intent(in) :: Resultxx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)     
      complex*16,intent(in) :: Resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resultyx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resultyy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in) :: twoSnucl,twoMzplimit
      integer,intent(in) :: outUnitno
      logical,intent(in) :: cartesian
      
      integer,intent(in) :: verbosity
      
      complex*16 Result11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)     
      complex*16 Result1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 Resultm11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 Resultm1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)

      integer :: twoMzp,twoMz,twoMzlimit
c     
c**********************************************************************
c     
      if (verbosity.eq.1000) continue
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     hgrie June/Oct 2014: use symmetry when specified, output to file
c
            write(*,*) "cartesian Q=",cartesian
            if (cartesian) then 
c     for output in CARTESIAN BASIS of photon polarisation   
c     let twoMzp only run over half of indices when symmetry is used: up to twoMzplimit
               do twoMzp=twoSnucl,twoMzplimit,-2
c     for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary since
c     ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time                  
                  if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                     twoMzlimit = 0
                  else
                     twoMzlimit = -twoSnucl
                  end if   
                  do twoMz=twoSnucl,twoMzlimit,-2
                     write (outUnitno,*) Resultxx(twoMzp,twoMz)
                     write (outUnitno,*) Resultxy(twoMzp,twoMz)
                     write (outUnitno,*) Resultyx(twoMzp,twoMz)
                     write (outUnitno,*) Resultyy(twoMzp,twoMz)
                  end do
               end do
c     hgrie June 2014:
c     if twoMzprime not trivial, use symmetry to get other MEs
c     hgrie May 2018: COULD STILL BE WRONG, AND DEFINITELY ONLY FOR 3HE -- NOT CHECKED WHEN CONVERTED to twoMz*!
               if ( twoMzplimit.eq.0 ) then
                  do twoMzp=twoSnucl,twoMzplimit,-2
                     do twoMz=twoSnucl,-twoSnucl,-2
                        write(*,*) "THIS RESULT UNCHECKED AFTER CONVERSION TO twoMz* -- LIKELY WRONG, DEFINITELY ONLY FOR 3He"
                        write (outUnitno,*) (-1)**((twoSnucl+1-twoMz)/2)*DCONJG(Resultxx(twoMzp,-twoMz)) ! "1" => any positive twoMzp, I think
                        write (outUnitno,*) (-1)**((twoSnucl+1-twoMz)/2)*DCONJG(Resultxy(twoMzp,-twoMz)) ! "1" => any positive twoMzp, I think
                        write (outUnitno,*) (-1)**((twoSnucl+1-twoMz)/2)*DCONJG(Resultyx(twoMzp,-twoMz)) ! "1" => any positive twoMzp, I think
                        write (outUnitno,*) (-1)**((twoSnucl+1-twoMz)/2)*DCONJG(Resultyy(twoMzp,-twoMz)) ! "1" => any positive twoMzp, I think
                     end do ! twoMz   
                  end do    ! twoMzp
               end if
c     for output in SPHERICAL BASIS of photon polarisation   
            else 
               call transcarttospher(Result11,Result1m1,Resultm11,Resultm1m1,
     &              twoSnucl,Resultxx,Resultxy,Resultyx,Resultyy,verbosity)
               if ( twoMzplimit.eq.0 ) then 
c     if twoMzprime not trivial, use symmetry to get other MEs
                  call usesymmetryroutine(Result11,Result1m1,Resultm11,Resultm1m1,twoSnucl,verbosity)
               end if   
c     output to file
               do twoMzp=twoSnucl,-twoSnucl,-2
                  do twoMz=twoSnucl,-twoSnucl,-2
                     write(outUnitno,*) Result11(twoMzp,twoMz)
                     write(outUnitno,*) Result1m1(twoMzp,twoMz)
                     write(outUnitno,*) Resultm11(twoMzp,twoMz)
                     write(outUnitno,*) Resultm1m1(twoMzp,twoMz)
                  end do
               end do
            end if
c     hgrie Aug 2020: if so wanted, output first independent MEs also to screen in a form that can directly be pasted into mathematica
            if (verbosity.ge.0) call outputtomath(Result11,Result1m1,Resultm11,Resultm1m1,twoSnucl,verbosity)
      return
      end                       ! outputroutine

      subroutine outputPiPhoto(outUnitno,cartesian,twoSnucl,twoMzplimit,
     &     Resultx,Resulty,Resultz,verbosity)
c     hgrie May 2018: new routines, outsourced from main.*.f
c
c     construct symmetry-partners if necessary and outpout in cartesian or spherical coordinates
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities 
c     
      complex*16,intent(in) :: Resultx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resulty(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resultz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer :: twoMzp,twoMz,twoMzlimit
      integer,intent(in) :: twoSnucl,twoMzplimit
      integer,intent(in) :: outUnitno
      logical,intent(in) :: cartesian
      
      integer,intent(in) :: verbosity
      
      character(len=25) ::str

      complex*16 :: sigmax(-1:1,-1:1)  ! (ms3p,ms3): sigma-x
      complex*16 :: sigmay(-1:1,-1:1) ! (ms3p,ms3): sigma-y
      complex*16 :: sigmaz(-1:1,-1:1)  ! (ms3p,ms3): sigma-z
      complex*16 :: SVec(3,-1:1,-1:1)
      character(39) formt
c     
c**********************************************************************
c     

       sigmax=0.d0
       sigmay=0.d0
       sigmaz=0.d0
      
       sigmax(1,-1)=dcmplx(1.d0,0)
       sigmax(-1,1)=dcmplx(1.d0,0)
 
       sigmay(1,-1)=dcmplx(0, -1.d0)
       sigmay(-1,1)=dcmplx (0, 1.d0) 
 
       sigmaz(1,1)=dcmplx(1.d0,0)
       sigmaz(-1,-1)=dcmplx(-1.d0,0)

       SVec(1,:,:)=sigmax
       SVec(2,:,:)=sigmay
       SVec(3,:,:)=sigmaz

c      formt= '(F0.3,SP,F0.3,"i")'
       formt='(F0.8,SP,F0.8,"i",A25,F0.3,SP,F0.3,"i")'
       formt='(F0.8,SP,F0.8,"i",A25,F0.0,SP,F0.0,"i")'
       write(outUnitno,*) ""
       write(*,*) ""

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(outUnitno,*) "x polarization - epsilon=<1,0,0>"
       write(*,*) "x polarization - epsilon=<1,0,0>"
       do twoMzp=twoSnucl,-twoSnucl,-2
          do twoMz=twoSnucl,-twoSnucl,-2
             write(outUnitno,formt) Resultx(twoMzp,twoMz),"sigmax(twoMzp,twoMz)=",sigmax(twoMzp,twoMz) !for (Mzp, Mz)=",twoMzp, twoMz
             write(*,formt) Resultx(twoMzp,twoMz),"sigmax(twoMzp,twoMz)=",sigmax(twoMzp,twoMz)
          end do
       end do
c      call printsum(outUnitno, sigmax, Resultx, twoSnucl,twoMzplimit)
       write(outUnitno,*) ""
       write(*,*) ""
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(outUnitno,*) "y polarization - epsilon=<0,1,0>"
       write(*,*) "y polarization - epsilon=<0,1,0>"
       do twoMzp=twoSnucl,-twoSnucl,-2
          do twoMz=twoSnucl,-twoSnucl,-2
             write(outUnitno,formt) Resulty(twoMzp,twoMz),"sigmay(twoMzp,twoMz)=",sigmay(twoMzp,twoMz) 
             write(*,formt) Resulty(twoMzp,twoMz),"sigmay(twoMzp,twoMz)=",sigmay(twoMzp,twoMz) 
          end do
       end do
c      call printsum(outUnitno, sigmay, Resulty, twoSnucl,twoMzplimit)

       write(outUnitno,*) ""
       write(*,*) ""
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       write(outUnitno,*) "z polarization - epsilon=<0,0,1>"
       write(*,*) "z polarization - epsilon=<0,0,1>"
       do twoMzp=twoSnucl,-twoSnucl,-2
          do twoMz=twoSnucl,-twoSnucl,-2
             write(outUnitno,formt) Resultz(twoMzp,twoMz),"sigmaz(twoMzp,twoMz)=",sigmaz(twoMzp,twoMz) 
             write(*,formt) Resultz(twoMzp,twoMz),"sigmaz(twoMzp,twoMz)=",sigmaz(twoMzp,twoMz) 
          end do
       end do

c      call printsum(outUnitno, sigmaz, Resultz, twoSnucl,twoMzplimit)
       write(outUnitno,*) ""
       write(*,*) ""
c     hgrie Aug 2020: if so wanted, output first independent MEs also to screen in a form that can directly be pasted into mathematica
            if (verbosity.ge.0) call outputtomathPiPhoto(Resultx,Resulty,Resultz, twoSnucl,verbosity)
      return
      end                       ! output PiPhoto

      subroutine printsum(outUnitno, sigma, resultxy, twoSnucl,twoMzplimit)

      implicit none
      integer,intent(in) :: twoSnucl,twoMzplimit
      integer,intent(in) :: outUnitno
      integer :: twoMzp,twoMz,twoMzlimit
      complex*16 :: sigma(-1:1,-1:1)
      complex*16 :: resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 total
      total=cmplx(0.d0,0.d0)

       do twoMzp=twoSnucl,-twoSnucl,-2
          do twoMz=twoSnucl,-twoSnucl,-2
             if (sigma(twoMzp,twoMz).ne.cmplx(0.d0,0.d0)) then
                 total=total+(resultxy(twoMzp,twoMz)/sigma(twoMzp,twoMz))
             end if
c            write(outUnitno,formt) Resultx(twoMzp,twoMz),"sigmax(twoMzp,twoMz)=",sigmax(twoMzp,twoMz) !, "for (Mzp, Mz)=",twoMzp, twoMz
c            write(*,formt) Resultx(twoMzp,twoMz),"sigmax(twoMzp,twoMz)=",sigmax(twoMzp,twoMz)
          end do
       end do
      if (AIMAG(total).ne.0.d0) then
        write(*,'(A5,F0.8,SP,F0.8,"i")') "sum=",total
        write(outUnitno,'(A4,F0.8,SP,F0.8,"i")') "sum=",total
      else
        write(*,'(A5,F0.8)') "sum=",Real(total)
        write(outUnitno,'(A4,F0.8)') "sum=",Real(total)
      end if
      end subroutine

      subroutine SingleDiagramOutput(cartesian,twoSnucl,twoMzplimit,
     &     DiagramX,label,verbosity)
c     hgrie May 2018: new routines, outsourced from main.*.f
c
c     construct symmetry-partners if necessary and outpout in cartesian or spherical coordinates
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities 
c     
      complex*16,intent(in) :: DiagramX(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      Character(len=40) :: label
      integer,intent(in) :: twoSnucl,twoMzplimit
      logical,intent(in) :: cartesian
      
      integer,intent(in) :: verbosity
      
      integer :: twoMzp,twoMz,twoMzlimit
      character(len=25) ::str
c     
c**********************************************************************
c     
c      label=trim(label)    
       write(*,*) ""

       write(*,*) label
       do twoMzp=twoSnucl,-twoSnucl,-2
          do twoMz=twoSnucl,-twoSnucl,-2
             write(*,*) DiagramX(twoMzp,twoMz)!,"for (Mzp, Mz)=",twoMzp,twoMz
          end do
       end do

       write(*,*) ""
            if (verbosity.ge.0) call outputtomathPiPhoto(Resultx,Resulty,twoSnucl,verbosity)
      return
      end                       ! SingleDiagramOutput
