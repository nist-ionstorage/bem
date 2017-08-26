c
c This a driver to test fastlap V2.0 as a callable procedure from fortran
c for computing a desingularized solution of a single time step of 
c Steve Scorpio's free surface problem.
c Input is from Steve's dumped data.
c
c NB you can precondition these problems with OL which has been extended to 
c handle non-square local problems, but it is not very effective.  Use the SP
c preconditioner which was designed for the desingularized case.
c
c Written in Modern Fortran by Korsmeyer, winter 1993.  (c) MIT 1993.
c
      implicit none
      integer bsize, fssize, size
      parameter (bsize=510, fssize=2486)
      parameter (size=bsize+fssize)
      double precision a(size,size), r, rx, ry, rz
      integer ipvt(size), info, idum
      double precision x(3,4,size), xnrm(3,size),
     $     xcoll(3,size), lhsvect(size), rhsvect(size) 
      integer type(size), shape(size), dtype(size), job, fljob
      integer rhstype(size), lhstype(size), rhsindex(4,size),
     $     lhsindex(4,size)
      integer NULL, POINT_SOURCE, CONSTANT_SOURCE, CONSTANT_DIPOLE, 
     $     LINEAR_SOURCE, LINEAR_DIPOLE
      data NULL, POINT_SOURCE, CONSTANT_SOURCE, CONSTANT_DIPOLE, 
     $     LINEAR_SOURCE, LINEAR_DIPOLE /0, 1, 11, 12, 21, 22/
      double precision scorp_sol(size)
      integer nummom, numlev, i, j, k, nlhs, nrhs, nsing
      double precision error, max_a_error, ave_a_error,
     $     max_r_error, ave_r_error
      character*1 dumchar
      integer fastlap
      external fastlap
      double precision tol
      integer numitr, maxitr
      data maxitr /1024/
c
c Query for run parameters and data file name.
c
c
      write(*,'(//''We solve a single time step of the Scorpio''/
     $     ''free-sruface problem  '')')
      write(*,'(//''You should use the SP preconditioner for''/
     $     ''this problem.  Set PRECOND = SP in mulGlobal.h. '')')
      write(*,'(''Select expansion order and tree depth:  ''$)')
      read(*,*) nummom, numlev
      write(*,'(''Select the iteration tolerance:  ''$)')
      read(*,*) tol
c
c Read desingularized source, collocation points, etc.
c (There are 510 body points and 2486 free-surface points.)
c
      open(unit=20,file='bdbvp-1.dat',status='old')
      open(unit=21,file='bdbvp-2.dat',status='old')
      open(unit=22,file='fsbvp.dat',status='old')
c
c Get past the headers.
c
      do i = 1,3
         read(20,'(a)') dumchar
         read(21,'(a)') dumchar
         read(22,'(a)') dumchar
      end do
      do i = 1,4
         read(20,'(a)') dumchar
         read(22,'(a)') dumchar
      end do
      print *, ' Reading',bsize+fssize,' source and field points.' 
      nsing = size
      nlhs = nsing
      nrhs = nsing
      do i = 1, bsize
         read(20,*) (xcoll(k,i),k=1,3),(x(k,1,i),k=1,3),
     $        rhsvect(i),scorp_sol(i)
         rhsindex(1,i) = i
         lhsindex(1,i) = i
         shape(i) = 1
         dtype(i) = 1
         rhstype(i) = POINT_SOURCE
         lhstype(i) = POINT_SOURCE
         read(21,*) (xnrm(k,i),k=1,3)
      end do
      do i = bsize+1, size
         read(22,*) (xcoll(k,i),k=1,3),(x(k,1,i),k=1,3),
     $        rhsvect(i),scorp_sol(i)
         rhsindex(1,i) = i
         lhsindex(1,i) = i
         shape(i) = 1
         dtype(i) = 0
         rhstype(i) = POINT_SOURCE
         lhstype(i) = POINT_SOURCE
      end do
c
c Call fastlap to compute a solution or field at the vector of field points.
c
      fljob = 2
      print *, "Calling fastlap."
      numitr = fastlap(nlhs,nrhs,nsing,x,shape,dtype,lhstype,rhstype,
     $     lhsindex,rhsindex,lhsvect,rhsvect,xcoll,xnrm,numlev,nummom,
     $     maxitr,tol,fljob)
      write(*,'(//i3,
     $     '' iterations knocked down residual to:'',e14.8)')
     $     numitr, tol
c
c Compute the solution by LU, and call it exact.  Or possibly 
c go around this block and read a previous solution.
c
      go to 100
      do j = 1,size
         do i = 1,bsize
            r = (xcoll(1,i) - x(1,1,j))**2
     $           + (xcoll(2,i) - x(2,1,j))**2
     $           + (xcoll(3,i) - x(3,1,j))**2
            rx = -(xcoll(1,i) - x(1,1,j)) / r**1.5
            ry = -(xcoll(2,i) - x(2,1,j)) / r**1.5
            rz = -(xcoll(3,i) - x(3,1,j)) / r**1.5
            a(i,j) = xnrm(1,i) * rx + xnrm(2,i) * ry + xnrm(3,i) * rz
         end do
         do i = bsize+1,size
            r = (xcoll(1,i) - x(1,1,j))**2
     $           + (xcoll(2,i) - x(2,1,j))**2
     $           + (xcoll(3,i) - x(3,1,j))**2
            a(i,j) = 1./sqrt(r)
         end do
      end do
      call DGETRF(size,size,a,size,ipvt,INFO)
      call DGETRS('N',size,1,a,size,ipvt,
     $     rhsvect,size,INFO)
 100  continue
      open(unit=24,file='oldLU.dat',status='old',form='formatted')
      do i = 1,size
         read(24,*) idum, rhsvect(i)
      end do
      close(24)
c
c Compute the error in the multi and the scorp solution.
c
      
      ave_a_error = 0.
      ave_r_error = 0.
      max_a_error = 0.
      max_r_error = 0.
      open(unit=25,file='solutions.dat',
     $     status='unknown',form='formatted')
      write(25,'(''Index   LU solution   Dom-Decomp   fastlap'')')
      do i = 1,nlhs
         write(25,999) i,rhsvect(i),scorp_sol(i),lhsvect(i)
         error = abs(rhsvect(i) - lhsvect(i))
         write(98,999) i,error,error/rhsvect(i)
         ave_a_error = ave_a_error + error
         ave_r_error = ave_r_error + error/rhsvect(i)
         max_a_error = max(max_a_error,error)
         max_r_error = max(max_r_error,error/rhsvect(i))
         if (error/rhsvect(i) .gt. 100) write(88,*) 'multi:  ',
     $        i, rhsvect(i), lhsvect(i)
      end do
 999  format(i6,3e20.10)
      ave_a_error = ave_a_error / float(nlhs)
      ave_r_error = ave_r_error / float(nlhs)
      write(*,'(//''Average absolute multi error ='',
     $     f13.8)') ave_a_error
      write(*,'(''Maximum absolute multi error ='',
     $     f13.8)') max_a_error
      write(*,'(//''Average relative multi error ='',
     $     f13.8)') ave_r_error
      write(*,'(''Maximum relative multi error ='',
     $     f13.8)') max_r_error
      ave_a_error = 0.
      ave_r_error = 0.
      max_a_error = 0.
      max_r_error = 0.
      do i = 1,nlhs
         error = abs(rhsvect(i) - scorp_sol(i))
         write(97,999) i,error,error/rhsvect(i)
         ave_a_error = ave_a_error + error
         ave_r_error = ave_r_error + error/rhsvect(i)
         max_a_error = max(max_a_error,error)
         max_r_error = max(max_r_error,error/rhsvect(i))
         if (error/rhsvect(i) .gt. 100) write(88,*) 'scorp:  ',
     $        i, rhsvect(i), scorp_sol(i)
      end do
      ave_a_error = ave_a_error / float(nlhs)
      ave_r_error = ave_r_error / float(nlhs)
      write(*,'(//''Average absolute scorpio error ='',
     $     f13.8)') ave_a_error
      write(*,'(''Maximum absolute scorpio error ='',
     $     f13.8)') max_a_error
      write(*,'(//''Average relative scorpio error ='',
     $     f13.8)') ave_r_error
      write(*,'(''Maximum relative scorpio error ='',
     $     f13.8)') max_r_error
      stop
      end










