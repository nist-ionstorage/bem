c
c This a driver to test fastlap V1.9 as a callable procedure from fortran
c for computing a desingularized solution of the problem of a sphere trans-
c lating in an infinite fluid.  The geometry of the problem is specified in 
c the file desing.in as is the solution.  The solution is the source strength.
c Since this depends on the desingularization distance, there is no exact 
c solution available.  We simply check the FastLap result against a result
c computed by an O(N^2) which has been around a while and we have some 
c trust in.  Output is written to desing.out.
c
c NB you cannot precondition these problems with the OL preconditioner
c as the local blocks will not be square.  You can use the SP preconditioner
c or no preconditioner.  For the test case of a sphere, the SP preconditioner
c makes matters worse, but desingularization is a funny thing in that 
c the more you do it, the more your conditioning deteriorates, so then you
c should come in and try to un-do it in some sense with a preconditioner?
c So PRECOND set to NONE in mulGlobal.h gives the best results for the 
c test input.  
c
c Written in Modern Fortran by Korsmeyer, winter 1993.  (c) MIT 1993.
c
      implicit none
      integer size
      parameter (size=2000)
      double precision x(3,4,size), poten(size), dbydnpoten(size),
     $     xcen(3,size), xcoll(3,size), xnrm(3,size), lhsvect(size), 
     $     rhsvect(size) 
      double precision d_theta
      integer type(size), shape(size), job, fljob
      integer rhstype(size), lhstype(size), rhsindex(4,size),
     $     lhsindex(4,size), dtype(size)
      integer NULL, POINT_SOURCE, CONSTANT_SOURCE, CONSTANT_DIPOLE, 
     $     LINEAR_SOURCE, LINEAR_DIPOLE
      data NULL, POINT_SOURCE, CONSTANT_SOURCE, CONSTANT_DIPOLE, 
     $     LINEAR_SOURCE, LINEAR_DIPOLE /0, 1, 11, 12, 21, 22/
      double precision exact_sol(size)
      integer nummom, numlev, i, j, k, nlhs, nrhs, nsing
      double precision error, max_error, ave_error
      data max_error, ave_error /2*0.0/
      integer fastlap
      external fastlap
      double precision tol
      integer numitr, maxitr
      data tol, maxitr /0.0001,36/
      double precision pi
      data pi /3.14159265/
c
c Query for run parameters and data file name.
c
c
      write(*,'(//''We use a desingularized single layer to solve''/
     $     ''the problem of a sphere translating in an infinite''/
     $     ''fluid. '')')
      write(*,'(//''You can only pre-condition this problem with ''/
     $     ''the SP preconditioner, but NONE gives better results. '')')
      write(*,'(''Select expansion order and tree depth:  ''$)')
      read(*,*) nummom, numlev
c
c Lines inserted to read desingularized source and collocation points.
c
      open(unit=20,file='desing.dat',status='old')
      read(20,*) nsing
      print *, ' Reading',nsing,' source and field points.' 
      nlhs = nsing
      nrhs = nsing
      do i = 1, nsing
         read(20,*) (xcoll(k,i),k=1,3),(x(k,1,i),k=1,3)
         rhsvect(i) =  -(1./2.) * xcoll(1,i)
         rhsindex(1,i) = i
         lhsindex(1,i) = i
         shape(i) = 1
         rhstype(i) = POINT_SOURCE
         lhstype(i) = POINT_SOURCE
         dtype(i) = 0
         xnrm(1,i) = 0.
         xnrm(2,i) = 0.
         xnrm(3,i) = 0.
      end do
c
c Call fastlap to compute a solution or field at the vector of field points.
c
      fljob = 2
      numitr = fastlap(nlhs,nrhs,nsing,x,shape,dtype,lhstype,rhstype,
     $     lhsindex,rhsindex,lhsvect,rhsvect,xcoll,xnrm,numlev,nummom,
     $     maxitr,tol,fljob)
      write(*,'(//i3,
     $     '' iterations knocked down residual to:'',e14.8)')
     $     numitr, tol
c
c Use the desing solution source strengths to compute the potential
c at a bunch of field points
c         
      nlhs = 100
      d_theta = pi / nlhs
      do i = 1, nlhs
         xcoll(1,i) = cos((i-1)*d_theta)
         xcoll(2,i) = 0.0
         xcoll(3,i) = sqrt(1.0-xcoll(1,i)*xcoll(1,i)) 
         exact_sol(i) =  -(1./2.)*xcoll(1,i)
      end do
      do i = 1,nrhs
         rhsvect(i) = lhsvect(i)
      end do
      fljob = 0
      numitr = fastlap(nlhs,nrhs,nsing,x,shape,dtype,lhstype,rhstype,
     $     lhsindex,rhsindex,lhsvect,rhsvect,xcoll,xnrm,numlev,nummom,
     $     maxitr,tol,fljob)
c
c Compute the error in the solution.
c
      do i = 1,nlhs
         error = abs(exact_sol(i) - lhsvect(i))
         ave_error = ave_error + error
         max_error = max(max_error,error)
      end do
      ave_error = ave_error / float(nlhs)
      write(*,'(//''Average absolute error on a meridian ='',
     $     f13.8)') ave_error
      write(*,'(''Maximum absolute error on a meridian ='',
     $     f13.8)') max_error
      stop
      end










