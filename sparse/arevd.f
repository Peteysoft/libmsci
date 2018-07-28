      subroutine sarevd(n, nev, ncv, v, eval, maxitr, tol, which, parm)
c
c
c     Example program to illustrate the idea of reverse communication
c     for a standard nonsymmetric eigenvalue problem.
c
c     We implement example one of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit square [0,1]x[0,1] with zero Dirichlet boundary 
c         condition.
c
c     ... OP = A  and  B = I.
c
c     ... Assume "call av (nx,x,y)" computes y = A*x.c
c
c     ... Use mode 1 of SNAUPD.
c
c\BeginLib
c
c\Routines called:
c     snaupd  ARPACK reverse communication interface routine.
c     sneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     av      Matrix vector multiplication routine that computes A*x.
c     tv      Matrix vector multiplication routine that computes T*x, 
c             where T is a tridiagonal matrix.  It is used in routine
c             av.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: ndrv1.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(ncv)
      Real
     &                  ax(n), d(ncv,3), resid(n), 
     &                  v(n,ncv), workd(3*n), 
     &                  workev(3*ncv), 
     &                  workl(3*ncv*ncv+6*ncv)
      real              eval(nev, 2)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Real
     &                  tol, sigmar, sigmai
      logical           first, rvec
      integer           lun
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                  zero
      parameter         (zero = 0.0E+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real
     &                  slapy2, snrm2
      external          slapy2, snrm2, saxpy
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   | 
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
c     open log file:
      lun=10
      open (unit=lun,file="ev.log")

      bmat  = 'I'
c      which = 'SM'
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in SNAUPD as           |  
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in SNAUPD to start the Arnoldi iteration. | 
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv 
c      tol    = zero 
c      tol    = 0.00001
      ido    = 0
      info   = 0
      ldv = n
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
c      maxitr = 300
      mode   = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) | 
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine SNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call snaupd ( ido, bmat, n, which, nev, tol, resid, 
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c
            call av (n, workd(ipntr(1)), workd(ipntr(2)), parm)
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
      end if 
c 
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in SNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using SNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the imaginary part is returned   |
c        | in the second column of D.  The corresponding |
c        | eigenvectors are returned in the first NEV    |
c        | columns of the two dimensional array V if     |
c        | requested.  Otherwise, an orthogonal basis    |
c        | for the invariant subspace corresponding to   |
c        | the eigenvalues in D is returned in V.        |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SNEUPD. |
c           %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else 
c
             first  = .true.
             nconv  = iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                if (d(j,2) .eq. zero)  then
c
c                  %--------------------%
c                  | Ritz value is real |
c                  %--------------------%
c
                   call av(n, v(1,j), ax, parm)
                   call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = snrm2(n, ax, 1)
                   d(j,3) = d(j,3) / abs(d(j,1))
c
                else if (first) then
c
c                  %------------------------%
c                  | Ritz value is complex. |
c                  | Residual of one Ritz   |
c                  | value of the conjugate |
c                  | pair is computed.      | 
c                  %------------------------%
c        
                   call av(n, v(1,j), ax, parm)
                   call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = snrm2(n, ax, 1)
                   call av(n, v(1,j+1), ax, parm)
                   call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                   call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                   d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                   d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call smout(lun, nconv, 3, d, ncv, -6,
     &            'Ritz values (Real,Imag) and relative residuals')

c            we want to return the eigenvalues in addition to the
c            vectors
             do i=1, nev
               eval(i, 1)=d(i, 1)
               eval(i, 2)=d(i, 2)
             end do

          end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit
     &                  Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
c
         write (lun, *), ' '
         write (lun, *), ' _NDRV1 '
         write (lun, *), ' ====== '
         write (lun, *), ' ' 
         write (lun, *), ' Size of the matrix is ', n
         write (lun, *), ' The number of Ritz values requested is ', nev
         write (lun, *), ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         write (lun, *), ' What portion of the spectrum: ', which
         write (lun, *), ' The number of converged Ritz values is ', 
     &              nconv 
         write (lun, *), ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         write (lun, *), ' The number of OP*x is ', iparam(9)
         write (lun, *), ' The convergence criterion is ', tol
         write (lun, *), ' '
c
      end if

      close(lun)
c
c     %---------------------------%
c     | Done with program sndrv1. |
c     %---------------------------%
c
 9000 continue
c
      end

      subroutine darevd(n, nev, ncv, v, eval, maxitr, tol, which, parm)
c
c
c     Example program to illustrate the idea of reverse communication
c     for a standard nonsymmetric eigenvalue problem.
c
c     We implement example one of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit square [0,1]x[0,1] with zero Dirichlet boundary 
c         condition.
c
c     ... OP = A  and  B = I.
c
c     ... Assume "call av (nx,x,y)" computes y = A*x.c
c
c     ... Use mode 1 of SNAUPD.
c
c\BeginLib
c
c\Routines called:
c     snaupd  ARPACK reverse communication interface routine.
c     sneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     av      Matrix vector multiplication routine that computes A*x.
c     tv      Matrix vector multiplication routine that computes T*x, 
c             where T is a tridiagonal matrix.  It is used in routine
c             av.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: ndrv1.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(ncv)
      Real*8
     &                  ax(n), d(ncv,3), resid(n), 
     &                  v(n,ncv), workd(3*n), 
     &                  workev(3*ncv), 
     &                  workl(3*ncv*ncv+6*ncv)
      real*8            eval(nev, 2)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Real*8
     &                  tol, sigmar, sigmai
      logical           first, rvec
      integer           lun
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real*8
     &                  zero
      parameter         (zero = 0.0E+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real*8
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   | 
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
c     open log file:
      lun=10
      open (unit=lun,file="ev.log")

      bmat  = 'I'
      which = 'SM'
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in SNAUPD as           |  
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in SNAUPD to start the Arnoldi iteration. | 
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv 
c      tol    = zero 
      tol    = 0.00001
      ido    = 0
      info   = 0
      ldv = n
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) | 
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine SNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, 
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c
            call dav (n, workd(ipntr(1)), workd(ipntr(2)), parm)
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
      end if 
c 
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in SNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using SNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the imaginary part is returned   |
c        | in the second column of D.  The corresponding |
c        | eigenvectors are returned in the first NEV    |
c        | columns of the two dimensional array V if     |
c        | requested.  Otherwise, an orthogonal basis    |
c        | for the invariant subspace corresponding to   |
c        | the eigenvalues in D is returned in V.        |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SNEUPD. |
c           %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else 
c
             first  = .true.
             nconv  = iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                if (d(j,2) .eq. zero)  then
c
c                  %--------------------%
c                  | Ritz value is real |
c                  %--------------------%
c
                   call dav(n, v(1,j), ax, parm)
                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = snrm2(n, ax, 1)
                   d(j,3) = d(j,3) / abs(d(j,1))
c
                else if (first) then
c
c                  %------------------------%
c                  | Ritz value is complex. |
c                  | Residual of one Ritz   |
c                  | value of the conjugate |
c                  | pair is computed.      | 
c                  %------------------------%
c        
                   call dav(n, v(1,j), ax, parm)
                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = snrm2(n, ax, 1)
                   call dav(n, v(1,j+1), ax, parm)
                   call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                   call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                   d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                   d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call dmout(lun, nconv, 3, d, ncv, -6,
     &            'Ritz values (Real,Imag) and relative residuals')

c            we want to return the eigenvalues in addition to the
c            vectors
             do i=1, nev
               eval(i, 1)=d(i, 1)
               eval(i, 2)=d(i, 2)
             end do

          end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit
     &                  Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
c
         write (lun, *), ' '
         write (lun, *), ' _NDRV1 '
         write (lun, *), ' ====== '
         write (lun, *), ' ' 
         write (lun, *), ' Size of the matrix is ', n
         write (lun, *), ' The number of Ritz values requested is ', nev
         write (lun, *), ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         write (lun, *), ' What portion of the spectrum: ', which
         write (lun, *), ' The number of converged Ritz values is ', 
     &              nconv 
         write (lun, *), ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         write (lun, *), ' The number of OP*x is ', iparam(9)
         write (lun, *), ' The convergence criterion is ', tol
         write (lun, *), ' '
c
      end if
c
c     %---------------------------%
c     | Done with program sndrv1. |
c     %---------------------------%
c
 9000 continue

      close (lun)
c
      end

