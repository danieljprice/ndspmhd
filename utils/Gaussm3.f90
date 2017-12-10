!     Last change:  TC   12 Sep 98    5:06 pm

!     Platform Information
!  1. Program Title: LF90
!     Command Line: lf90 <name> -nwrap -chk -vax -xref -wo -g
!     Working Directory: c:\lf90\bin

!     Files:
!     c:\lf90\bin\gaussm2.lib
!     c:\lf90\bin\gaussm2.mod
!     c:\wormyf90\pz1d\gaussm2.f90

!     VERSIONS: gaussm3
! 	This is based on gaussm2 - but with major changes that WEIGHTS & ABSCISSAS are no longer common variables
!     They are obtained from gauleg as local variables. User must specify what order of integration (# Gaussian Pts)
!     Consequently, all Gaussian Quadrature routines are modified with CALL gauleg.


!
!     9/9/98 - All ALLOCATED variables (pointers & allocatables) have been deallocated except for
!		WEIGhts and ABSCissas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTES
! 1. This MODULE can be compiled on its own like header file.
!    Then a *.lib file is created. To use functions in module, 
! 	  must have "USE module_name" statement.
! 2. FATAL -- Specification of assumed size may be done for last dimension only.
!        a(*,*) wrong      a(9,*) correct
! 3. But -- Specification of assumed shape can be done in ANY dimension 
!        a(:,:)  O.K.
! 4. "Module" is different from "include". Non-module files which have functions
!    can be included when compiling main file by using "INCLUDE" inside the main 
!    file.
! 5. Dummy arguments cannot be ALLOCATABLE. !!BUT!! Dummy arguments can have 
!    assumed shapes such as a(:,:)
! 6. INTENT must be used with dummy argument only.
! 7. Array must not have both ALLOCATABLE and POINTER attributes, but
!    pointers can be ALLOCATED.
! 8. POINTER array (eg. ans) must be declared with a deferred shape specification
!    (see "Dynamic Arrays" in the Lahey Fortran 90 Language Reference).
!    ans(3,4)  wrong          ans(:,:) correct
! 9. When dummy arguments are pointers, their actual argument must be pointers
!    too; and vice versa.
! 10. Intrinsic function MAXLOC(find location of maximum element e.g. in a list)
!    	returns an ARRAY.
! 11. Intrinsic function MATMUL for Matrix Multiplication can handle REAL of 
!     selected_kind (eg. double precision)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  1. Topics: Gaussian elimination solution of linear equations & Gaussian Qudrature!
! 	2. The variables in this module are GLOBALLY COMMON to all PROGRAM UNITS that "use"s
! this module and all INTERNAL subprograms within this module. The common variables are the
! Gauss-Legendre abscissas and weights for the Gaussian Qudrature.
! 	3. The user must first call "gauleg" from the main program to set up all the abscissas
! and weights.

module gaussm3

                        ! Common Global variables within module !
  	implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!   PRIVATE
   REAL (dbp) :: newv
   REAL(dbp)  :: EPS, M_PI
   PARAMETER (EPS=3.0d-15)       	!EPS is the relative precision
   PARAMETER (M_PI=3.141592654d0)      ! Pi value

!   PUBLIC :: newv, EPS, M_PI, n, xabsc, weig, dbp, qgss2d

   INTERFACE            

   END INTERFACE

   CONTAINS
!* This module has the following INTERNAL FUNCTIONS:
!* gauleg, qgauss, qgss3d, qgss2d, gsselm, identity_matrix
!* This module has the following INTERNAL SUBROUTINES:
!* linear_solver
!* They can call on each other without first specifying their type
!* NO INTERFACE nor EXTERNAL is required since they are INTERNAL functions

!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
	SUBROUTINE  gauleg(ngp, xabsc, weig)

      implicit none
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	   do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop

   End subroutine gauleg

!********************************************************************************
!*     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
   recursive function qgauss(func, x1, x2, ngp) RESULT(intgrl)
     	implicit none
		REAL(dbp)  intgrl, x1, x2, func
     	REAL(dbp)  xm, xl
     	INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

     	intgrl = 0.0d0
	   xm = 0.5 * (x2 + x1)
   	xl = 0.5 * (x2 - x1)
	  	do j = 1, ngp
     	   intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) )
     	END do

		intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   END function qgauss


   recursive function qgaussmat1(func, x1, x2, ngp, frow, fcol) RESULT(intgrl)
     	implicit none
     	REAL(dbp), INTENT(IN) :: x1, x2
      INTEGER :: frow, fcol
		REAL(dbp) :: intgrl(frow, fcol), tmpm(frow,fcol)
      REAL(dbp) ::  func
     	REAL(dbp) ::  xm, xl, arg
     	INTEGER j

      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

     	intgrl(:,:) = 0.0d0
      tmpm(:,:) = 0.0d0
	   xm = 0.5 * (x2 + x1)
   	xl = 0.5 * (x2 - x1)
	  	do j = 1, ngp
         arg =  xm + xl*xabsc(j)
      PRINT *, 'szhgd ds'
	        PRINT *,arg
         tmpm = func(arg)
     	   intgrl = intgrl + weig(j) * tmpm
     	END do

		intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   END function qgaussmat1



!**********************************************************************************
!*           Generic 2D Gaussian Quadraure Routines                               *
!**********************************************************************************
!* Use this function to calculate 2D integral by Gaussian Quadrature
!* Must supply boundary conditions x1,x2, y1(x), y2(x)
!* and the function to be integrated over f(x,y)
   RECURSIVE function qgss2d(origfn, xx1, xx2, yf1, yf2, ngp) RESULT(inth)
     	implicit none                               ! returns integral in inth
		REAL(dbp)  inth, xx1, xx2, yf1, yf2, origfn
     	REAL(dbp)  xm, xl, xtmp
      INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)



   interface
      real(dbp) function origfn(xp,yp) RESULT(vfun2d)     ! Original Function's Interface
                import
                implicit none
   		REAL(dbp) :: xp, yp
		end function origfn
      
      real(dbp) function yf1(x) RESULT(vy1)
                import
      	implicit none
      	REAL(dbp) :: x
      end  function yf1
      real(dbp) function yf2(x) RESULT(vy2)
                import
      	implicit none
      	REAL(dbp) :: x
      end  function yf2
	end interface

      call gauleg(ngp, xabsc, weig)

     	inth = 0.0d0
	   xm = 0.5 * (xx2 + xx1)
   	xl = 0.5 * (xx2 - xx1)

	  	do j = 1, ngp
     	   xtmp = xm + xl*xabsc(j)                ! Gauss-Legendre Abcissas
     	   inth = inth + weig(j) * qgssgy()
     	END do

		inth = inth * xl;    !Scale the answer to the range of integration  */

	CONTAINS
   	RECURSIVE function qgssgy() RESULT(intg)
     		implicit none                                ! returns integral in intg
			REAL(dbp)  intg
	     	REAL(dbp)  ym, yl, ytmp                ! all undeclared variables are
   	  	INTEGER j                                    !   COOMON with HOST

    	 	intg = 0.0d0
		   ym = 0.5 * (yf2(xtmp) + yf1(xtmp))
  		 	yl = 0.5 * (yf2(xtmp) - yf1(xtmp))

		  	do j = 1, ngp
     		   ytmp = ym + yl*xabsc(j)                ! Gauss-Legendre Abcissas
     	   	intg = intg + weig(j) * origfn(xtmp,ytmp)
	     	END do
                               PRINT *, 'Hallo'
			intg = intg * yl;    !Scale the answer to the range of integration  */
   	END function qgssgy

   END FUNCTION  qgss2d


!**********************************************************************************
!*           Generic 2D FORTRAN90 Gaussian Quadraure Routines                     * Newer than qgss2d!!!
!**********************************************************************************
!* Use this function to calculate 2D integral by Gaussian Quadrature
!* Must supply boundary conditions y1(x), y2(x)
!* and the function to be integrated over qgf902d(x,y)

!   recursive function qgss2df90(qgf902d, x1, x2) RESULT(ss)
!     	implicit none
!      REAL(dbp), INTENT(IN) :: x1, x2
! 	  	REAL(dbp) :: ss
!
!      INTERFACE
!      ! User supplied 2d function for Gaussian Quadrature Integration (FORTRAN90)
!      function qgf902d(x,y)
!         implicit none
!   		INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!         REAL(dbp), INTENT(IN) :: x
!         REAL(dbp), DIMENSION(:), INTENT(IN) :: y
!         REAL(dbp), DIMENSION(SIZE(y)) :: qgf902d
!      end function qgf902d
!      function y1d2(x)
!         implicit none
!   		INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!         REAL(dbp), INTENT(IN) :: x
!         REAL(dbp) :: y1d2
!      end function y1d2
!      function y2d2(x)
!         implicit none
!   		INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!         REAL(dbp), INTENT(IN) :: x
!         REAL(dbp) :: y2d2
!      end function y2d2
!
!
!      END INTERFACE
!
!      ss = qgausf90(h,x1,x2)                                   ! Integrate
!
!      CONTAINS
!
!	function h(x)
!   	implicit none
!      REAL(dbp), DIMENSION(:), INTENT(IN) :: x
!      REAL(dbp), DIMENSION(SIZE(x)) :: h
!      INTEGER :: ii
!      do ii = 1, SIZE(x)
!         xsav=x(ii)
!         h(ii) = qgausf90(g, y1d2(xsav), y2d2(xsav))           ! Integrate
!      end do
!   end function h
!
!   function g(y)
!   	implicit none
!      REAL(dbp), DIMENSION(:), INTENT(IN) :: y
!      REAL(dbp), DIMENSION(SIZE(y)) :: g
!      g = qgf902d(xsav, y)
!   end function g
!
!
!   RECURSIVE function qgausf90(fdum, a,b) RESULT(qgausf90r) ! Integral of f from a to b
!      implicit none
!      REAL(dbp), INTENT(IN) :: a,b
!      REAL(dbp) :: qgausf90r
!      REAL(dbp) :: xm, xl           ! x_middle,  x_lower
!
!      interface
!         function  fdum(x)
!            implicit NONE
!		   	INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!      		REAL(dbp), DIMENSION(:), INTENT(IN) :: x
!      		REAL(dbp), DIMENSION(SIZE(x)) :: fdum
!         END function  fdum
!     	end interface
!
!	   xm = 0.5d0 * (b + a)
!   	xl = 0.5d0 * (b - a)
!      qgausf90r = xl* SUM(weig(:)*fdum(xm + (xabsc(:)*xl)) )
!   end function qgausf90


!   END function qgss2df90



!**********************************************************************************
!*           Generic 3D Gaussian Quadraure Routines                               *
!**********************************************************************************
!* Use this function to calculate 3D integral by Gaussian Quadrature
!* Must supply boundary conditions x1,x2, y1(x), y2(x), z1(x,y), z2(x,y)
!* and the function to be integrated over f(x,y,z)
   RECURSIVE function qgss3d(origfn, xx1, xx2, yf1, yf2, zf1, zf2, ngp) RESULT(inth)
     	implicit none                               ! returns integral in inth
		REAL(dbp)  inth, xx1, xx2, yf1, yf2, zf1, zf2, origfn
     	REAL(dbp)  xm, xl, xtmp, ytmp
      INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)


   interface
      real(dbp) function origfn(xp,yp,zp) RESULT(vfun2d)     ! Original Function's Interface
      import
	   	implicit none
   		REAL(dbp)  :: xp, yp, zp
		end function origfn
      real(dbp) function yf1(x) RESULT(vy1)
      import
      	implicit none
      	REAL(dbp)  :: x
      end  function yf1
      real(dbp) function yf2(x) RESULT(vy2)
      import
      	implicit none
      	REAL(dbp) :: x
      end  function yf2
      real(dbp) function zf1(x,y) RESULT(vz1)
      import
      	implicit none
      	REAL(dbp)  ::  x, y
      end  function zf1
      real(dbp) function zf2(x,y) RESULT(vz2)
      import
      	implicit none
      	REAL(dbp)  :: x, y
      end  function zf2
   end interface

      call gauleg(ngp, xabsc, weig)


     	inth = 0.0d0
	   xm = 0.5 * (xx2 + xx1)
   	xl = 0.5 * (xx2 - xx1)

	  	do j = 1, ngp
     	   xtmp = xm + xl*xabsc(j)                ! Gauss-Legendre Abcissas
     	   inth = inth + weig(j) * qgssgy()
     	END do

		inth = inth * xl;    !Scale the answer to the range of integration  */

	CONTAINS
	   RECURSIVE function qgssgy() RESULT(intg)
   	  	implicit none                                ! returns integral in intg
			REAL(dbp)  intg
	     	REAL(dbp)  ym, yl		                  ! all undeclared variables are
   	  	INTEGER j                                    !   COOMON with HOST

	     	intg = 0.0d0
		   ym = 0.5 * (yf2(xtmp) + yf1(xtmp))
   		yl = 0.5 * (yf2(xtmp) - yf1(xtmp))

		  	do j = 1, ngp
   	  	   ytmp = ym + yl*xabsc(j)                ! Gauss-Legendre Abcissas
     		   intg = intg + weig(j) * qgssfz()
	     	END do

			intg = intg * yl;    !Scale the answer to the range of integration  */
	   END function qgssgy

  	 	RECURSIVE function qgssfz() RESULT(intf)
     		implicit none                                ! returns integral in intg
			REAL(dbp)  intf
   	  	REAL(dbp)  zm, zl, ztmp                ! all undeclared variables are
     		INTEGER j                                    !   COOMON with HOST

	     	intf = 0.0d0
		   zm = 0.5 * (zf2(xtmp,ytmp) + zf1(xtmp,ytmp))
   		zl = 0.5 * (zf2(xtmp,ytmp) - zf1(xtmp,ytmp))

		  	do j = 1, ngp
   	  	   ztmp = zm + zl*xabsc(j)                ! Gauss-Legendre Abcissas
     		   intf = intf + weig(j) * origfn(xtmp,ytmp,ztmp)
	     	END do

			intf = intf * zl;    !Scale the answer to the range of integration  */
	   END function qgssfz


   END FUNCTION  qgss3d



!*  ************************************************************ */
!*    Solving a system of linear equations by Gauss elimination
!*      1st  matrix element -> (1,1)
!*      Last matrix element -> (row,row)
!*      plus an augmented column in the (row+1)th column
!*      Input: a matrix of coefficients plus an augmented column - a(row,row+1)
!*      Output: an array of solutions - Gsselm(row)                    */
!*  ************************************************************ */
function Gsselm(a,row)
	implicit none
	INTEGER, INTENT(IN) :: row
	REAL(dbp) , INTENT(IN OUT)  ::  a(:,:)   	!Assume shape (:)
	REAL(dbp) , DIMENSION(row) :: Gsselm

	INTEGER i,j,k
	INTEGER, DIMENSION(2) :: shap
	REAL(dbp) , ALLOCATABLE :: swap_ik(:)
	REAL(dbp)  :: tmp

	ALLOCATE (swap_ik(row+1))

! 	Initialise
	swap_ik(:) = 0.0d0            ! Whole vector initialized to zero
	tmp = 0.0d0

! Check dimensions of input matrix
	shap = SHAPE(a)
	if ( (shap(1) .NE. row) .OR.  (shap(2) .NE. row+1) ) then
   	stop 'something wrong in shape'
	end if


!/*   Gaussian Elimination - Row Reduction of matrix  */
	do k=1, row-1                             ! total of row-1 operations

!/*  Pivotal strategy - SWAP rows to make pivotal element a[k][k] have the
!    greatest magnitude in its column. This prevents unnecessary division by
!    a small number.                           */
   	do i = k+1, row
      	if ( (dabs(a(i,k))-dabs(a(k,k))).gt.eps  ) then
         	do j = k, row+1                     !/* If pivotal element is not */
            	swap_ik(j) = a(k,j)              !/* the highest then  */
          	   a(k,j) = a(i,j)                  !/* swap i'th and k'th rows */
            	a(i,j) = swap_ik(j)
         	end do 		!j-loop
         end if
   	end do 				!i-loop


!/*   If the Matrix is SINGULAR then EXIT program      */
  		IF ( dabs(a(k,k)) < EPS ) then
   		print *,'After swapping rows to make the pivotal element become the'
      	print *,'highest magnitude element in its column, its magnitude is'
	      print *,'still extremely small.'
	      print *,'Hence this is a SINGULAR MATRIX - no unique solution or '
	      print *,'check that the input dimensions are correct.'
!   	   call break()
	   END if
   


!/*      Perform row-reduction with pivotal element a[k][k]     */
		do i = k+1, row
			do j = row+1, k, -1			!/* starting from end of column */
	     	 	a(i,j) = a(i,j) - a(k,j) / a(k,k) * a(i,k)
			end DO 							!/* end of j loop     */
		end do 	 							!/* end of 2nd i loop */

	end DO 									!/* end of k loop     */
!  At this point, the bottom triangle is Zero


!/*   Back Substitution - Solutions of equations   */
	Gsselm(row) = a(row,row+1) / a(row,row)
	do k = row-1, 1, -1
   	tmp = 0.0d0
   	do j = k+1, row
      	tmp = tmp + a(k,j)*Gsselm(j)
	   end do 							!j-loop
   	Gsselm(k) = ( a(k,row+1) - tmp ) / a(k,k)
	end do 								!k-loop

   deallocate (swap_ik)
	RETURN
END function Gsselm



!*  ************************************************************ */
!  A function that produces an Identity Matrix of dimension (n,n) !
!*  ************************************************************ */
function identity_matrix(nn)
   implicit none
   INTEGER, INTENT(IN) :: nn
   REAL(dbp) , DIMENSION(nn,nn) :: identity_matrix
									! This function returns an n X n array.
   INTEGER :: j            ! local loop index
   
   identity_matrix = 0.0   ! Set each element to zero

   DO j = 1,nn                      ! Change value of each element
      identity_matrix(j,j) = 1.0    ! along diagonal to zero
   END DO
end function identity_matrix




!*  ************************************************************
!*    Solving MULTIPLE systems of linear equations by Gauss-Jordan Elimination
!*    Row operations are performed until Reduced-Row Echelon form is obtained;
!*    i.e. original matrix of coefficients becomes unit matrix.

!* USAGE: Multiple Linear Solver / Inverse Matrix / Determinant
!*                To find the inverse matrix of aaa
!* call Linear_Solver(aaa, DBLE(identity_matrix(ii)), ii, ii, xxx, determ)
!*                TO find determinant ONLY without solving equations
!* -> leave "b" as undefined and set "soln = 0"

!* 1. INPUT: a(row,row) -> Matrix of coefficients of system.
!* 2. INPUT: b(row,soln) -> Vectors(matrix) of solutions.
!*    row = # of equations,      soln = # of solution vectors
!* 3. LOCAL: c(row,row+soln) -> Augmented matrix.
!* 4. OUTPUT:ans(row,soln) -> #soln unknown(answer) vectors
!* 5. DBLE(...) -> to make into double precision
!* 6. Must open file UNIT=60 from main calling routine to write to file
!*  ************************************************************
subroutine Linear_Solver(a,b,row,soln,ans,detmt)
	implicit none
	INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)
	REAL(dbp) , INTENT(OUT) :: detmt

	INTEGER i,j, k(1),p          	! need to specify p with dimension one to
                                 ! use MAXLOC intrinsic function
	REAL(dbp) , POINTER  :: swap_ij(:), c(:,:)
	REAL(dbp)  :: tmp

!	CHARACTER :: pauss
!   REAL(dbp)  :: iin

	ALLOCATE (c(row,row+soln), swap_ij(row+soln))
	detmt = 1.0d0

! Constructing ONE single augmented matrix (actually c is a pointer)
  	c(1:row, 1:row) = a(1:row, :)
   c(1:row, row+1 : row+soln) = b(1:row, :)

!PRINT *,'MAKING AUG whole'
!do i = 1,row
!    PRINT "(8F9.5)",c(i,:)
!end do
!PRINT *,'end'


! Begin Row-Reduction Operations
	do j = 1, row                 ! for the j-th column

!PRINT *,'bef col'
!do iin = 1, row
!	PRINT "(8F9.5)",c(iin,:)
!end do
!PRINT *,'end'
!READ *,pauss

      !WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
		!do p = 1,row
		!	WRITE(60,500) (c(p,i), i = 1,row)
      !end do
		!WRITE(60,510)

   	k = j - 1 + MAXLOC( dabs(c(j:row , j)) )        ! assigning to the k-array
	   p = k(1)


! Find row p that has the largest element in column j
! among the elements below the diagonal element.

!/*  Pivotal strategy - SWAP rows to make pivotal element c[j][j] have the
!    greatest magnitude in its column. This prevents unnecessary division by
!    a small number.                           */
   	if (p .ne. j) then
   		swap_ij(:) = c(j,:);    c(j,:) = c(p,:);   c(p,:) = swap_ij(:)
	      detmt = -1.0d0 * detmt
   	end if
! Determinant change signs if rows are swapped.


! If after swapping rows the diagonal element(now having the largest value)
! is still very small then Matrix is singular
	   if ( (dabs(c(j,j)) < 1.0d-10*eps).AND.(j.ne.row)) then
   	   PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
         PRINT *, 'diagonal element is ', c(j,j)
      	call exit(10)
	  	end if

!PRINT *,'swap col'
!do iin = 1, row
!    PRINT "(8F9.5)",c(iin,:)
!end do
!PRINT *,'end'
!READ *,pauss


! Divide the j-th row by leading diagonal
   	tmp = c(j,j);     c(j,:) = c(j,:) / tmp
   
! Finding Determinant as the factors of the diagoanals
	   detmt = detmt * tmp
      WRITE(60,*) j, tmp,detmt

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
   	do i = 1, row
      	if (i .ne. j) then
      		tmp = c(i,j)
         	c(i,:) = c(i,:) - tmp * c(j,:)
      	end if
	   end do

!PRINT *,'zer col'
!do iin = 1, row
!    PRINT "(8F9.5)",c(iin,:)
!END do
!READ *,pauss


	end do      ! j loop

!PRINT *,'whole'
!do j = 1, row
!	PRINT "(8F9.5)",c(j,:)
!end do

   ans = c(1:row , row+1:row+soln)
	deallocate (c, swap_ij)

!   WRITE(60,*)'Linear Solver Subroutine: Original matrix'
!	do j = 1, row
!		WRITE(60,500) (a(j,i), i = 1,row)
!	end do
!   WRITE(60,510)

!   WRITE(60,*)'Linear Solver Subroutine: Solution matrix'
!	do j = 1, row
!		WRITE(60,500) (b(j,i), i = 1,soln)
!	end do
!   WRITE(60,510)

!   WRITE(60,*)'Result from Linear Solver Subroutine'
!	do j = 1, row
!		WRITE(60,500) (ans(j,i), i = 1,soln)
!	end do
!
!   WRITE(60,510)
!   WRITE(60,*)'determinant is ',detmt


500  FORMAT(300g15.4)
510  FORMAT(//)

END subroutine Linear_Solver


!**************************************************************************
! Features of LinSolv2:
! LinSolv2 like Linear_Solver but uses row & column pivoting
! Theorems of DETERMINANTS (from Anton & Rorres, 1987).
!  1. If A is any square matrix that contains a row of zeros, then det(A) = 0
!  2. If A is a nXn triangular matrix(bottom triangles are zeros) then det(A) is the
!     product of the diagonal entries.
!  3. Theorem regarding row operations:
!  a) If a row of original matrix A is divided by constant k and the resultant matrix is A'
!	 	then det(A) = k det(A')
!  b) If A' is the matrix that results from interchanging two rows of matrix A 
!		then det(A) = - det(A')
!  c) Adding/Subtracting a multiple of one row to another row does not affect the determinant.
!  4. If A is any square matrix then det(A) = det(A').
!  5. If A and B are square matrices of the same size then det(AB) = det(A)det(B)
!  6. A square matrix A is invertible if and only if det(A)<>0. i.e. det(A^-1) = 1/(det(A)).
!  7. If A is an invertible matrix then A^-1 = adj(A)/det(A).
!  Postulate:
!  1. Changing columns (i.e. changing the order of variables) does not affect the determinant.
!  Pivoting:- rows and columns pivoted independently
!     1. For every row of the submatrix, find the largest element along the columns and
!        divide that particular row with it - Normalize each row.
!     2. For the j-th row at the j-th operation, swap columns to have the largest element of 
!        the j-th row at the (1,1) position of the submatrix.
!     3. Look for the largest element in the j-th column, during j-th operation and swap rows.
!  Note this is different from IMPLICIT pivoting from step 2 onwards where it is decided based
!  on the largest element along column j, whether to swap rows or not. When it is decided 
! 	which row to swap, then all operations are performed on the original matrix before step 1 
!  is applied.

subroutine LinSolv2(a,b,row,soln,ans,detmt)
	implicit none
	INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)
	REAL(dbp) , INTENT(OUT) :: detmt

	INTEGER i,j, k(1),p          	! need to specify p with dimension one to
                                 ! use MAXLOC intrinsic function
   INTEGER :: imm, memoryvec(row), iw,ix,iy,iz, itmp
	REAL(dbp) , POINTER  :: swap_ij(:), c(:,:)
	REAL(dbp)  :: tmp
   REAL(dbp), ALLOCATABLE :: ctmp(:,:)

	ALLOCATE (c(row,row+soln), swap_ij(row+soln))
   ALLOCATE (ctmp(row,row+soln))
	detmt = 1.0d0

! ERROR Check #1
! Checking for rows of zeros (Theorem 1)
   do i = 1, row
      tmp = 0.0d0
      do j = 1, row 
         tmp = tmp + dabs( a(i,j) )
      end do
      if (tmp.lt.eps*eps) then
         PRINT *, 'Error the sum of row ',i,' is less than ',eps*eps
         PRINT *, 'High possibility that matrix is singular!!!'
         call EXIT(10)
      end if
   end do

! Constructing ONE single augmented matrix (actually c is a pointer)
	do i = 1, row
   	c(i, 1:row) = a(i,:)
	   c(i, row+1 : row+soln) = b(i,:)
	end do

! Initializing memory vector for the index/position for solution 
! Must be used if there is column swapping
   do i = 1, row
      memoryvec(i) = i
   end do


! Begin Row-Reduction Operations
	do j = 1, row                 ! for the j-th column

         WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
		do p = 1,row
			WRITE(60,500) (c(p,i), i = 1,row+soln)
      end do
		WRITE(60,510)

!  PIVOTING - of the submatrix for the j-th operation
!  i.e. for every row, normalize the entire row by the element with the largest magnitude
!  ONLY do so provided that largest value in a row is < 1.0
  		do i = j, row
      	tmp = MAXval( dabs(c(i,j:row)) )
      	if ( (tmp.lt.1.0d0) .AND. (tmp.gt.eps) ) then
            do p = j, row+soln  			! the preceeding columns in this row are ZEROS anyway!
               c(i,p) = c(i,p) / tmp
            end do
      		detmt = detmt * tmp
	      end if
   	end do         ! i loop


            ! Finding location of max value along row j
		k = j - 1 + MAXLOC( dabs(c(j , j:row)) )        ! assigning to the k-array
      imm = k(1)                                     	! Swap row/col  imm, j

      if (imm.ne.j) then
	      do iw = 1, row
 				if (iw.eq.imm) then
    				ix=j
       		elseif (iw.eq.j) then
			     	ix=imm
   	    	else
			     	ix=iw
       		ENDIF
				do iy = 1, row
   	 			if (iy.eq.imm) then
      	 			iz=j
     				elseif (iy.eq.j) then
       				iz=imm
	    			else
   	    			iz=iy
     				ENDIF
    				ctmp(ix,iz) = c(iw,iy)     ! Building up temporary matrix to represent
	 			end do     ! end loop iy      ! the matrix where the row/col is changed imm<->i
            do iy = row+1 , row+soln
               ctmp(ix,iy) = c(iw,iy)
            end do
			end do    	  ! end for iw
      	itmp = memoryvec(j)             	! Keep track of swapped variables/solutions
	      memoryvec(j) = memoryvec(imm)
   	   memoryvec(imm) = itmp
! Actual swapping of row/col imm<->j
	      c(:,:) = ctmp(:,:)
      end if


! Find row p that has the largest element in column j
! among the elements below the diagonal element.

!  Pivotal strategy - SWAP rows to make pivotal element c[j][j] have the
!  greatest magnitude in its column. This prevents unnecessary division by
!  a small number.
   	k = j - 1 + MAXLOC( dabs(c(j:row , j)) )        ! assigning to the k-array
	   p = k(1)

   	if (p .ne. j) then
   		swap_ij(:) = c(j,:);    c(j,:) = c(p,:);   c(p,:) = swap_ij(:)
	      detmt = -1.0d0 * detmt
   	end if
! Determinant change signs if rows are swapped (Theorem 3b.)



! ERROR Check #2
! If after swapping rows the diagonal element(now having the largest value)
! is still very small then Matrix is singular
	   if ( (dabs(c(j,j)) < 1.0d-10*eps).AND.(j.ne.row)) then
   	   PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
         PRINT *, 'diagonal element is ', c(j,j)
      	call exit(10)
	  	end if

! ERROR Check #3
! If at the j-th row, all other elements are zero and element(j,j) is very small
! then might be singular or inconsistent matrix
      tmp = 0.0d0
      do i = 1, j-1
         tmp = tmp + dabs(c(j,i))
      end do
      do i = j+1, row
         tmp = tmp + dabs(c(j,i))
      end do                 
      if ( (tmp.lt.eps) .AND. (dabs(c(j,j)).lt.100*eps) ) then
   	   PRINT *, 'ERROR: Matrix is Singular/Inconsistent. Found at ',j,'-th row operation'
         PRINT *, 'Diagonal element is too small ', c(j,j)
         PRINT *, 'And all other elements in this row are almost zero'
      	call exit(10)
      end if

! Divide the j-th row by leading diagonal
   	tmp = c(j,j)
      c(j,j) = 1.0d0
      DO i = j+1, row+soln
        	c(j,i) = c(j,i) / tmp
      END do  



! Finding Determinant as the factors of the diagoanals (Theorem 3a.)
	   detmt = detmt * tmp

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
   	do i = 1, row                                      ! 1 0 ~ ~
      	if (i .ne. j) then                              ! 0 1 ~ ~
      		tmp = c(i,j)                                 ! 0 0 1 ~
            c(i,j) = 0.0d0                               ! 0 0 ~ ~
            write (60,*) i, tmp
            do p = j+1, row+soln
               c(i,p) = c(i,p) - tmp * c(j,p)
            end do
      	end if
	   end do


	end do      ! j loop

   ans(memoryvec(:),:) = c(:,row+1:row+soln)

	deallocate (c, ctmp, swap_ij)

   WRITE(60,*)'Linear Solver Subroutine: Original matrix'
	do j = 1, row
		WRITE(60,500) (a(j,i), i = 1,row)
	end do
   WRITE(60,510)

   WRITE(60,*)'Linear Solver Subroutine: Solution matrix'
	do j = 1, row
		WRITE(60,500) (b(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'Result from Linear Solver Subroutine'
	do j = 1, row
		WRITE(60,500) (ans(j,i), i = 1,soln)
	end do
   WRITE(60,510)
   WRITE(60,*)'determinant is ',detmt

500  FORMAT(300g15.4)
510  FORMAT(//)
520  FORMAT(30i3)

END subroutine LinSolv2



!**************************************************************************
! Features of LinSolv3
! LinSolv3 like LinSolv2 but uses IMPLICIT pivoting (See Numerical Recipes)
!  Postulate:
!  1. Changing columns (i.e. changing the order of variables) does not affect the determinant.
!  Pivoting:- rows and columns pivoted independently
!     1. For every row of the submatrix, find the largest element along the columns and
!        divide that particular row with it - Normalize each row. Actually, all that is needed
!        in this step is to collect the elements in the j-th column under "pivcol(:)"
!  	2.	If the largest element in "pivcol" is NOT in row j then swap rows with row j.

subroutine LinSolv3(a,b,row,soln,ans,detmt)
	implicit none
	INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)
	REAL(dbp) , INTENT(OUT) :: detmt

	INTEGER i,j, k(1),p          	! need to specify p with dimension one to
                                 ! use MAXLOC intrinsic function
	REAL(dbp) , POINTER  :: swap_ij(:), c(:,:)
	REAL(dbp)  :: tmp
   REAL(dbp), ALLOCATABLE ::  pivcol(:)


	ALLOCATE (c(row,row+soln), swap_ij(row+soln))
	detmt = 1.0d0

! ERROR Check #1
! Checking for rows of zeros (Theorem 1)
   do i = 1, row
      tmp = 0.0d0
      do j = 1, row 
         tmp = tmp + dabs( a(i,j) )
      end do
      if (tmp.lt.eps*eps) then
         PRINT *, 'Error the sum of row ',i,' is less than ',eps*eps
         PRINT *, 'High possibility that matrix is singular!!!'
         call EXIT(10)
      end if
   end do

! Constructing ONE single augmented matrix (actually c is a pointer)
	do i = 1, row
   	c(i, 1:row) = a(i,:)
	   c(i, row+1 : row+soln) = b(i,:)
	end do

   ALLOCATE ( pivcol(row) )
   pivcol(:) = 0.0d0


! Begin Row-Reduction Operations
	do j = 1, row                 ! for the j-th column


      WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
		do p = 1,row
			WRITE(60,500) (c(p,i), i = 1,row+soln)
      end do
		WRITE(60,510)

!  IMPLICIT PIVOTING (step1)- of the submatrix for the j-th operation
!  i.e. for every row, normalize the entire row by the element with the largest magnitude
!  Actually just want to know about the j-th column after normalization
  		do i = j, row
      	tmp = MAXval( dabs(c(i,j:row)) )
         pivcol(i) = c(i,j) / tmp
   	end do         ! i loop


! (step2) Find row p that has the largest element in column j
! among the elements below the diagonal element
! as if the rows are normalized.
!  Pivotal strategy - SWAP rows to make pivotal element c[j][j] have the
!  greatest magnitude in its column. This prevents unnecessary division by
!  a small number.
   	k = j - 1 + MAXLOC( dabs(pivcol(j:row)) )        ! assigning to the k-array
	   p = k(1)

   	if (p .ne. j) then
   		swap_ij(:) = c(j,:);    c(j,:) = c(p,:);   c(p,:) = swap_ij(:)
	      detmt = -1.0d0 * detmt
   	end if
! Determinant change signs if rows are swapped (Theorem 3b.)



! ERROR Check #2
! If after swapping rows the diagonal element(now having the largest value)
! is still very small then Matrix is singular
	   if ( (dabs(c(j,j)) < 1.0d-10*eps).AND.(j.ne.row)) then
   	   PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
         PRINT *, 'diagonal element is ', c(j,j)
      	call exit(10)
	  	end if

! ERROR Check #3
! If at the j-th row, all other elements are zero and element(j,j) is very small
! then might be singular or inconsistent matrix
      tmp = 0.0d0
      do i = 1, j-1
         tmp = tmp + dabs(c(j,i))
      end do
      do i = j+1, row
         tmp = tmp + dabs(c(j,i))
      end do                 
      if ( (tmp.lt.eps) .AND. (dabs(c(j,j)).lt.100*eps) ) then
   	   PRINT *, 'ERROR: Matrix is Singular/Inconsistent. Found at ',j,'-th row operation'
         PRINT *, 'Diagonal element is too small ', c(j,j)
         PRINT *, 'And all other elements in this row are almost zero'
      	call exit(10)
      end if

! Divide the j-th row by leading diagonal
   	tmp = c(j,j)
      c(j,j) = 1.0d0
      DO i = j+1, row+soln
        	c(j,i) = c(j,i) / tmp
      END do  


! Finding Determinant as the factors of the diagoanals (Theorem 3a.)
	   detmt = detmt * tmp

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
   	do i = 1, row                                      ! 1 0 ~ ~
      	if (i .ne. j) then                              ! 0 1 ~ ~
      		tmp = c(i,j)                                 ! 0 0 1 ~
            c(i,j) = 0.0d0                               ! 0 0 ~ ~
				do p = j+1, row+soln
               c(i,p) = c(i,p) - tmp * c(j,p)
            end do
      	end if
	   end do

      WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
		do p = 1,row
			WRITE(60,500) (c(p,i), i = 1,row+soln)
      end do
		WRITE(60,510)


	end do      ! j loop

   ans = c(1:row , row+1:row+soln)
   deallocate (c, pivcol, swap_ij)

   WRITE(60,*)'Linear Solver Subroutine: Original matrix'
	do j = 1, row
		WRITE(60,500) (a(j,i), i = 1,row)
	end do
   WRITE(60,510)

   WRITE(60,*)'Linear Solver Subroutine: Solution matrix'
	do j = 1, row
		WRITE(60,500) (b(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'Result from Linear Solver Subroutine'
	do j = 1, row
		WRITE(60,500) (ans(j,i), i = 1,soln)
	end do
   WRITE(60,510)
   WRITE(60,*)'determinant is ',detmt

500  FORMAT(300g15.4)
510  FORMAT(//)
520  FORMAT(30i3)

END subroutine LinSolv3


!**************************************************************************
! Features of LinSolv4:
! LinSolv4 like Linear_Solver but uses FULL (simultaneous row & column) pivoting
!  Postulate:
!  1. Changing columns (i.e. changing the order of variables) does not affect the determinant.
!  Pivoting:-
!     1. For the submatrix led by the element at position (j,j), find the largest element 
!        as the pivot. Record its position in prow, pcol
!     2. Swap rows and cols to make the pivot at position (j,j)

subroutine LinSolv4(a,b,row,soln,ans,detmt)
	implicit none
	INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)
	REAL(dbp) , INTENT(OUT) :: detmt

	INTEGER i,j, p          	! need to specify p with dimension one to
                                 ! use MAXLOC intrinsic function
   INTEGER :: memoryvec(row), iw,ix,iy,iz, itmp, prow, pcol
	REAL(dbp) , POINTER  :: swap_ij(:), c(:,:)
	REAL(dbp)  :: tmp
   REAL(dbp), ALLOCATABLE :: ctmp(:,:)

	ALLOCATE (c(row,row+soln), swap_ij(row+soln))
   ALLOCATE (ctmp(row,row+soln))
	detmt = 1.0d0

! ERROR Check #1
! Checking for rows of zeros (Theorem 1)
   do i = 1, row
      tmp = 0.0d0
      do j = 1, row 
         tmp = tmp + dabs( a(i,j) )
      end do
      if (tmp.lt.eps*eps) then
         PRINT *, 'Error the sum of row ',i,' is less than ',eps*eps
         PRINT *, 'High possibility that matrix is singular!!!'
         call EXIT(10)
      end if
   end do


! Constructing ONE single augmented matrix (actually c is a pointer)
	do i = 1, row
   	c(i, 1:row) = a(i,:)
	   c(i, row+1 : row+soln) = b(i,:)
	end do



   do i = 1, row
   	do j = 1, row + soln
 !  		if (c(i,j) .lt. 1.0d-2*eps) then
 !           c(i,j) = 0.0d0
 !	     	end if
      end do
   end do



! Initializing memory vector for the index/position for solution 
! Must be used if there is column swapping
   do i = 1, row
      memoryvec(i) = i
   end do


! Begin Row-Reduction Operations
	do j = 1, row                 ! for the j-th column

         WRITE(60,*)'Current Matrix: Original matrix ',j,'-th operation'
		do p = 1,row
			WRITE(60,500) (c(p,i), i = 1,row+soln)
      end do
		WRITE(60,510)

!  PIVOTING - of the submatrix for the j-th operation
!  Search for Pivot from whole submatrix led by (j,j)
      tmp = 0.0d0;      prow = j; 		pcol = j
  		do i = j, row
         do p = j, row
            if ( tmp .lt. dabs(c(i,p)) ) then         ! tmp < c(i,p)
               tmp = dabs(c(i,p))
               prow = i;   pcol = p
				end if
        	end do
   	end do         ! i loop

! Swap columns NB swapping column must be accompanied by swapping row.
! Swap row/col  pcol <->  j
      if (pcol.ne.j) then
         WRITE(60,*) 'col change from - to ',pcol,j
	      do iw = 1, row
 				if (iw.eq.pcol) then
    				ix=j
       		elseif (iw.eq.j) then
			     	ix=pcol
   	    	else
			     	ix=iw
       		ENDIF
				do iy = 1, row
   	 			if (iy.eq.pcol) then
      	 			iz=j
     				elseif (iy.eq.j) then
       				iz=pcol
	    			else
   	    			iz=iy
     				ENDIF
    				ctmp(ix,iz) = c(iw,iy)     ! Building up temporary matrix to represent
	 			end do     ! end loop iy      ! the matrix where the row/col is changed imm<->i
            do iy = row+1 , row+soln
               ctmp(ix,iy) = c(iw,iy)
            end do
			end do    	  ! end for iw
      	itmp = memoryvec(j)             	! Keep track of swapped variables/solutions
	      memoryvec(j) = memoryvec(pcol)
   	   memoryvec(pcol) = itmp
! Actual swapping of row/col imm<->j
	      c(:,:) = ctmp(:,:)
! Row position of pivot might also change after changing columns/rows
         if (prow.eq.j) then
            prow = pcol
         ELSEIF (prow.eq.pcol) then
           	prow = j
        	end if
      end if               ! end column swap
                           ! ASSUMING swapping col/row do NOT change determinant



! Swap rows
   	if (prow .ne. j) then
   		swap_ij(:) = c(j,:);    c(j,:) = c(prow,:);   c(prow,:) = swap_ij(:)
	      detmt = -1.0d0 * detmt
   	end if
! Determinant change signs if rows are swapped (Theorem 3b.)


! ERROR Check #2
! If after PIVOTING the diagonal element(now having the largest value)
! is still very small then Matrix is singular
! provided c(j,j) is not the last diagonal element
	   if ( (dabs(c(j,j)) < 1.0d-12*eps).AND.(j.ne.row)) then
   	   PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
         PRINT *, 'diagonal element is ', c(j,j)
      	call exit(10)
	  	end if

! ERROR Check #3
! If at the j-th row, all other elements are zero and element(j,j) is very small
! then might be singular or inconsistent matrix
      tmp = 0.0d0
      do i = 1, j-1
         tmp = tmp + dabs(c(j,i))
      end do
      do i = j+1, row
         tmp = tmp + dabs(c(j,i))
      end do                 
      if ( (tmp.lt.eps) .AND. (dabs(c(j,j)).lt.1.0d-12*eps) ) then
   	   PRINT *, 'ERROR: Matrix is Singular/Inconsistent. Found at ',j,'-th row operation'
         PRINT *, 'Diagonal element is too small ', c(j,j)
         PRINT *, 'And all other elements in this row are almost zero'
      	call exit(10)
      end if

! Divide the j-th row by leading diagonal
   	tmp = c(j,j)
      c(j,j) = 1.0d0
      DO i = j+1, row+soln
        	c(j,i) = c(j,i) / tmp
      END do  


! Finding Determinant as the factors of the diagoanals (Theorem 3a.)
	   detmt = detmt * tmp

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
   	do i = 1, row                                      ! 1 0 ~ ~
      	if (i .ne. j) then                              ! 0 1 ~ ~
      		tmp = c(i,j)                                 ! 0 0 1 ~
            c(i,j) = 0.0d0                               ! 0 0 ~ ~
            write (60,*) i, tmp
            do p = j+1, row+soln
               c(i,p) = c(i,p) - tmp * c(j,p)
!               if (dabs(c(i,p)).lt.1.0d-16) c(i,p) = 0.0d0
            end do
      	end if
	   end do


	end do      ! j loop

   ans(memoryvec(:),:) = c(:,row+1:row+soln)

	deallocate (c, ctmp, swap_ij)

   WRITE(60,*)'Linear Solver Subroutine: Original matrix'
	do j = 1, row
		WRITE(60,500) (a(j,i), i = 1,row)
	end do
   WRITE(60,510)

   WRITE(60,*)'Linear Solver Subroutine: Solution matrix'
	do j = 1, row
		WRITE(60,500) (b(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'Result from Linear Solver Subroutine'
	do j = 1, row
		WRITE(60,500) (ans(j,i), i = 1,soln)
	end do
   WRITE(60,510)
   WRITE(60,*)'determinant is ',detmt

500  FORMAT(300g15.4)
510  FORMAT(//)
520  FORMAT(30i3)

END subroutine LinSolv4




!****************************************************************************
!***  Subroutine that rearrange the rows and columns of matrix kmat      ***!
!****************************************************************************
! Consider a a system of equation Ax=b
! "kmat" is the augmented matrix consisting of [A|b] i.e. matrix A and coumn b
! Example:  x1 x2 x3 x4 x5 x6 x7 x8
! know about x3, x8, x6, x7 --- put these as the first four elements
! In this example: "totpnt" = 8, "known" =4, "posn"({1,2,3,4}) = {3,8,6,7}
! Generate binary sequence {00100111}   aim to get {11110000}
! "xvec" is the original un arranged x-vector --- {1,2,3,4,5,6,7,8}
! After arrangement, xvec would be different and this info is outputted.
! "augcol" is the number of augmeted columns --- eg. augcol=1 if b has one column
subroutine arrange(kmat, posn, xvec, totpnt, known, augcol)
   implicit none
!Passed Variables
   INTEGER, INTENT(IN) :: totpnt, known, augcol
   INTEGER, INTENT(IN) :: posn(known)
   INTEGER, INTENT(OUT) :: xvec(totpnt)
   REAL(dbp) :: kmat(totpnt,totpnt+augcol)

!Local Variables
   INTEGER :: ii,ij, xvectmp
	INTEGER, ALLOCATABLE :: binseq(:)
!	REAL(dbp), ALLOCATABLE :: atmp(:,:)
   REAL(dbp), ALLOCATABLE :: atmprowi(:), atmprowj(:), atmpcoli(:), atmpcolj(:)
   REAL(dbp) :: aii, aij, aji, ajj

! Initial Data
!	totpnt = 8
!	known = 4

!   ALLOCATE (atmp(totpnt,totpnt+augcol))
   ALLOCATE ( atmprowi(totpnt), atmprowj(totpnt), atmpcoli(totpnt+augcol), atmpcolj(totpnt+augcol) )
   ALLOCATE (binseq(totpnt) )

!	posn(1) = 3
!	posn(2) = 8
!	posn(3) = 6
!	posn(4) = 7


  	binseq(:) = 0              	!/* initialize binary sequence  */
	do ii =1, totpnt
     	xvec(ii) = ii               !/* initialize original sequence of solution vector  */
	ENDDO

!/* creating a binary sequence representing knowns and unknowns  */
! 1 = known    &  0 = unknown
! Eg. Originally {00100111}   -> require  {11110000}
     	binseq(posn(1:known)) = 1

	do ii = 1, totpnt
  		if (binseq(ii) .eq. 0) then
         do ij = totpnt, ii+1, -1
           	if (binseq(ij) .eq. 1) then
!/* 1. swap row/col  ii<>ij  */
               PRINT *, 'swap i <-> j ',ii,ij
!				   atmp(:,:) = kmat(:,:)         ! assign the whole matrix
!				   atmp(ii,:) = kmat(ij,:)       ! swap rows  ii<->ij
!				   atmp(ij,:) = kmat(ii,:)
!				   atmp(:,ii) = kmat(:,ij)       ! swap cols  ii<->ij
!				   atmp(:,ij) = kmat(:,ii)
! 					atmp(ii, ii) = kmat(ij, ij)   ! swap the 4 elements
!				   atmp(ii, ij) = kmat(ij, ii)
!				   atmp(ij, ii) = kmat(ii, ij)
!				  	atmp(ij, ij) = kmat(ii, ii)
               atmpcoli(:) = kmat(ii,:)       ! swap rows  ii<->ij
               atmpcolj(:) = kmat(ij,:)       ! swap rows  ii<->ij
               atmprowi(:) = kmat(:,ii)       ! swap rows  ii<->ij
               atmprowj(:) = kmat(:,ij)       ! swap rows  ii<->ij
 					aii = kmat(ii, ii)   ! swap the 4 elements
				   aij = kmat(ii, ij)
				   aji = kmat(ij, ii)
				  	ajj = kmat(ij, ij)
               kmat(ii,:)   =    atmpcolj(:)
               kmat(ij,:)   =    atmpcoli(:)
               kmat(:,ii)   =    atmprowj(:)
               kmat(:,ij)   =    atmprowi(:)
               kmat(ii, ii) = 	ajj
               kmat(ii, ij) =    aji
               kmat(ij, ii) =    aij
               kmat(ij, ij) =   	aii
!/*  2. After swapping matrix -> Change Sequence   */
             	binseq(ii) = 1
               binseq(ij) = 0
!/*  3. Assign the swapped matrix (atmp) back to its original name (kmat) */
!					kmat(1:totpnt, 1:totpnt+augcol) = atmp(1:totpnt, 1:totpnt+augcol)
!/* 4. Reorganize solution vector  */
              	xvectmp = xvec(ii)
              	xvec(ii) = xvec(ij)
              	xvec(ij) = xvectmp
!/*  5. Skip out of ij loop and continue ii loop  */
              	EXIT              ! this will exit the innermost do loop
         	ENDIF  ! backward scan ij
         ENDDO 	! loop ij
 		ENDIF   	! forward scan ii
	ENDDO      !  loop ii
         
!   deallocate (atmp)
	 deallocate (binseq)
    DEALLOCATE ( atmprowi, atmprowj, atmpcoli, atmpcolj )

end subroutine arrange




!*******************************************************************************
!  ******** LU Decomposition ************
!*******************************************************************************
subroutine LUDecomp(a,b,row,soln,ans,detmt)
	implicit none
	INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)
	REAL(dbp) , INTENT(OUT) :: detmt

   INTEGER :: i, j, k, l, q(1), p
   REAL(dbp), ALLOCATABLE :: low(:,:), upp(:,:), dia(:,:)
   REAL(dbp), POINTER :: atmp(:,:), btmp(:,:), swap(:), check(:,:)

	detmt = 1.0d0

   ALLOCATE ( low(row,row), upp(row,row), dia(row,soln) )
   ALLOCATE ( atmp(row,row), btmp(row,soln), swap(row+soln), check(row,soln) )
   low(:,:) = 0.0d0; 		upp(:,:) = 0.0d0;    	dia(:,:) = 0.0d0;
   ans(:,:) = 0.0d0
   atmp(:,:) = 0.0d0; 		btmp(:,:) = 0.0d0;      swap(:) = 0.0d0
   check(:,:) = 0.0d0


   atmp(:,:) = a(:,:); 		btmp(:,:) = b(:,:)

   do j = 1,row         ! Upper Matrix has diagonals of ONES
      upp(j,j) = 1.0d0
   end do


! Do the j-th operation to DECOMPOSE    L.U <- A
   do j = 1, row
      write (60,*) 'Operation # ', j

!               ! 1. Finding diagonals for LOWER MATRIX
!      low(j,j) = atmp(j,j)
!      do k = 1, j-1
!         low(j,j) = low(j,j) - low(j,k) * upp(k,j)
!      end do

               ! 2. Finding LOWER MATRIX
      do i = j, row
         low(i,j) = atmp(i,j)
         do k = 1, j-1
            low(i,j) = low(i,j) - low(i,k) * upp(k,j)
         end do
      end do

!  Pivotal strategy - SWAP rows to make pivotal element c[j][j] have the
!  greatest magnitude in its column. This prevents unnecessary division by
!  a small number.
   	q = j - 1 + MAXLOC( dabs(low(j:row,j)) )        ! assigning to the q-array
	   p = q(1)

   	if (p .ne. j) then
   		swap(1:row) = atmp(j,:);   atmp(j,:) = atmp(p,:); 	atmp(p,:) = swap(1:row)
         swap(1:row) = low(j,:);    low(j,:) = low(p,:);   	low(p,:) = swap(1:row)
			swap(1:soln) = btmp(j,:);  btmp(j,:) = btmp(p,:);   	btmp(p,:) = swap(1:soln)

	      detmt = -1.0d0 * detmt
   	end if


! Determinant change signs if rows are swapped (Theorem 3b.)


                	! 1. Finding diagonals for LOWER MATRIX
!      low(j,j) = atmp(j,j)
!      do k = 1, j-1
!         low(j,j) = low(j,j) - low(j,k) * upp(k,j)
!      end do

               ! 2. Finding LOWER MATRIX
      do i = j, row
         low(i,j) = atmp(i,j)
         do k = 1, j-1
            low(i,j) = low(i,j) - low(i,k) * upp(k,j)
         end do
      end do

               ! 3. Finding UPPER MATRIX
      do i = j+1, row
         upp(j,i) = atmp(j,i) / low(j,j)
         do k = 1, j-1
            upp(j,i) = upp(j,i) - low(j,k)*upp(k,i)/low(j,j)
         end do
      end do

   end do      ! end j-loop i.e. Completed LU Decomposition

! Finding Dia Matrix - Forward Substitution
   do l = 1, soln
	   do i = 1, row
   	   dia(i,l) = btmp(i,l) / low(i,i)
      	do k = 1, i-1
         	dia(i,l) = dia(i,l) - low(i,k) * dia(k,l) / low(i,i)
	      end do
   	end do
   end do

! Finding ANSWER Matrix - Backward Substitution
  	ans(:,:) = 0.0d0
  	do l = 1, soln
  		do i = row, 1, -1
      	ans(i,l) = dia(i,l)
      	do j = i+1, row
         	ans(i,l) = ans(i,l) - upp(i,j) * ans(j,l)
	      end do
   	end do
   END DO

   check(:,:) = MATMUL(a,ans)

   WRITE(60,*)'LU Solver Subroutine: Original matrix'
	do j = 1, row
		WRITE(60,500) (a(j,i), i = 1,row)
	end do
   WRITE(60,510)

   WRITE(60,*)'LU Solver Subroutine: Solution matrix'
	do j = 1, row
		WRITE(60,500) (b(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'LU Matrix'
	do j = 1, row
		WRITE(60,500) (low(j,i), i = 1,row), (upp(j,i), i = 1,row)
	end do
   WRITE(60,510)

   WRITE(60,*)'Dia Matrix'
	do j = 1, row
		WRITE(60,500) (dia(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'Result from LU Solver Subroutine'
	do j = 1, row
		WRITE(60,500) (ans(j,i), i = 1,soln)
	end do
   WRITE(60,510)

   WRITE(60,*)'Check - Compare with original solution matrix'
	do j = 1, row
		WRITE(60,500) (check(j,i), i = 1,soln)
	end do
   WRITE(60,510)



   WRITE(60,*)'determinant is ',detmt

   deallocate (low, upp, dia, atmp, btmp, swap, check)

500  FORMAT(300g15.4)
510  FORMAT(//)
520  FORMAT(30i3)


END subroutine LUDecomp



!******************************************************************************
!  Given a matrix A, with logical dimensions M by N and physical dimensions MP by NP, this
!  routine computes its singular value decomposition, A = U.W.V'. The matrix U replaces 
!  A on output. The diagonal matrix of singular values W is output as a vector W. The matrix 
!  V (not the Transpose V') is output as V. M must be greater or equal to N; if it is smaller,
!  then A should be filled up to square with zero rows.
!******************************************************************************
!  Verified: working on 6 Feb 1998
!
! NOTES: 
! 1. When used to solve Linear Equations with n equations and n unknowns,
!     mp=np=m=n.
! 2. When finding inverses, n = soln and b = IdentityMatrix(n)
!
! Modifications to Original Program:
! 1. Equalties/Inequalities are replaced by Differences and compared with EPS
! 2. The Matrix U in U.W.V' is actually stored in "a"
!
! DEBUG:
! 1.  "IMPLICIT NONE" and "REAL(DBP)"
! 2.  Parameters might need to be larger
!******************************************************************************
SUBROUTINE SVDCMP(a,m,n,mp,np,w,v)
   implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
   
   INTEGER :: m,n,mp,np, i,j,l,k, its, nm, jj
   REAL(dbp) :: g, sscale, anorm, s, f,h,x,z,y,c

!   REAL(dbp), parameter :: EPS = 3.0d-15

   INTEGER, PARAMETER :: nmax = 1000        !Maximum anticipated value of N
   REAL(dbp) :: A(mp,np), w(np), v(np,np), rv1(nmax)


!   PRINT *, 'Enter the tolerance or precision required'
!   READ *, eps
   PRINT *, 'Precision chosen as ',eps

   if (m.lt.n) then
   	PRINT *, 'You must augment A with extra zero rows'
      call exit(10)
   ENDIF 

            !Householder Reduction to bidiagonal form
!(see Forsythe,Malcolm,Moler, "Computer Methods for Mathematical Computations"
   g=0.0d0
   sscale = 0.0d0
   anorm = 0.0d0
   do i = 1,n
      l = i + 1
      rv1(i) = sscale*g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if (i.le.m) then
         do k = i,m
            sscale = sscale + dABS(a(k,i))
         end do       ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = i,m
               a(k,i) = a(k,i) / sscale
               s = s + a(k,i)*a(k,i)
            end do    ! k loop
            f = a(i,i)
            g = - SIGN(SQRT(s),f)
            h = f*g - s
            a(i,i) = f - g
            if (i.ne.n) then
               do j = l,n
                  s = 0.0d0
                  do k = i,m
                     s = s + a(k,i)*a(k,j)
                  end do      ! k loop
                  f = s / h
                  do k = i, m 
                     a(k,j) = a(k,j) + f*a(k,i)
                  end do   ! k loop
               end do      ! j loop
            end if
            do k = i, m 
               a(k,i) = sscale * a(k,i)
            end do         ! k loop
         end if
      end if

      w(i) = sscale * g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if ((i.le.m).AND.(i.ne.n)) then
         do k = l, n
            sscale = sscale + dABS(a(i,k))
         end do         ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = l, n
               a(i,k) = a(i,k) /sscale
               s = s + a(i,k) * a(i,k)
            end do      ! k loop 
            f = a(i,l) 
            g = - SIGN(SQRT(s),f)
            h = f * g - s
            a(i,l) = f - g
            do k = l, n
               rv1(k) = a(i,k) / h
            end do      ! k loop
            if (i.ne.m) then
               do j = l, m 
                  s = 0.0d0
                  do k = l, n 
                     s = s + a(j,k)*a(i,k)
                  end do   ! k loop
                  do k = l, n 
                     a(j,k) = a(j,k) + s*rv1(k)
                  end do   ! k loop
               end do      ! j loop
            end if
				do k = l, n
               a(i,k) = sscale * a(i,k)
           	end do
         end if
      end if
      anorm = MAX(anorm, (dABS(w(i)) + dABS(rv1(i))))
   end do

! Accumulation of right-hand Transformations
   do i = n, 1, -1 
      if (i.lt.n) then
!         if (g.ne.0.0d0) then
			if ( dabs(g-0.0d0).gt.EPS ) then
            do j = l, n       ! Double division to avoid possible overflow
               v(j,i) = (a(i,j) / a(i,l)) / g
            end do      ! j loop
            do j = l, n
               s = 0.0d0
               do k = l, n
                  s = s + a(i,k)*v(k,j)
               end do   ! k loop
               do k = l, n
                  v(k,j) = v(k,j) + s * v(k,i)
               end do   ! k loop
           	end do      ! j loop
         end if
         do j = l, n 
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
         end do
      end if
      v(i,i) = 1.0d0
      g = rv1(i)
      l = i
   end do

! Accumulation of left-hand Transformations
   do i = n, 1, -1
      l = 1 + i
      g = w(i)
      if (i.lt.n) then
         do j = l, n
            a(i,j) = 0.0d0
         end do
      end if
!      if (g.ne.0.0d0) then
      if ( dabs(g-0.0d0).gt.EPS ) then
         g = 1.0d0 / g
         if (i.ne.n) then
            do j = l,n 
               s = 0.0d0
               do k = l, m
                  s = s + a(k,i)*a(k,j)
               end do   ! k loop
               f = (s/a(i,i)) * g
               do k = i, m 
                  a(k,j) = a(k,j) + f * a(k,i)
               end do   ! k loop
            end do      ! j loop
         end if
         do j = i, m 
            a(j,i) = a(j,i) * g
         end do         ! j loop
      else
         do j = i, m
            a(j,i) = 0.0d0
         end do         ! j loop
      end if
      a(i,i) = a(i,i) + 1.0d0
   end do               ! i loop

! Diagonalization of the bidigonal form
   do k = n, 1, -1                  !Loop over singular values
      do its = 1,30                 !Loop over allowed iterations
         do l = k, 1, -1            !Test for splitting
            nm = l - 1              ! Note that rv1(1) is always zero
!           if ( (dABS(rv1(l))+anorm) .eq. anorm ) GO TO 2
!          	if ( (dABS(w(nm))+anorm) .eq. anorm ) GO TO 1
            if ( dabs((dABS(rv1(l))+anorm) - anorm).lt.eps ) GO TO 2
          	if ( dabs((dABS(w(nm))+anorm) - anorm).lt.eps ) GO TO 1
         end do      !  l loop

1        c = 0.0d0                  ! Cancellation of rv1(l), if l>1 :
         s = 1.0d0
         do i = l, k
            f = s * rv1(i)
!            if ( (dABS(f)+anorm) .ne. anorm ) then
            if ( dabs( (dABS(f)+anorm) - anorm) .GT. eps ) then

               g = w(i)
               h = SQRT(f*f + g*g)
               w(i) = h
               h = 1.0d0 / h
               c = g * h
               s = -f * h
               do j = 1, m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
               end do   ! j loop
            end if
         end do         ! i loop
2        z = w(k) 
         if (l .eq. k) then         ! convergence
				if (z .lt. 0.0d0) then  ! Singular value is made non-negative
               w(k) = -z
               do j = 1,n
                  v(j,k) = -v(j,k)
               end do         ! j loop
	    		end if
            GO TO 3
         end if
         if (its.eq.30) then
         	PRINT*, 'No Convergence in 30 iterations'
            call exit(10)
         ENDIF
         x = w(l)          ! Shift from bottom 2-by-2 minor
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0d0*h*y)
         g = SQRT(f*f + 1.0d0)
         f = ( (x-z)*(x+z) + h*((y/(f+SIGN(g,f))) - h) ) / x

! Next   QR Transformation
         c = 1.0d0
         s = 1.0d0
         do j = l, nm
          	i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = SQRT(f*f + h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, n
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) = (x*c) + (z*s)
               v(jj,i) = -(x*s) + (z*c)
           	end do
            z = SQRT(f*f + h*h)
            w(j) = z
!            if (z.ne.0.0d0) then
            if (  dabs(z-0.0d0).gt.eps  ) then
               z = 1.0d0 / z
               c = f*z
               s = h*z
            end if
            f = (g*c) + (y*s)
            x = -(g*s) + (y*c)
            do jj = 1, m
               y = a(jj,j)
               z = a(jj,i)
               a(jj,j) = (y*c) + (z*s)
               a(jj,i) = -(y*s) + (z*c)
            end do
         end do         ! j loop
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
      end do            ! its loop
3  continue
   end do               ! k loop

   return
END SUBROUTINE svdcmp



!*********************************************************************
! Added by C.Chee Feb 1998
! NOTE:
! 1. Given matrix B, find answer matrix ANS where A.ANS = B.
!     and "soln" is the number of columns of Solution Matrix B.
! DEBUG:
! 1.  "IMPLICIT NONE" and "REAL(DBP)"
!*********************************************************************
subroutine svdsolve(a,b,row,soln,ans)
	implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)

   INTEGER, INTENT(IN) :: row,soln
	REAL(dbp) , INTENT(IN) :: a(row,row), b(row,soln)
	REAL(dbp) , INTENT(OUT) :: ans(row,soln)

   REAL(dbp), allocatable :: usvd(:,:), wsvd(:), vsvd(:,:), winvutb(:,:)
   INTEGER :: ii,ij
   REAL(dbp) :: wmax, wmin, factor1
   CHARACTER :: ask2*80


   REAL(dbp) :: tmpmat(row,row), wmat(row,row), wmatinv(row, row)

   ALLOCATE ( usvd(row,row), wsvd(row), vsvd(row,row), winvutb(row,soln) )



   usvd(:,:) = a(:,:)
   vsvd(:,:) = 0.0d0
   wsvd(:) = 0.0d0

! Produces a SINGULAR VALUE DECOMPOSITION of matrix a
   call SVDCMP(usvd, row,row, row,row, wsvd, vsvd)

   wmat = identity_matrix(row)
   do ij = 1,row
   	wmat(ij,ij) = wsvd(ij)
   end do
   tmpmat = a - MATMUL(usvd, MATMUL(wmat,TRANSPOSE(vsvd)) )
   WRITE(60,*) 'Here is it'
   do ii = 1, row
   	WRITE(60,333) (tmpmat(ii,ij), ij=1,row)
   ENDDO
333 FORMAT(300g15.3)



! KEY TO SVD - If Diagonal elements are too small then MAKE IT ZERO!!!!!!!
   PRINT *, 'Solution by Singular Value Decomposition'
   PRINT *, 'Is matrix expected to be close to singular (y/n) ?'
   READ(*, *) ask2

   if (ask2 .eq. 'y') then
	   wmax = MAXVAL(dabs(wsvd(:)))
   	PRINT *, 'enter precision for svd (wmax = ',wmax, ')'
	   READ *, factor1
   	wmin = factor1 * eps * wmax          ! 10 times the set precision
! Actually DON'T set w(i,i) = 0, the fact that it's equivalent to zero lies in the
! next few lines in Solving for Ax=b or x=V.[1/w(i,i)].U`.b
		PRINT *, 'w-min threshold   Precision'
   	PRINT *, wmin, eps
   else
      wmin = 0.0d0
   end if


   wmatinv(:,:) = 0.0d0
   do ij = 1,row                       ! modified 9/9/98
     	if ( dabs(wsvd(ij)) .ge.wmin) 	wmatinv(ij,ij) = 1.0/wsvd(ij)
   end do
   winvutb = MATMUL(wmatinv, MATMUL(TRANSPOSE(usvd), b))

   ans(:,:) = 0.0d0
   ans(:,:) = MATMUL(vsvd,winvutb)

   deallocate (usvd, wsvd, vsvd, winvutb)

   return
END subroutine svdsolve




end module gaussm3


