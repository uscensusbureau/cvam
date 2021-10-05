!#####################################################################
module math_funcs
   ! Various mathematical functions that build on those available
   ! from the R API
   use program_constants
   use error_handler
   use math_R
   implicit none
   private ! by default
   ! list public functions and subroutines
   public llincgamma, llincgamma_d1, llincgamma_d2
   ! private parameters
   character(len=*), parameter :: modname = "math"
   !##################################################################
contains
   !##################################################################
   integer(our_int) function llincgamma(a, x, ans, err) result(answer)
      ! logarithm of lower incomplete gamma function
      ! Required arguments
      !    a = argument to gamma function (must be positive)
      !    x = upper limit of integration (must be positive)
      !    ans = result
      !    err = error handler object
      implicit none
      ! required args
      real(kind=our_dble), intent(in) :: a
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: log_cdf, log_gamma
      character(len=*), parameter :: subname = "llincgamma"
      ! begin
      answer = RETURN_FAIL
      if( ( a <= 0.D0 ) .or. ( x <= 0.D0 ) ) goto 100
      if( pgamma_R(x, a, 1.D0, log_cdf, err, logp = .true. ) &
           == RETURN_FAIL ) goto 800
      if( lgamma_R(a, log_gamma, err ) == RETURN_FAIL ) goto 800
      ans = log_cdf + log_gamma
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument a or x is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function llincgamma
   !##################################################################
   integer(our_int) function llincgamma_deriv_1(a, x, h, ans, err) &
        result(answer) 
      ! First derivative of the log-lower incomplete gamma function
      ! with respect to a, approximated by two-point central difference
      ! formula
      ! Required arguments
      !    a = argument to gamma function (not checked, assumed positive)
      !    x = upper limit of integration (not checked, assumed positive)
      !    h = step size (not checked for validity)
      !    ans = result
      !    err = error handler object
      implicit none
      ! required args
      real(kind=our_dble), intent(in) :: a
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(in) :: h
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! optional
      ! locals
      real(kind=our_dble) :: f_plus, f_minus
      character(len=*), parameter :: subname = "llincgamma_deriv_1"
      ! begin
      answer = RETURN_FAIL
      if( llincgamma(a+h, x, f_plus, err) == RETURN_FAIL ) goto 800 
      if( llincgamma(a-h, x, f_minus, err) == RETURN_FAIL ) goto 800 
      ans = ( f_plus - f_minus ) / ( 2.D0 * h )
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function llincgamma_deriv_1
   !##################################################################
   integer(our_int) function llincgamma_deriv_2(a, x, h, ans, err) &
        result(answer) 
      ! Second derivative of the log-lower incomplete gamma function
      ! with respect to a, approximated by three-point central difference
      ! formula
      ! Required arguments
      !    a = argument to gamma function (must be positive, not checked)
      !    x = upper limit of integration (must be positive, not checked)
      !    h = step size (not checked)
      !    ans = result
      !    err = error handler object
      ! Optional arguments
      implicit none
      ! required args
      real(kind=our_dble), intent(in) :: a
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(in) :: h
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! optional
      ! locals
      real(kind=our_dble) :: f_plus, f, f_minus
      character(len=*), parameter :: subname = "llincgamma_deriv_2"
      ! begin
      answer = RETURN_FAIL
      if( llincgamma(a+h, x, f_plus, err) == RETURN_FAIL ) goto 800 
      if( llincgamma(a, x, f, err) == RETURN_FAIL ) goto 800 
      if( llincgamma(a-h, x, f_minus, err) == RETURN_FAIL ) goto 800 
      ans = ( f_plus - 2.D0*f + f_minus ) / ( h**2 )
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function llincgamma_deriv_2
   !##################################################################
   integer(our_int) function llincgamma_d1(a, x, ans, err, h, richardson) &
        result(answer) 
      ! Numerical 1st derivative of the log-lower incomplete gamma function
      ! with respect to a
      ! Required arguments
      !    a = argument to gamma function (must be positive)
      !    x = upper limit of integration (must be positive)
      !    ans = result
      !    err = error handler object
      ! Optional arguments
      !    h = step size (default is min(0.1, a/2)
      !    richardson = if .true., uses richardson extrapolation
      !       (default is .true.)
      implicit none
      ! required args
      real(kind=our_dble), intent(in) :: a
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! optional
      real(kind=our_dble), intent(in), optional :: h
      logical, intent(in), optional :: richardson
      ! locals
      real(kind=our_dble) :: h_local, dA, dB
      logical :: richardson_local
      character(len=*), parameter :: subname = "llincgamma_d1"
      ! begin
      answer = RETURN_FAIL
      if( ( a <= 0.D0 ) .or. ( x <= 0.D0 ) ) goto 100
      if( present(h) ) then
         h_local = h
      else
         h_local = min( 0.1D0, a/2.D0 )
      end if
      if( h_local <= 0.D0 ) goto 150 
      if( h_local > a ) goto 200 
      if( present(richardson) ) then
         richardson_local = richardson
      else
         richardson_local = .true.
      end if
      if( richardson_local ) then
         if( llincgamma_deriv_1( a, x, h_local, dA, err ) &
              == RETURN_FAIL ) goto 800
         if( llincgamma_deriv_1( a, x, h_local/2.D0, dB, err ) &
              == RETURN_FAIL ) goto 800
         ans = ( 4.D0 * dB - dA ) / 3.D0
      else
         if( llincgamma_deriv_1( a, x, h_local, ans, err ) &
              == RETURN_FAIL ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument a or x is non-positive" )
      goto 800
150   call err_handle(err, 1, &
           comment = "Step size h is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Step size h is too large, exceeds a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function llincgamma_d1
   !##################################################################
   integer(our_int) function llincgamma_d2(a, x, ans, err, h, richardson) &
        result(answer) 
      ! Numerical 2nd derivative of the log-lower incomplete gamma function
      ! with respect to a
      ! Required arguments
      !    a = argument to gamma function (must be positive)
      !    x = upper limit of integration (must be positive)
      !    ans = result
      !    err = error handler object
      ! Optional arguments
      !    h = step size (default is min(0.1, a/2)
      !    richardson = if .true., uses richardson extrapolation
      !       (default is .true.)
      implicit none
      ! required args
      real(kind=our_dble), intent(in) :: a
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! optional
      real(kind=our_dble), intent(in), optional :: h
      logical, intent(in), optional :: richardson
      ! locals
      real(kind=our_dble) :: h_local, dA, dB
      logical :: richardson_local
      character(len=*), parameter :: subname = "llincgamma_d2"
      ! begin
      answer = RETURN_FAIL
      if( ( a <= 0.D0 ) .or. ( x <= 0.D0 ) ) goto 100
      if( present(h) ) then
         h_local = h
      else
         h_local = min( 0.1D0, a/2.D0 )
      end if
      if( h_local <= 0.D0 ) goto 150 
      if( h_local > a ) goto 200 
      if( present(richardson) ) then
         richardson_local = richardson
      else
         richardson_local = .true.
      end if
      if( richardson_local ) then
         if( llincgamma_deriv_2( a, x, h_local, dA, err ) &
              == RETURN_FAIL ) goto 800
         if( llincgamma_deriv_2( a, x, h_local/2.D0, dB, err ) &
              == RETURN_FAIL ) goto 800
         ans = ( 4.D0 * dB - dA ) / 3.D0
      else
         if( llincgamma_deriv_2( a, x, h_local, ans, err ) &
              == RETURN_FAIL ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument a or x is non-positive" )
      goto 800
150   call err_handle(err, 1, &
           comment = "Step size h is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Step size h is too large, exceeds a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function llincgamma_d2
   !##################################################################
 end module math_funcs
!#####################################################################
