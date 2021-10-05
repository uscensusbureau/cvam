!#####################################################################
module math_R
   ! Access to R's C API mathematical functions through F77 and C, 
   ! using the functions in mathRF77.f and mathRC.c
   use program_constants
   use error_handler
   implicit none
   private ! by default
   ! list public functions and subroutines
   public get_randgen_state_R, put_randgen_state_R, &
        dunif_R, punif_R, qunif_R, runif_R, &
        dnorm_R, pnorm_R, qnorm_R, rnorm_R, &
        dchisq_R, pchisq_R, qchisq_R, rchisq_R, &
        dt_R, pt_R, qt_R, rt_R, &
        dbeta_R, pbeta_R, qbeta_R, rbeta_R, &
        dgamma_R, pgamma_R, qgamma_R, rgamma_R, &
        dbinom_R, dbinom_raw_R, pbinom_R, qbinom_R, rbinom_R, &
        dpois_R, dpois_raw_R, ppois_R, qpois_R, rpois_R, &
        dnbinom_R, pnbinom_R, qnbinom_R, rnbinom_R, &
        lgamma_R, digamma_R, trigamma_R, &
        log1p_R, log1pmx_R, log1pexp_R, expm1_R, lgamma1p_R
   ! private parameters
   character(len=*), parameter :: modname = "math_R"
   !##################################################################
contains
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function get_randgen_state_R(err) result(answer)
      !	must be called before random number generation
      implicit none
      ! Required arguments
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_randgen_state_R"
      ! begin
      answer = RETURN_FAIL
      call getRandgenStateRF77()
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function get_randgen_state_R
   !##################################################################
   integer(our_int) function put_randgen_state_R(err) result(answer)
      !	must be called after random number generation
      implicit none
      ! Required arguments
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "put_randgen_state_R"
      ! begin
      answer = RETURN_FAIL
      call putRandgenStateRF77()
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function put_randgen_state_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dunif_R(x, ans, err, &
      a, b, logd) result(answer)
      !	uniform probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   a = lower limit of support (default is 0.D0)
      !   b = upper limit of support (default is 1.D0)
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: a, b
      logical, intent(in), optional :: logd
      ! locals
      real(kind=our_dble) :: a_local, b_local
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dunif_R"
      ! begin
      answer = RETURN_FAIL
      if( present(a) ) then
         a_local = a
      else
         a_local = 0.D0
      end if
      if( present(b) ) then
         b_local = b
      else
         b_local = 1.D0
      end if
      if( b_local <= a_local ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dunifRF77(x, a_local, b_local, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Upper limit b does not exceed lower limit a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dunif_R
   !##################################################################
   integer(our_int) function punif_R(q, ans, err, &
      a, b, lowertail, logp) result(answer)
      !	uniform cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   a = lower limit of support (default is 0.D0)
      !   b = upper limit of support (default is 1.D0)
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: a, b
      logical, intent(in), optional :: lowertail, logp
      ! locals
      real(kind=our_dble) :: a_local, b_local
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "punif_R"
      ! begin
      answer = RETURN_FAIL
      if( present(a) ) then
         a_local = a
      else
         a_local = 0.D0
      end if
      if( present(b) ) then
         b_local = b
      else
         b_local = 1.D0
      end if
      if( b_local <= a_local ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call punifRF77(q, a_local, b_local, &
           lowertail_local, logp_local, ans)
      ! normalexit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Upper limit b does not exceed lower limit a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function punif_R
   !##################################################################
   integer(our_int) function qunif_R(p, ans, err, &
      a, b, lowertail, logp ) result(answer)
      !	uniform quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   a = lower limit of support (default is 0.D0)
      !   b = upper limit of support (default is 1.D0)
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: a, b
      logical, intent(in), optional :: lowertail, logp
      ! locals
      real(kind=our_dble) :: a_local, b_local
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qunif_R"
      ! begin
      answer = RETURN_FAIL
      if( present(a) ) then
         a_local = a
      else
         a_local = 0.D0
      end if
      if( present(b) ) then
         b_local = b
      else
         b_local = 1.D0
      end if
      if( b_local <= a_local ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qunifRF77(p, a_local, b_local, lowertail_local, &
           logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Upper limit b does not exceed lower limit a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function qunif_R
   !##################################################################
   integer(our_int) function runif_R(ans, err, &
        a, b) result(answer)
      !	uniform random variate
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      real(kind=our_dble), intent(in), optional :: a, b
      ! locals
      real(kind=our_dble) :: a_local, b_local
      character(len=*), parameter :: subname = "runif_R"
      ! begin
      answer = RETURN_FAIL
      if( present(a) ) then
         a_local = a
      else
         a_local = 0.D0
      end if
      if( present(b) ) then
         b_local = b
      else
         b_local = 1.D0
      end if
      if( b_local <= a_local ) goto 100
      call runifRF77(a_local, b_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Upper limit b does not exceed lower limit a" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function runif_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dnorm_R(x, ans, err, &
      mu, sigma, logd) result(answer)
      !	normal probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   mu = mean (default is 0.D0)
      !   sigma = standard deviation (default is 1.D0)
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: mu, sigma
      logical, intent(in), optional :: logd
      ! locals
      real(kind=our_dble) :: mu_local, sigma_local
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dnorm_R"
      ! begin
      answer = RETURN_FAIL
      if( present(mu) ) then
         mu_local = mu
      else
         mu_local = 0.D0
      end if
      if( present(sigma) ) then
         if( sigma <= 0.D0 ) goto 100 
         sigma_local = sigma
      else
         sigma_local = 1.D0
      end if
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dnormRF77(x, mu_local, sigma_local, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument sigma is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dnorm_R
   !##################################################################
   integer(our_int) function pnorm_R(q, ans, err, &
      mu, sigma, lowertail, logp) result(answer)
      !	normal cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   mu = mean (default is 0.D0)
      !   sigma = standard deviation (default is 1.D0)
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: mu, sigma
      logical, intent(in), optional :: lowertail, logp
      ! locals
      real(kind=our_dble) :: mu_local, sigma_local
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pnorm_R"
      ! begin
      answer = RETURN_FAIL
      if( present(mu) ) then
         mu_local = mu
      else
         mu_local = 0.D0
      end if
      if( present(sigma) ) then
         if( sigma <= 0.D0 ) goto 100 
         sigma_local = sigma
      else
         sigma_local = 1.D0
      end if
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pnormRF77(q, mu_local, sigma_local, &
           lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument sigma is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pnorm_R
   !##################################################################
   integer(our_int) function qnorm_R(p, ans, err, &
      mu, sigma, lowertail, logp ) result(answer)
      !	normal quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   mu = mean (default is 0.D0)
      !   sigma = standard deviation (default is 1.D0)
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      real(kind=our_dble), intent(in), optional :: mu, sigma
      logical, intent(in), optional :: lowertail, logp
      ! locals
      real(kind=our_dble) :: mu_local, sigma_local
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qnorm_R"
      ! begin
      answer = RETURN_FAIL
      if( present(mu) ) then
         mu_local = mu
      else
         mu_local = 0.D0
      end if
      if( present(sigma) ) then
         if( sigma <= 0.D0 ) goto 100 
         sigma_local = sigma
      else
         sigma_local = 1.D0
      end if
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qnormRF77(p, mu_local, sigma_local, lowertail_local, &
           logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument sigma is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function qnorm_R
   !##################################################################
   integer(our_int) function rnorm_R(ans, err, &
        mu, sigma) result(answer)
      !	normal random variate
      ! Required args
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   mu = mean (default is 0.D0)
      !   sigma = standard deviation (default is 1.D0)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      real(kind=our_dble), intent(in), optional :: mu, sigma
      ! locals
      real(kind=our_dble) :: mu_local, sigma_local
      character(len=*), parameter :: subname = "rnorm_R"
      ! begin
      answer = RETURN_FAIL
      if( present(mu) ) then
         mu_local = mu
      else
         mu_local = 0.D0
      end if
      if( present(sigma) ) then
         if( sigma <= 0.D0 ) goto 100 
         sigma_local = sigma
      else
         sigma_local = 1.D0
      end if
      call rnormRF77(mu_local, sigma_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument sigma is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function rnorm_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dchisq_R(x, df, ans, err, &
      logd) result(answer)
      !	chisquare probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dchisq_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dchisqRF77(x, df, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dchisq_R
   !##################################################################
   integer(our_int) function pchisq_R(q, df, ans, err, &
      lowertail, logp) result(answer)
      !	chisquare cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pchisq_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pchisqRF77(q, df, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pchisq_R
   !##################################################################
   integer(our_int) function qchisq_R(p, df, ans, err, &
      lowertail, logp ) result(answer)
      !	chisquare quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qchisq_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qchisqRF77(p, df, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qchisq_R
   !##################################################################
   integer(our_int) function rchisq_R(df, ans, err) result(answer)
      !	chisquare random variate
      ! Required args
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rchisq_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      call rchisqRF77(df, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rchisq_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dt_R(x, df, ans, err, &
      logd) result(answer)
      !	Student's t probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dt_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dtRF77(x, df, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dt_R
   !##################################################################
   integer(our_int) function pt_R(q, df, ans, err, &
      lowertail, logp) result(answer)
      !	Student's t cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pt_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call ptRF77(q, df, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pt_R
   !##################################################################
   integer(our_int) function qt_R(p, df, ans, err, &
      lowertail, logp ) result(answer)
      !	Student's t quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qt_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qtRF77(p, df, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qt_R
   !##################################################################
   integer(our_int) function rt_R(df, ans, err) result(answer)
      !	Student's t random variate
      ! Required args
      !   df = degrees of freedom
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rt_R"
      ! begin
      answer = RETURN_FAIL
      if( df <= 0.D0 ) goto 100
      call rtRF77(df, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument df is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rt_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dbeta_R(x, shape1, shape2, ans, err, &
      logd) result(answer)
      !	beta probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   shape1, shape2 = shape parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, shape1, shape2
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dbeta_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape1 <= 0.D0 ) .or. ( shape2 <= 0.D0 ) ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dbetaRF77(x, shape1, shape2, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape1 or shape2 is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dbeta_R
   !##################################################################
   integer(our_int) function pbeta_R(q, shape1, shape2, ans, err, &
      lowertail, logp) result(answer)
      !	beta cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   shape1, shape2 = shape parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(in) :: shape1, shape2
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pbeta_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape1 <= 0.D0 ) .or. ( shape2 <= 0.D0 ) ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pbetaRF77(q, shape1, shape2, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape1 or shape2 is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pbeta_R
   !##################################################################
   integer(our_int) function qbeta_R(p, shape1, shape2, ans, err, &
      lowertail, logp ) result(answer)
      !	beta quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   shape1, shape2 = shape parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(in) :: shape1, shape2
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qbeta_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape1 <= 0.D0 ) .or. ( shape2 <= 0.D0 ) ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qbetaRF77(p, shape1, shape2, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape1 or shape2 is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qbeta_R
   !##################################################################
   integer(our_int) function rbeta_R(shape1, shape2, ans, err) result(answer)
      !	beta random variate
      ! Required args
      !   shape1, shape2 = shape parameters
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: shape1, shape2
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rbeta_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape1 <= 0.D0 ) .or. ( shape2 <= 0.D0 ) ) goto 100
      call rbetaRF77(shape1, shape2, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape1 or shape2 is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rbeta_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dgamma_R(x, shape, scale, ans, err, &
      logd) result(answer)
      !	gamma probability density function
      ! Required args
      !   x = value at which density is to be computed
      !   shape, scale = shape, scale parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, shape, scale
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dgamma_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape <= 0.D0 ) .or. ( scale <= 0.D0 ) ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dgammaRF77(x, shape, scale, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape or scale is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dgamma_R
   !##################################################################
   integer(our_int) function pgamma_R(q, shape, scale, ans, err, &
      lowertail, logp) result(answer)
      !	gamma cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   shape, scale = shape, scale parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(in) :: shape, scale
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pgamma_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape <= 0.D0 ) .or. ( scale <= 0.D0 ) ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pgammaRF77(q, shape, scale, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape or scale is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pgamma_R
   !##################################################################
   integer(our_int) function qgamma_R(p, shape, scale, ans, err, &
      lowertail, logp ) result(answer)
      !	gamma quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   shape, scale = shape, scale parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(in) :: shape, scale
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qgamma_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape <= 0.D0 ) .or. ( scale <= 0.D0 ) ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qgammaRF77(p, shape, scale, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape or scale is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qgamma_R
   !##################################################################
   integer(our_int) function rgamma_R(shape, scale, ans, err) result(answer)
      !	gamma random variate
      ! Required args
      !   shape, scale = shape, scale parameters
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: shape, scale
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rgamma_R"
      ! begin
      answer = RETURN_FAIL
      if( ( shape <= 0.D0 ) .or. ( scale <= 0.D0 ) ) goto 100
      call rgammaRF77(shape, scale, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument shape or scale is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rgamma_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dbinom_R(x, n, prob, ans, err, &
      logd) result(answer)
      !	binom probability mass function, gives warning if x and n
      ! are not integers
      ! Required args
      !   x = value at which mass is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dbinomRF77(x, n, prob, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dbinom_R
   !##################################################################
   integer(our_int) function dbinom_raw_R(x, n, probp, probq, ans, err, &
      logd) result(answer)
      !	binom probability mass function, works for noninteger n, p
      ! Required args
      !   x = value at which mass is to be computed
      !   n, probp = n, prob parameters
      !   probq = 1.D0 - probp
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, n, probp, probq
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dbinom_raw_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( probp < 0.D0 ) .or. ( probp > 1.D0 ) ) goto 200
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dbinomrawRF77(x, n, probp, probq, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dbinom_raw_R
   !##################################################################
   integer(our_int) function pbinom_R(q, n, prob, ans, err, &
      lowertail, logp) result(answer)
      !	binom cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pbinomRF77(q, n, prob, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pbinom_R
   !##################################################################
   integer(our_int) function qbinom_R(p, n, prob, ans, err, &
      lowertail, logp ) result(answer)
      !	binom quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qbinomRF77(p, n, prob, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qbinom_R
   !##################################################################
   integer(our_int) function rbinom_R(n, prob, ans, err) result(answer)
      !	binom random variate
      ! Required args
      !   n, probp = n, prob parameters
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200

      call rbinomRF77(n, prob, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rbinom_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dpois_R(x, lambda, ans, err, &
      logd) result(answer)
      !	poisson probability mass function, gives zero for noninteger x
      ! Required args
      !   x = value at which density is to be computed
      !   lambda = mean
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, lambda
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dpois_R"
      ! begin
      answer = RETURN_FAIL
      if( lambda <= 0.D0 ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dpoisRF77(x, lambda, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument lambda is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dpois_R
   !##################################################################
   integer(our_int) function dpois_raw_R(x, lambda, ans, err, &
      logd) result(answer)
      !	poisson probability mass function, continuous version
      ! Required args
      !   x = value at which density is to be computed
      !   lambda = mean
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, lambda
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dpois_raw_R"
      ! begin
      answer = RETURN_FAIL
      if( lambda <= 0.D0 ) goto 100
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dpoisrawRF77(x, lambda, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument lambda is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function dpois_raw_R
   !##################################################################
   integer(our_int) function ppois_R(q, lambda, ans, err, &
      lowertail, logp) result(answer)
      !	poisson cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   lambda = mean
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q
      real(kind=our_dble), intent(in) :: lambda
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "ppois_R"
      ! begin
      answer = RETURN_FAIL
      if( lambda <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call ppoisRF77(q, lambda, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument lambda is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function ppois_R
   !##################################################################
   integer(our_int) function qpois_R(p, lambda, ans, err, &
      lowertail, logp ) result(answer)
      !	poisson quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   lambda = mean
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p
      real(kind=our_dble), intent(in) :: lambda
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qpois_R"
      ! begin
      answer = RETURN_FAIL
      if( lambda <= 0.D0 ) goto 100
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qpoisRF77(p, lambda, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument lambda is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qpois_R
   !##################################################################
   integer(our_int) function rpois_R(lambda, ans, err) result(answer)
      !	poisson random variate
      ! Required args
      !   lambda = mean
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: lambda
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rpois_R"
      ! begin
      answer = RETURN_FAIL
      if( lambda <= 0.D0 ) goto 100
      call rpoisRF77(lambda, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument lambda is non-positive" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rpois_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function dnbinom_R(x, n, prob, ans, err, &
      logd) result(answer)
      !	negative binom probability mass function, gives warning and 
      ! returns zero if x is not integer
      ! Required args
      !   x = value at which mass is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   logd = if .true., log density is computed (defaults to .false.)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: logd
      ! locals
      integer(kind=our_int) :: logd_local
      character(len=*), parameter :: subname = "dnbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(logd) ) then
         if( logd ) then
            logd_local = 1
         else
            logd_local = 0
         end if
      else
         logd_local = 0
      end if
      call dnbinomRF77(x, n, prob, logd_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function dnbinom_R
   !##################################################################
   integer(our_int) function pnbinom_R(q, n, prob, ans, err, &
      lowertail, logp) result(answer)
      !	neg binom cumulative distribution function
      ! Required args
      !   q = value at which cdf is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is computed; if 0, upper
      !      (default is T)
      !   logp = if T, log-probability is computed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: q, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "pnbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call pnbinomRF77(q, n, prob, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
   end function pnbinom_R
   !##################################################################
   integer(our_int) function qnbinom_R(p, n, prob, ans, err, &
      lowertail, logp ) result(answer)
      !	neg binom quantile function
      ! Required args
      !   p = prob at which quantile is to be computed
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      ! Optional args
      !   lowertail = if T, lower tail is assumed; if F, upper
      !      (default is T)
      !   logp = if T, log-probability is assumed
      !      (default is F)
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: p, n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! Optional arguments
      logical, intent(in), optional :: lowertail, logp
      ! locals
      integer(kind=our_int) :: lowertail_local, logp_local
      character(len=*), parameter :: subname = "qnbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      if( present(lowertail) ) then
         if( lowertail ) then
            lowertail_local = 1
         else
            lowertail_local = 0
         end if
      else
         lowertail_local = 1
      end if
      if( present(logp) ) then
         if( logp ) then
            logp_local = 1
         else
            logp_local = 0
         end if
      else
         logp_local = 0
      end if
      call qnbinomRF77(p, n, prob, lowertail_local, logp_local, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function qnbinom_R
   !##################################################################
   integer(our_int) function rnbinom_R(n, prob, ans, err) result(answer)
      !	neg binom random variate
      ! Required args
      !   n, prob = n, prob parameters
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: n, prob
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "rnbinom_R"
      ! begin
      answer = RETURN_FAIL
      if( n <= 0.D0 ) goto 100
      if( ( prob < 0.D0 ) .or. ( prob > 1.D0 ) ) goto 200
      call rnbinomRF77(n, prob, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
100   call err_handle(err, 1, &
           comment = "Argument n is non-positive" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Argument prob is outside of [0,1]" )
      goto 800
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function rnbinom_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function lgamma_R(x, ans, err ) result(answer)
      !	logarithm of (absolute value of) gamma function
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "lgamma_R"
      ! begin
      answer = RETURN_FAIL
      call lgammaRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function lgamma_R
   !##################################################################
   integer(our_int) function digamma_R(x, ans, err ) result(answer)
      !	first deriv of log gamma function
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "digamma_R"
      ! begin
      answer = RETURN_FAIL
      call digammaRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function digamma_R
   !##################################################################
   integer(our_int) function trigamma_R(x, ans, err ) result(answer)
      !	second deriv of log gamma function
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "trigamma_R"
      ! begin
      answer = RETURN_FAIL
      call trigammaRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function trigamma_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
   integer(our_int) function log1p_R(x, ans, err ) result(answer)
      !	Computes log(1 + x) (log 1 plus x), accurately even for 
      ! small x, i.e., abs(x) << 1
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "log1p_R"
      ! begin
      answer = RETURN_FAIL
      call log1pRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function log1p_R
   !##################################################################
   integer(our_int) function log1pmx_R(x, ans, err ) result(answer)
      !	Computes log(1 + x) - x, accurately even for 
      ! small x, i.e., abs(x) << 1
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "log1pmx_R"
      ! begin
      answer = RETURN_FAIL
      call log1pmxRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function log1pmx_R
   !##################################################################
   integer(our_int) function log1pexp_R(x, ans, err ) result(answer)
      !	Computes log(1 + exp(x)), accurately even for 
      ! large x, i.e., x > 720
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "log1pexp_R"
      ! begin
      answer = RETURN_FAIL
      call log1pexpRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function log1pexp_R
   !##################################################################
   integer(our_int) function expm1_R(x, ans, err ) result(answer)
      !	Computes exp(x) - 1 (exp x minus 1 ), accurately even for
      ! small x, i.e. abs(x) << 1
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "expm1_R"
      ! begin
      answer = RETURN_FAIL
      call expm1RF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function expm1_R
   !##################################################################
   integer(our_int) function lgamma1p_R(x, ans, err ) result(answer)
      !	Computes log(gamma(x + 1)), accurately even for
      ! small x, i.e. 0 < x < 0.5
      ! Required args
      !   x = value at which function is to be computed
      !   ans = result
      !   err =  error handler object
      implicit none
      ! Required arguments
      real(kind=our_dble), intent(in) :: x
      real(kind=our_dble), intent(out) :: ans
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "lgamma1p_R"
      ! begin
      answer = RETURN_FAIL
      call lgamma1pRF77(x, ans)
      ! normal exit
      answer = RETURN_SUCCESS
      ! added to suppress warning about unused dummy argument err
      if( answer == RETURN_FAIL ) goto 800
      return
      ! error traps
800   call err_handle(err, 2,  whichsub = subname, whichmod = modname )
    end function lgamma1p_R
   !##################################################################
   !##################################################################
   !##################################################################
   !##################################################################
end module math_R
!#####################################################################
