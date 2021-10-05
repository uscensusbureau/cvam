!#####################################################################
module cvam_engine
   ! data structures and computational routines for CVAM modeling
   use program_constants
   use error_handler
   use dynalloc
   use quick_sort
   use tabulate
   use matrix_methods
   use math_R
   use math_funcs
   implicit none
   private ! by default
   ! declare public types (contents are still private)
   public :: workspace_type_cvam
   ! declare public parameters
   public :: model_type_str_len
   ! declare public procedures
   public :: nullify_workspace_type_cvam, run_cvam_model, &
        run_cvam_estimate_em, run_cvam_predict_em, &
        run_cvam_impute_freq, run_cvam_impute_microdata, &
        run_cvam_lik, run_mlogit, run_mlogit_loglik_derivs, &
        run_lcprev_loglik_derivs
   ! parameters private to this module
   real(kind=our_dble), parameter :: &
        log_huge = log( huge( real(0,kind=our_dble) ) ), &
        log_tiny = log( tiny( real(0,kind=our_dble) ) )
   character(len=*), parameter :: modname = "cvam_engine"
   ! public parameters
   integer(kind=our_int), parameter :: model_type_str_len = 30, &
        method_str_len = 30, estimate_type_str_len = 30
   ! private data types
   !##################################################################
   type :: workspace_type_int_array_1d
      ! for building ragged arrays
      integer(kind=our_int), allocatable :: vec(:)
   end type workspace_type_int_array_1d
   !##################################################################
   type :: workspace_type_int_array_2d
      type(workspace_type_int_array_1d), allocatable :: vec(:)
   end type workspace_type_int_array_2d
   !##################################################################
   type :: workspace_type_int_array_3d
      type(workspace_type_int_array_2d), allocatable :: vec(:)
   end type workspace_type_int_array_3d
   !##################################################################
   type :: workspace_type_cvam_estimate
      ! stores information pertaining to a set of estimated cell
      ! probabilities (marginal and/or conditional)
      ! the structure of this workspace is very similar to that of
      ! a submodel
      integer(kind=our_int) :: estimate_type_int = 0  ! 
      character(len=estimate_type_str_len) :: estimate_type = ""
      integer(kind=our_int) :: estimate_number = 0 ! unique id
      integer(kind=our_int) :: nvar = 0 ! number of variables in estimate
      !### these arrays have length nvar
      ! number of levels this variable (length=nvar):
      integer(kind=our_int), allocatable :: nlev(:) 
      ! which position (column number) in the full model (input data)
      integer(kind=our_int), allocatable :: model_posn(:) 
      ! variable treated as fixed in this estimate? (length=nvar):
      logical, allocatable :: fixed(:)
      !### these arrays have length ncol_input_data
      ! is a model variable present in this estimate?
      logical, allocatable :: var_present(:)
      ! posn of model variable in estimate (0 if absent)
      integer(kind=our_int), allocatable :: estimate_posn(:) 
      integer(kind=our_int) :: nvar_present = 0 ! same as nvar
      integer(kind=our_int) :: nvar_absent = 0 ! ncol_input_data - nvar
      ! objects pertaining to the table that cross-classifies
      ! observations by base levels of vars in this estimate 
      integer(kind=our_int) :: ncells = 0
      integer(kind=our_int), allocatable :: cumprod(:)  ! length nvar
      real(kind=our_dble), allocatable :: prob(:) ! length ncells
      logical, allocatable :: str_zero(:) ! length ncells
      ! this one is allocated if model_type = "log-linear":
      real(kind=our_dble), allocatable :: SEs(:) ! length ncells
      real(kind=our_dble), allocatable :: xpi(:,:) ! ncells x p
      ! objects for cycling around the cells of the table corresponding
      ! to variables  being conditioned on (fixed), and those
      ! that are not (random); these are counters
      integer(kind=our_dble), allocatable :: var(:) ! length nvar
      integer(kind=our_dble) :: cell = 0 ! cell number
      integer(kind=our_dble) :: cell_fixed_part = 0 
      integer(kind=our_dble) :: cell_random_part = 0
      logical :: begin_cycle_fixed = .false.
      logical :: cycle_done_fixed = .false.
      logical :: begin_cycle_random = .false.
      logical :: cycle_done_random = .false.
      ! integer vector that shows, for each cell in the complete-data
      ! table, the corresponding cell of the estimate table
      integer(kind=our_int), allocatable :: cells(:)
      ! storage in main workspace
      integer(kind=our_int) :: first_in_estimates = 0
      integer(kind=our_int) :: last_in_estimates = 0
   end type workspace_type_cvam_estimate
   !##################################################################
   type :: workspace_type_cvam_predict
      ! stores information pertaining to predicting one or more
      ! variables
      ! the structure of this workspace is very similar to that of
      ! a submodel, except that there is no such thing as a fixed
      ! variable
      integer(kind=our_int) :: nvar = 0 ! number of vars to predict
      !### these arrays have length nvar
      ! number of levels this variable (length=nvar):
      integer(kind=our_int), allocatable :: nlev(:) 
      ! which position (column number) in the full model (input data)
      integer(kind=our_int), allocatable :: model_posn(:) 
      ! variable treated as fixed in this predict? (length=nvar):
      logical, allocatable :: fixed(:)  ! will all be .false.
      !### these arrays have length ncol_input_data
      ! is a model variable present in this estimate?
      logical, allocatable :: var_present(:)
      ! posn of model variable in predict (0 if absent)
      integer(kind=our_int), allocatable :: predict_posn(:) 
      integer(kind=our_int) :: nvar_present = 0 ! same as nvar
      integer(kind=our_int) :: nvar_absent = 0 ! ncol_input_data - nvar
      ! objects pertaining to the table that cross-classifies
      ! observations by base levels of vars in this estimate 
      integer(kind=our_int) :: ncells = 0
      integer(kind=our_int), allocatable :: cumprod(:)  ! length nvar
      real(kind=our_dble), allocatable :: prob(:) ! length ncells
      logical, allocatable :: str_zero(:) ! length ncells
      ! objects for cycling around the cells of the table corresponding
      ! to variables  being conditioned on (fixed), and those
      ! that are not (random); these are counters
      integer(kind=our_dble), allocatable :: var(:) ! length nvar
      integer(kind=our_dble) :: cell = 0 ! cell number
      integer(kind=our_dble) :: cell_fixed_part = 0 
      integer(kind=our_dble) :: cell_random_part = 0
      logical :: begin_cycle_fixed = .false.
      logical :: cycle_done_fixed = .false.
      logical :: begin_cycle_random = .false.
      logical :: cycle_done_random = .false.
      ! integer vector that shows, for each cell in the complete-data
      ! table, the corresponding cell of the estimate table
      integer(kind=our_int), allocatable :: cells(:)
   end type workspace_type_cvam_predict
   !##################################################################
   type :: workspace_type_cvam
      ! stores information for fitting a cvam model
      !  input_data(:,:) = integer( nrow_input_data, ncol_input_data)
      !     = matrix of coarsened data patterns
      !  input_data_freq = real( nrow_input_data )
      !     = vector of frequencies associated with input_data
      !  n_base_levels = integer( ncol_input_data )
      !  n_coarse_levels = integer( ncol_input_data )
      !  n_levels = integer( ncol_input_data )
      !  fixed_in_all = logical( ncol_input_data )
      !  mapping = ragged 3d integer array to hold the mappings from
      !     the levels in input_data to base levels
      private
      character(len=method_str_len) :: method = ""
      character(len=model_type_str_len) :: model_type = ""
      integer(kind=our_int) :: nrow_input_data = 0
      integer(kind=our_int) :: ncol_input_data = 0
      integer(kind=our_int), allocatable :: input_data(:,:)
      integer(kind=our_int), allocatable :: input_data_freq_int(:)
      real(kind=our_dble), allocatable :: input_data_freq(:)
      integer(kind=our_int), allocatable :: n_base_levels(:)
      integer(kind=our_int), allocatable :: n_coarse_levels(:)
      integer(kind=our_int), allocatable :: n_levels(:)
      logical, allocatable :: fixed(:)
      logical, allocatable :: fixed_tmp(:)
      type(workspace_type_int_array_3d) :: mapping
      ! objects pertaining to the table that cross-classifies
      ! observations by base levels of all vars (complete-data table)
      integer(kind=our_int) :: nvar = 0 ! same as ncol_input_data
      integer(kind=our_int) :: ncells = 0
      integer(kind=our_int), allocatable :: nlev(:)  !  same as n_base_levels
      integer(kind=our_int), allocatable :: cumprod(:)  ! length nvar
      integer(kind=our_int), allocatable :: mvcode(:)  ! length nvar
      ! objects used when cycling through all the cells of the 
      ! complete-data table that contribute to a given row
      ! of input_data
      logical :: begin_cycle = .false. ! set to .true. to start a new cycle
      logical :: cycle_done = .false. ! becomes .true. at end of the cycle
      integer(kind=our_int), allocatable :: var(:) ! length nvar,
      !   stores the base levels of the current cell
      integer(kind=our_int), allocatable :: var_index(:) ! length nvar,
      !   stores the index of each base_level in var in its mapping set
      logical, allocatable :: var_max(:) ! length nvar
      integer(kind=our_int) :: cell = 0 ! the cell of table
      ! information specific to the saturated or log-linear model
      integer(kind=our_int) :: n = 0 ! same as ncells, if log-linear
      integer(kind=our_int) :: p = 0 ! columns of model_matrix
      real(kind=our_dble), allocatable :: model_matrix(:,:) ! (n,p)
      real(kind=our_dble), allocatable :: offset(:) ! length n
      logical, allocatable :: str_zero(:) ! length ncells
      ! objects for cycling around the cells of the table corresponding
      ! to variables whose margins are fixed by the model, and those
      ! whose margins are random (not fixed)
      integer(kind=our_dble) :: cell_fixed_part = 0 
      integer(kind=our_dble) :: cell_random_part = 0
      logical :: begin_cycle_fixed = .false.
      logical :: cycle_done_fixed = .false.
      logical :: begin_cycle_random = .false.
      logical :: cycle_done_random = .false.
      ! prior stuff
      real(kind=our_dble) :: flatten = 0.D0
      real(kind=our_dble) :: ridge = 0.D0
      integer(kind=our_int) :: nrow_prior_data = 0
      integer(kind=our_int), allocatable :: prior_data(:,:)
      real(kind=our_dble), allocatable :: prior_data_freq(:)
      integer(kind=our_int), allocatable :: prior_data_freq_int(:)
      ! control parameters
      real(kind=our_dble) :: start_val_jitter = 0.D0
      real(kind=our_dble) :: crit_EM = 0.D0
      real(kind=our_dble) :: crit_NR = 0.D0
      real(kind=our_dble) :: crit_boundary = 0.D0
      integer(kind=our_int) :: iter_max_EM = 0
      integer(kind=our_int) :: iter_max_NR = 0
      logical :: start_val_use = .false.
      character(len=method_str_len) :: start_val_default = ""
      logical :: exclude_all_na = .false.
      logical :: omit_data = .false.
      integer(kind=our_int) :: iter_mcmc_nominal = 0
      integer(kind=our_int) :: iter_mcmc = 0
      integer(kind=our_int) :: burn_mcmc = 0
      integer(kind=our_int) :: thin_mcmc = 0
      integer(kind=our_int) :: impute_every = 0
      logical :: save_prob_series = .false.
      character(len=method_str_len) :: type_mcmc = ""
      integer(kind=our_int) :: stuck_limit = 0
      integer(kind=our_int) :: iter_approx_bayes = 0
      logical :: impute_approx_bayes = .false.
      real(kind=our_dble) :: df_da = 0.D0
      real(kind=our_dble) :: step_size_da = 0.D0
      real(kind=our_dble) :: scale_fac_da = 0.D0
      real(kind=our_dble) :: df_rwm = 0.D0
      real(kind=our_dble) :: scale_fac_rwm = 0.D0
      ! data use objects
      logical, allocatable :: input_data_use(:)
      logical, allocatable :: input_data_use_tmp(:)
      logical, allocatable :: prior_data_use(:)
      real(kind=our_dble) :: total_freq_use_prior = 0.D0
      integer(kind=our_int) :: total_freq_use_data_int = 0
      ! objects for model fitting
      real(kind=our_dble), allocatable :: freq(:) ! length ncells
      real(kind=our_dble), allocatable :: freq_mean(:) ! length ncells
      integer(kind=our_int), allocatable :: freq_int(:) ! length ncells
      real(kind=our_dble), allocatable :: prob(:) ! length ncells
      real(kind=our_dble), allocatable :: prob_new(:) ! length ncells
      ! these are never allocated if model_type is "saturated"
      real(kind=our_dble), allocatable :: beta(:) ! length p
      real(kind=our_dble), allocatable :: beta_new(:) ! length p
      real(kind=our_dble), allocatable :: beta_null(:) ! length p
      real(kind=our_dble), allocatable :: beta_hat(:) ! length p
      real(kind=our_dble), allocatable :: mu(:) ! length n
      real(kind=our_dble), allocatable :: mu_old(:) ! length n
      real(kind=our_dble), allocatable :: logmu(:) ! length n
      integer(kind=our_int) :: ncells_nonzero = 0
      integer(kind=our_int) :: degrees_of_freedom = 0
      ! workspaces that are never allocated if model_type is "saturated"
      real(kind=our_dble), allocatable :: wknA(:)
      real(kind=our_dble), allocatable :: wkpA(:)
      real(kind=our_dble), allocatable :: wkpB(:)
      real(kind=our_dble), allocatable :: beta_tmp(:)
      real(kind=our_dble), allocatable :: mu_tmp(:)
      real(kind=our_dble), allocatable :: logmu_tmp(:)
      integer(kind=our_int), allocatable :: freq_int_tmp(:)
      real(kind=our_dble), allocatable :: freq_tmp(:)
      real(kind=our_dble), allocatable :: wkppA(:,:)
      real(kind=our_dble), allocatable :: wkppB(:,:)
      real(kind=our_dble), allocatable :: lambda(:)
      ! objects for storing loglikelihood etc
      real(kind=our_dble) :: logprior = 0.D0
      real(kind=our_dble) :: loglik = 0.D0
      real(kind=our_dble) :: logP = 0.D0
      real(kind=our_dble) :: start_logP = 0.D0
      real(kind=our_dble) :: logP_mstep = 0.D0
      real(kind=our_dble), allocatable :: loglik_vec(:)
      real(kind=our_dble), allocatable :: logP_vec(:)
      real(kind=our_dble), allocatable :: score_mstep(:)
      real(kind=our_dble), allocatable :: hess_mstep(:,:)
      integer(kind=our_int) :: iter_mstep = 0
      logical :: converged_mstep = .false.
      real(kind=our_dble) :: max_diff_mstep = 0.D0
      integer(kind=our_int) :: iter = 0
      logical :: converged = .false.
      real(kind=our_dble) :: max_diff = 0.D0
      ! objects for computing score and Hessian of observed logP
      ! under a log-linear model, these are never allocated if
      ! model_type is saturated
      real(kind=our_dble), allocatable :: score(:)
      real(kind=our_dble), allocatable :: hess(:,:)
      real(kind=our_dble), allocatable :: Mrow(:)    ! length n
      real(kind=our_dble), allocatable :: fhat(:)    ! length n
      real(kind=our_dble), allocatable :: bigFhat(:) ! length n
      real(kind=our_dble), allocatable :: Mx(:,:) ! (n,p)
      real(kind=our_dble), allocatable :: cfac_vhat(:,:)
      real(kind=our_dble), allocatable :: vhat_beta(:,:)
      real(kind=our_dble), allocatable :: vhat_beta_rwm(:,:)
      real(kind=our_dble), allocatable :: dvec(:)
      logical :: vhat_ok = .false.
      ! objects used when cycling through all the cells of the 
      ! complete-data table that contribute to a given row
      ! of input_data
      logical :: begin_cycle_2 = .false. ! set to .true. to start a new cycle
      logical :: cycle_done_2 = .false. ! becomes .true. at end of the cycle
      integer(kind=our_int), allocatable :: var_2(:) ! length nvar,
      !   stores the base levels of the current cell
      integer(kind=our_int), allocatable :: var_index_2(:) ! length nvar,
      !   stores the index of each base_level in var in its mapping set
      logical, allocatable :: var_max_2(:) ! length nvar
      integer(kind=our_int) :: cell_2 = 0 ! the cell of table
      ! objects for processing and storing estimates
      integer(kind=our_int) :: n_estimates = 0
      type(workspace_type_cvam_estimate), allocatable :: estimates(:)
      integer(kind=our_int) :: n_packed_estimates = 0
      real(kind=our_dble), allocatable :: packed_estimates(:)
      real(kind=our_dble), allocatable :: packed_estimates_mean(:)
      ! this one is allocated if model_type = "log-linear"
      real(kind=our_dble), allocatable :: packed_SEs(:)
      type(workspace_type_cvam_predict) :: predict
      ! objects for storing the results from mcmc
      ! Note: to save space, we don't create these inside the
      ! workspace; that could change in the future
      integer(kind=our_int) :: series_length = 0
      integer(kind=our_int) :: n_impute = 0
      integer(kind=our_int) :: series_length_prob = 0
      integer(kind=our_int) :: beta_accept_count = 0
      integer(kind=our_int) :: beta_current_reject_run = 0
      real(kind=our_dble) :: beta_accept_rate = 0.D0
      real(kind=our_dble), allocatable :: mh_ratios_beta(:)
      logical, allocatable :: mh_accept_beta(:)
      real(kind=our_dble), allocatable :: beta_can(:)  ! candidate
      real(kind=our_dble), allocatable :: beta_center(:)
      real(kind=our_dble), allocatable :: beta_scale(:,:)
      real(kind=our_dble), allocatable :: beta_scale_inv(:,:)
      real(kind=our_dble), allocatable :: beta_scale_inv_sqrt(:,:)
      real(kind=our_dble), allocatable :: beta_scale_sqrt(:,:)
      real(kind=our_dble), allocatable :: beta_mean(:)
      real(kind=our_dble), allocatable :: beta_cov_mat(:,:)
      real(kind=our_dble), allocatable :: prob_mean(:)
      integer(kind=our_int) :: iter_past_burn_in = 0
      integer(kind=our_int) :: store_count = 0
      integer(kind=our_int) :: imp_count = 0
      logical :: store_this_iter = .false.
      logical :: imp_this_iter = .false.
   end type workspace_type_cvam
   !###################################################################
   contains
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_int_array_1d( &
        work, err ) result( answer )
      implicit none
      ! args
      type(workspace_type_int_array_1d), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_int_array_1d"
      integer(kind=our_int) :: status
      ! begin
      answer = RETURN_FAIL
      if( allocated( work%vec ) ) then
         deallocate( work%vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function nullify_workspace_type_int_array_1d
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_int_array_2d( &
        work, err ) result( answer )
      implicit none
      ! args
      type(workspace_type_int_array_2d), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_int_array_2d"
      integer(kind=our_int) :: i, status
      ! begin
      answer = RETURN_FAIL
      if( allocated( work%vec ) ) then
         do i = 1, size(work%vec)
            if( nullify_workspace_type_int_array_1d( work%vec(i), err ) &
                 == RETURN_FAIL ) goto 800
         end do
         deallocate( work%vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function nullify_workspace_type_int_array_2d
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_int_array_3d( &
        work, err ) result( answer )
      implicit none
      ! args
      type(workspace_type_int_array_3d), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_int_array_3d"
      integer(kind=our_int) :: i, status
      ! begin
      answer = RETURN_FAIL
      if( allocated( work%vec ) ) then
         do i = 1, size(work%vec)
            if( nullify_workspace_type_int_array_2d( work%vec(i), err ) &
                 == RETURN_FAIL ) goto 800
         end do
         deallocate( work%vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function nullify_workspace_type_int_array_3d
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_cvam_estimate( &
        work, err ) result( answer )
      implicit none
      type(workspace_type_cvam_estimate), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_cvam_estimate"
      integer(kind=our_int) :: status
      ! begin
      answer = RETURN_FAIL
      work%estimate_type_int = 0
      work%estimate_type = ""
      work%estimate_number = 0
      work%nvar = 0
      if( allocated( work%nlev ) ) then
         deallocate( work%nlev, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%model_posn ) ) then
         deallocate( work%model_posn, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fixed ) ) then
         deallocate( work%fixed, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_present ) ) then
         deallocate( work%var_present, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%estimate_posn ) ) then
         deallocate( work%estimate_posn, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%nvar_present = 0
      work%nvar_absent = 0
      work%ncells = 0
      if( allocated( work%cumprod ) ) then
         deallocate( work%cumprod, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prob ) ) then
         deallocate( work%prob, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%str_zero ) ) then
         deallocate( work%str_zero, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%SEs ) ) then
         deallocate( work%SEs, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%xpi ) ) then
         deallocate( work%xpi, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var ) ) then
         deallocate( work%var, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%cell = 0
      work%cell_fixed_part = 0
      work%cell_random_part = 0
      work%begin_cycle_fixed = .false.
      work%cycle_done_fixed = .false.
      work%begin_cycle_random = .false.
      work%cycle_done_random = .false.
      if( allocated( work%cells ) ) then
         deallocate( work%cells, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%first_in_estimates = 0
      work%last_in_estimates = 0
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function nullify_workspace_type_cvam_estimate
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_cvam_predict( &
        work, err ) result( answer )
      implicit none
      type(workspace_type_cvam_predict), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_cvam_predict"
      integer(kind=our_int) :: status
      ! begin
      work%nvar = 0
      if( allocated( work%nlev ) ) then
         deallocate( work%nlev, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%model_posn ) ) then
         deallocate( work%model_posn, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fixed ) ) then
         deallocate( work%fixed, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_present ) ) then
         deallocate( work%var_present, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%predict_posn ) ) then
         deallocate( work%predict_posn, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%nvar_present = 0
      work%nvar_absent = 0
      work%ncells = 0
      if( allocated( work%cumprod ) ) then
         deallocate( work%cumprod, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prob ) ) then
         deallocate( work%prob, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%str_zero ) ) then
         deallocate( work%str_zero, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var ) ) then
         deallocate( work%var, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%cell = 0
      work%cell_fixed_part = 0
      work%cell_random_part = 0
      work%begin_cycle_fixed = .false.
      work%cycle_done_fixed = .false.
      work%begin_cycle_random = .false.
      work%cycle_done_random = .false.
      if( allocated( work%cells ) ) then
         deallocate( work%cells, stat=status )
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function nullify_workspace_type_cvam_predict
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_cvam( &
        work, err ) result( answer )
      implicit none
      ! args
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_cvam"
      integer(kind=our_int) :: i, status
      ! begin
      answer = RETURN_FAIL
      work%method = ""
      work%model_type = ""
      work%nrow_input_data = 0
      work%ncol_input_data = 0
      if( allocated( work%input_data ) ) then
         deallocate( work%input_data, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%input_data_freq_int ) ) then
         deallocate( work%input_data_freq_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%input_data_freq ) ) then
         deallocate( work%input_data_freq, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%n_base_levels ) ) then
         deallocate( work%n_base_levels, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%n_coarse_levels ) ) then
         deallocate( work%n_coarse_levels, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%n_levels ) ) then
         deallocate( work%n_levels, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fixed ) ) then
         deallocate( work%fixed, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fixed_tmp ) ) then
         deallocate( work%fixed_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( nullify_workspace_type_int_array_3d( work%mapping, err ) &
           == RETURN_FAIL ) goto 800
      work%nvar = 0
      work%ncells = 0
      if( allocated( work%nlev ) ) then
         deallocate( work%nlev, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cumprod ) ) then
         deallocate( work%cumprod, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mvcode ) ) then
         deallocate( work%mvcode, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%begin_cycle = .false.
      work%cycle_done = .false.
      if( allocated( work%var ) ) then
         deallocate( work%var, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_index ) ) then
         deallocate( work%var_index, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_max ) ) then
         deallocate( work%var_max, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%cell = 0
      work%n = 0
      work%p = 0
      if( allocated( work%model_matrix ) ) then
         deallocate( work%model_matrix, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%offset ) ) then
         deallocate( work%offset, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%str_zero ) ) then
         deallocate( work%str_zero, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%cell_fixed_part = 0 
      work%cell_random_part = 0
      work%begin_cycle_fixed = .false.
      work%cycle_done_fixed = .false.
      work%begin_cycle_random = .false.
      work%cycle_done_random = .false.
      work%flatten = 0.D0
      work%ridge = 0.D0
      work%nrow_prior_data = 0
      if( allocated( work%prior_data ) ) then
         deallocate( work%prior_data, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prior_data_freq ) ) then
         deallocate( work%prior_data_freq, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prior_data_freq_int ) ) then
         deallocate( work%prior_data_freq_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%start_val_jitter = 0.D0
      work%crit_EM = 0.D0
      work%crit_NR = 0.D0
      work%crit_boundary = 0.D0
      work%iter_max_EM = 0
      work%iter_max_NR = 0
      work%start_val_use = .false.
      work%start_val_default = ""
      work%exclude_all_na = .false.
      work%omit_data = .false.
      work%iter_mcmc_nominal = 0
      work%iter_mcmc = 0
      work%burn_mcmc = 0
      work%thin_mcmc = 0
      work%impute_every = 0
      work%save_prob_series = .false.
      work%type_mcmc = ""
      work%stuck_limit = 0
      work%iter_approx_bayes = 0
      work%impute_approx_bayes = .false.
      work%df_da = 0.D0
      work%step_size_da = 0.D0
      work%scale_fac_da = 0.D0
      work%df_rwm = 0.D0
      work%scale_fac_rwm = 0.D0
      if( allocated( work%input_data_use ) ) then
         deallocate( work%input_data_use, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%input_data_use_tmp ) ) then
         deallocate( work%input_data_use_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prior_data_use ) ) then
         deallocate( work%prior_data_use, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%total_freq_use_prior = 0.D0
      work%total_freq_use_data_int = 0
      if( allocated( work%freq ) ) then
         deallocate( work%freq, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%freq_mean ) ) then
         deallocate( work%freq_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%freq_int ) ) then
         deallocate( work%freq_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prob ) ) then
         deallocate( work%prob, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prob_new ) ) then
         deallocate( work%prob_new, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta ) ) then
         deallocate( work%beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_new ) ) then
         deallocate( work%beta_new, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_null ) ) then
         deallocate( work%beta_null, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_hat ) ) then
         deallocate( work%beta_hat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mu ) ) then
         deallocate( work%mu, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mu_old ) ) then
         deallocate( work%mu_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%logmu ) ) then
         deallocate( work%logmu, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%ncells_nonzero = 0
      work%degrees_of_freedom = 0
      if( allocated( work%wknA ) ) then
         deallocate( work%wknA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkpA ) ) then
         deallocate( work%wkpA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkpB ) ) then
         deallocate( work%wkpB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_tmp ) ) then
         deallocate( work%beta_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mu_tmp ) ) then
         deallocate( work%mu_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%logmu_tmp ) ) then
         deallocate( work%logmu_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%freq_int_tmp ) ) then
         deallocate( work%freq_int_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%freq_tmp ) ) then
         deallocate( work%freq_tmp, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkppA ) ) then
         deallocate( work%wkppA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkppB ) ) then
         deallocate( work%wkppB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%lambda ) ) then
         deallocate( work%lambda, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%logprior = 0.D0
      work%loglik = 0.D0
      work%logP = 0.D0
      work%start_logP = 0.D0
      work%logP_mstep = 0.D0
      if( allocated( work%loglik_vec ) ) then
         deallocate( work%loglik_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%logP_vec ) ) then
         deallocate( work%logP_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%score_mstep ) ) then
         deallocate( work%score_mstep, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%hess_mstep ) ) then
         deallocate( work%hess_mstep, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%iter_mstep = 0
      work%converged_mstep = .false.
      work%max_diff_mstep = 0.D0
      work%iter = 0
      work%converged = .false.
      work%max_diff = 0.D0
      if( allocated( work%score ) ) then
         deallocate( work%score, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%hess ) ) then
         deallocate( work%hess, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%Mrow ) ) then
         deallocate( work%Mrow, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fhat ) ) then
         deallocate( work%fhat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%bigFhat ) ) then
         deallocate( work%bigFhat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%Mx ) ) then
         deallocate( work%Mx, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cfac_vhat ) ) then
         deallocate( work%cfac_vhat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%vhat_beta ) ) then
         deallocate( work%vhat_beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%vhat_beta_rwm ) ) then
         deallocate( work%vhat_beta_rwm, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%dvec ) ) then
         deallocate( work%dvec, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%vhat_ok = .false.
      ! normal exit
      answer = RETURN_SUCCESS
      return
      work%begin_cycle_2 = .false.
      work%cycle_done_2 = .false.
      if( allocated( work%var_2 ) ) then
         deallocate( work%var_2, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_index_2 ) ) then
         deallocate( work%var_index_2, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%var_max_2 ) ) then
         deallocate( work%var_max_2, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%cell_2 = 0
      work%n_estimates = 0
      if( allocated( work%estimates ) ) then
         do i = 1, size( work%estimates )
            if( nullify_workspace_type_cvam_estimate( work%estimates(i), &
                 err ) == RETURN_FAIL ) goto 800
         end do
         deallocate( work%estimates, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%n_packed_estimates = 0
      if( allocated( work%packed_estimates ) ) then
         deallocate( work%packed_estimates, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%packed_estimates_mean ) ) then
         deallocate( work%packed_estimates_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%packed_SEs ) ) then
         deallocate( work%packed_SEs, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( nullify_workspace_type_cvam_predict( work%predict, &
           err ) == RETURN_FAIL ) goto 800
      work%series_length = 0
      work%n_impute = 0
      work%series_length_prob = 0
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      if( allocated( work%mh_ratios_beta ) ) then
         deallocate( work%mh_ratios_beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mh_accept_beta ) ) then
         deallocate( work%mh_accept_beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_can ) ) then
         deallocate( work%beta_can, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_center ) ) then
         deallocate( work%beta_center, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale ) ) then
         deallocate( work%beta_scale, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_inv ) ) then
         deallocate( work%beta_scale_inv, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_inv_sqrt ) ) then
         deallocate( work%beta_scale_inv_sqrt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_sqrt ) ) then
         deallocate( work%beta_scale_sqrt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_mean ) ) then
         deallocate( work%beta_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_cov_mat ) ) then
         deallocate( work%beta_cov_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prob_mean ) ) then
         deallocate( work%prob_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      work%store_this_iter = .false.
      work%imp_this_iter = .false.
      ! error trapa
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
   end function nullify_workspace_type_cvam
   !##################################################################
   integer(kind=our_int) function run_cvam_model( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        flatten, ridge, prior_data, prior_data_freq_int, &
        ! elements of ctrl_real
        start_val_jitter, crit_EM, crit_NR, crit_boundary, &
        ! elements of ctrl_int
        iter_max_EM, iter_max_NR, &
        start_val_use_int, start_val_default_int, exclude_all_na_int, &
        omit_data_int, &
        estimate_info, estimate_var_info, &
        ! mcmc info
        iter_mcmc, burn_mcmc, thin_mcmc, impute_every, &
        save_prob_series_int, type_mcmc_int, stuck_limit, &
        iter_approx_bayes, impute_approx_bayes_int, &
        df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
        ! inout arguments
        prob, beta, beta_hat, vhat_beta_rwm, &
        ! workspaces
        work, err, &
        ! outputs
        iter, converged, max_diff, loglik, logP, lambda, &
        freq, freq_mean, freq_int, &
        score, vhat_beta, prob_mean, beta_mean, beta_cov_mat, &
        total_freq_use_prior, total_freq_use_data_int, &
        degrees_of_freedom, &
        packed_estimates, packed_estimates_mean, packed_SEs, &
        beta_series, prob_series, logp_series, imputed_freq_int, &
        packed_estimates_series, &
        n_iter_actual, n_sample_actual, n_imp_actual, mh_accept_rate, &
        start_logP &
        ) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: flatten
      real(kind=our_dble), intent(in) :: ridge
      integer(kind=our_int), intent(in) :: prior_data(:,:)
      integer(kind=our_int), intent(in) :: prior_data_freq_int(:)
      real(kind=our_dble), intent(in) :: start_val_jitter
      real(kind=our_dble), intent(in) :: crit_EM
      real(kind=our_dble), intent(in) :: crit_NR
      real(kind=our_dble), intent(in) :: crit_boundary
      integer(kind=our_int), intent(in) :: iter_max_EM
      integer(kind=our_int), intent(in) :: iter_max_NR
      integer(kind=our_int), intent(in) :: start_val_use_int
      integer(kind=our_int), intent(in) :: start_val_default_int
      integer(kind=our_int), intent(in) :: exclude_all_na_int
      integer(kind=our_int), intent(in) :: omit_data_int
      integer(kind=our_int), intent(in) :: estimate_info(:,:)
      integer(kind=our_int), intent(in) :: estimate_var_info(:,:)
      integer(kind=our_int), intent(in) :: iter_mcmc
      integer(kind=our_int), intent(in) :: burn_mcmc
      integer(kind=our_int), intent(in) :: thin_mcmc
      integer(kind=our_int), intent(in) :: impute_every
      integer(kind=our_int), intent(in) :: save_prob_series_int
      integer(kind=our_int), intent(in) :: type_mcmc_int
      integer(kind=our_int), intent(in) :: stuck_limit
      integer(kind=our_int), intent(in) :: iter_approx_bayes
      integer(kind=our_int), intent(in) :: impute_approx_bayes_int
      real(kind=our_dble), intent(in) :: df_da
      real(kind=our_dble), intent(in) :: step_size_da
      real(kind=our_dble), intent(in) :: scale_fac_da
      real(kind=our_dble), intent(in) :: df_rwm
      real(kind=our_dble), intent(in) :: scale_fac_rwm
      ! declare inouts
      real(kind=our_dble), intent(inout) :: prob(:)
      real(kind=our_dble), intent(inout) :: beta(:)
      real(kind=our_dble), intent(inout) :: beta_hat(:)
      real(kind=our_dble), intent(inout) :: vhat_beta_rwm(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare outputs
      integer(kind=our_int), intent(out) :: iter
      logical, intent(out) :: converged
      real(kind=our_dble), intent(out) :: max_diff
      real(kind=our_dble), intent(out) :: loglik(:)
      real(kind=our_dble), intent(out) :: logP(:)
      real(kind=our_dble), intent(out) :: lambda(:)
      real(kind=our_dble), intent(out) :: freq(:)
      real(kind=our_dble), intent(out) :: freq_mean(:)
      integer(kind=our_int), intent(out) :: freq_int(:)
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: vhat_beta(:,:)
      real(kind=our_dble), intent(out) :: prob_mean(:)
      real(kind=our_dble), intent(out) :: beta_mean(:)
      real(kind=our_dble), intent(out) :: beta_cov_mat(:,:)
      real(kind=our_dble), intent(out) :: total_freq_use_prior
      integer(kind=our_int), intent(out) :: total_freq_use_data_int
      integer(kind=our_int), intent(out) :: degrees_of_freedom
      real(kind=our_dble), intent(out) :: packed_estimates(:)
      real(kind=our_dble), intent(out) :: packed_estimates_mean(:)
      real(kind=our_dble), intent(out) :: packed_SEs(:)
      real(kind=our_dble), intent(out) :: beta_series(:,:)
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      real(kind=our_dble), intent(out) :: mh_accept_rate
      real(kind=our_dble), intent(out) :: start_logP
      ! declare locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_cvam_model"
      ! begin
      answer = RETURN_FAIL
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( input_data, input_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( put_prior_into_workspace( flatten, ridge, prior_data, &
           prior_data_freq_int, work, err ) == RETURN_FAIL ) goto 800
      if( put_ctrl_into_workspace( start_val_jitter, crit_EM, crit_NR, &
           crit_boundary, &
           iter_max_EM, iter_max_NR, start_val_use_int, &
           start_val_default_int, exclude_all_na_int, &
           omit_data_int, &
           iter_mcmc, burn_mcmc, thin_mcmc, impute_every, &
           save_prob_series_int, type_mcmc_int, stuck_limit, &
           iter_approx_bayes, impute_approx_bayes_int, &
           df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_data_and_prior_use_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      total_freq_use_prior = work%total_freq_use_prior
      total_freq_use_data_int = work%total_freq_use_data_int
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      degrees_of_freedom = work%degrees_of_freedom
      if( structural_zero_check( work, err ) == RETURN_FAIL ) goto 800
      if( create_estimate_objects( estimate_info, &
           estimate_var_info, work, err ) == RETURN_FAIL ) goto 800
      if( create_mcmc_objects( beta_series, prob_series, &
           logp_series, imputed_freq_int, packed_estimates_series, &
           work, err ) == RETURN_FAIL ) goto 800
      if( work%model_type == "saturated" ) then
         if( start_val_saturated( prob, work, err ) &
              == RETURN_FAIL ) goto 800
         prob(:) = work%prob_new(:)
         work%prob(:) = work%prob_new(:)
      else if( work%model_type == "log-linear" ) then
         if( start_val_log_linear( beta, work, err ) &
              == RETURN_FAIL ) goto 800
         prob(:) = work%prob_new(:)
         beta(:) = work%beta_new(:)
         work%beta(:) = work%beta_new(:)
      else
         goto 100
      end if
      if( work%model_type == "saturated" ) then
         if( work%method == "EM" ) then
            if( run_em_saturated( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "MCMC" ) then
               if( run_da_saturated( prob_series, &
                    logp_series, imputed_freq_int, &
                    packed_estimates_series, &
                    n_iter_actual, n_sample_actual, n_imp_actual, &
                    work, err ) &
                    == RETURN_FAIL ) goto 800
         else
            goto 60
         end if
      else if( work%model_type == "log-linear" ) then
         if( work%method == "EM" ) then
            if( run_em_log_linear( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "MCMC" ) then
            if( work%type_mcmc == "DA" ) then
               if( run_da_log_linear( beta_series, prob_series, &
                    logp_series, imputed_freq_int, &
                    packed_estimates_series, &
                    n_iter_actual, n_sample_actual, n_imp_actual, &
                    work, err ) &
                    == RETURN_FAIL ) goto 800
            else if( work%type_mcmc == "RWM" ) then
               if( size( vhat_beta_rwm, 1 ) /= work%p ) goto 33
               if( size( vhat_beta_rwm, 2 ) /= work%p ) goto 33
               work%vhat_beta_rwm(:,:) = vhat_beta_rwm(:,:)
               if( run_rwm_log_linear( beta_series, prob_series, &
                    logp_series, imputed_freq_int, &
                    packed_estimates_series, &
                    n_iter_actual, n_sample_actual, n_imp_actual, &
                    work, err ) &
                    == RETURN_FAIL ) goto 800
            else
               goto 55
            end if
         else if( work%method == "approxBayes" ) then
            if( size( beta_hat ) /= work%p ) goto 34
            work%beta_hat(:) = beta_hat(:)
            if( size( vhat_beta_rwm, 1 ) /= work%p ) goto 33
            if( size( vhat_beta_rwm, 2 ) /= work%p ) goto 33
            work%vhat_beta_rwm(:,:) = vhat_beta_rwm(:,:)
            if( run_approx_bayes_log_linear( beta_series, prob_series, &
                 logp_series, imputed_freq_int, &
                 packed_estimates_series, &
                 n_iter_actual, n_sample_actual, n_imp_actual, &
                 work, err ) &
                 == RETURN_FAIL ) goto 800
         else
            goto 60
         end if
      else
         goto 100
      end if 
      !
      start_logP = work%start_logP
      !
      if( work%method == "EM" ) then
         prob(:) = work%prob_new(:)
         if( work%model_type == "log-linear" ) beta(:) = work%beta_new(:)
      else if( work%method == "MCMC" ) then
         prob(:) = work%prob(:)
         if( work%model_type == "log-linear" ) then
            beta(:) = work%beta(:)
            work%beta_new(:) = work%beta(:) ! to prepare for lambda
         end if
      else if( work%method == "approxBayes" ) then
         prob(:) = work%prob(:)
         beta(:) = work%beta(:)
         work%beta_new(:) = work%beta(:) ! to prepare for lambda
      else
         goto 60
      end if
      iter = work%iter
      converged = work%converged
      max_diff = work%max_diff
      if( size(loglik) < work%iter ) goto 10
      loglik( 1:work%iter ) = work%loglik_vec( 1:work%iter )
      if( size(logP) < work%iter ) goto 15
      logP( 1:work%iter ) = work%logP_vec( 1:work%iter )
      if( work%model_type == "log-linear" ) then
         if( size(lambda) /= work%p ) goto 20
         if( compute_lambda_from_beta( work%beta_new, work%lambda, &
              work, err ) == RETURN_FAIL ) goto 800
         lambda(:) = work%lambda(:)
      end if
      if( size(freq) /= work%ncells ) goto 25
      freq(:) = work%freq(:)
      if( size(freq_int) /= work%ncells ) goto 26
      freq_int(:) = work%freq_int(:)
      if( work%model_type == "log-linear" ) then
         if( size(score) /= work%p ) goto 30
         score(:) = work%score(:)
      end if
      ! return parameter and uncertainty estimates other than beta and prob
      if( work%method == "EM" ) then
         if( work%model_type == "log-linear" ) then
            if( size(vhat_beta, 1) /= work%p ) goto 35
            if( size(vhat_beta, 2) /= work%p ) goto 35
            vhat_beta(:,:) = work%vhat_beta(:,:)
         end if
      else  if( work%method == "MCMC" ) then
         if( size( prob_mean ) /= work%ncells ) goto 36
         prob_mean(:) = work%prob_mean(:)
         if( size( freq_mean ) /= work%ncells ) goto 37
         freq_mean(:) = work%freq_mean(:)
         if( work%model_type == "log-linear" ) then
            if( size(beta_mean) /= work%p ) goto 38
            if( size(beta_cov_mat, 1) /= work%p ) goto 39
            if( size(beta_cov_mat, 2) /= work%p ) goto 39
            beta_mean(:) = work%beta_mean(:)
            beta_cov_mat(:,:) = work%beta_cov_mat(:,:)
         end if
      else  if( work%method == "approxBayes" ) then
         if( size( prob_mean ) /= work%ncells ) goto 36
         prob_mean(:) = work%prob_mean(:)
         if( size( freq_mean ) /= work%ncells ) goto 37
         freq_mean(:) = work%freq_mean(:)
         if( size(beta_mean) /= work%p ) goto 38
         if( size(beta_cov_mat, 1) /= work%p ) goto 39
         if( size(beta_cov_mat, 2) /= work%p ) goto 39
         beta_mean(:) = work%beta_mean(:)
         beta_cov_mat(:,:) = work%beta_cov_mat(:,:)
      else
         goto 60
      end if
      !
      if( work%method == "MCMC" ) then
         if( work%model_type == "log-linear" ) then
            mh_accept_rate = work%beta_accept_rate
         else
            mh_accept_rate = 1.D0
         end if
      end if
      if( work%method == "approxBayes" ) mh_accept_rate = 1.D0
      !
      if( work%method == "EM" ) then
         if( size( packed_estimates ) > 0 ) then
            if( size( packed_estimates ) /= work%n_packed_estimates ) goto 40
            if( compute_estimates( work%prob_new, work, err ) &
                 == RETURN_FAIL ) goto 800
            packed_estimates(:) = work%packed_estimates(:)
         end if
         if( ( size( packed_SEs ) > 0 ) .and. &
              ( work%model_type == "log-linear" ) ) then
            if( size( packed_SEs ) /= work%n_packed_estimates ) goto 45
            if( compute_estimate_SEs( work%prob_new, work, err ) &
                 == RETURN_FAIL ) goto 800
            packed_SEs(:) = work%packed_SEs(:)
         end if
      else if( work%method == "MCMC" ) then
         if( work%n_estimates > 0 ) then
            if( size( packed_estimates_mean ) /= &
                 work%n_packed_estimates ) goto 40
            packed_estimates_mean(:) = work%packed_estimates_mean(:)
         end if
      else if( work%method == "approxBayes" ) then
         if( work%n_estimates > 0 ) then
            if( size( packed_estimates_mean ) /= &
                 work%n_packed_estimates ) goto 40
            packed_estimates_mean(:) = work%packed_estimates_mean(:)
         end if
      else
         goto 60
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
10    call err_handle(err, 1, &
            comment = "Argument loglik has insufficient size" )
      goto 800
15    call err_handle(err, 1, &
            comment = "Argument logP has insufficient size" )
      goto 800
20    call err_handle(err, 1, &
            comment = "Argument lambda has incorrect size" )
      goto 800
25    call err_handle(err, 1, &
            comment = "Argument freq has incorrect size" )
      goto 800
26    call err_handle(err, 1, &
            comment = "Argument freq_int has incorrect size" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument score has incorrect size" )
      goto 800
33    call err_handle(err, 1, &
            comment = "Argument vhat_beta_rwm has incorrect size" )
      goto 800
34    call err_handle(err, 1, &
            comment = "Argument beta_hat has incorrect size" )
      goto 800
35    call err_handle(err, 1, &
            comment = "Argument vhat_beta has incorrect size" )
      goto 800
36    call err_handle(err, 1, &
            comment = "Argument prob_mean has incorrect size" )
      goto 800
37    call err_handle(err, 1, &
            comment = "Argument freq_mean has incorrect size" )
      goto 800
38    call err_handle(err, 1, &
            comment = "Argument beta_mean has incorrect size" )
      goto 800
39    call err_handle(err, 1, &
            comment = "Argument beta_cov_mat has incorrect size" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument packed_estimates has incorrect size" )
      goto 800
45    call err_handle(err, 1, &
            comment = "Argument packed_SEs has incorrect size" )
      goto 800
55    call err_handle(err, 1, &
            comment = "Value of type_mcmc not recognized" )
      goto 800
60    call err_handle(err, 1, &
            comment = "Value of method not recognized" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Value of model_type not recognized" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_model
   !##################################################################
   integer(kind=our_int) function put_model_type_into_workspace( &
        model_type_int, work, err) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      character(len=*), parameter :: &
           subname = "put_model_type_into_workspace"
      ! begin
      answer = RETURN_FAIL
      if( model_type_int == 1 ) then
         work%model_type = "saturated"
      else if( model_type_int == 2 ) then
         work%model_type = "log-linear"
      else
         goto 100
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
100   call err_handle(err, 1, &
            comment = "Value of model_type not recognized" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function put_model_type_into_workspace
   !##################################################################
   integer(kind=our_int) function put_method_into_workspace( &
        method_int, work, err) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: method_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      character(len=*), parameter :: &
           subname = "put_method_into_workspace"
      ! begin
      answer = RETURN_FAIL
      if( method_int == 1 ) then
         work%method = "EM"
      else if( method_int == 2 ) then
         work%method = "MCMC"
      else if( method_int == 3 ) then
        if( work%model_type == "saturated" ) goto 50
        work%method = "approxBayes"
      else
         goto 100
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
50    call err_handle(err, 1, &
           comment = "Approximate Bayes cannot be run on a saturated model")
      goto 800
100   call err_handle(err, 1, &
           comment = "Value of method not recognized" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function put_method_into_workspace
   !##################################################################
   integer(kind=our_int) function put_data_into_workspace( &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        work, err ) &
        result(answer)
      ! run basic checks on input data and store in workspace, then
      ! create workspace objects
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: status, posn, i, j, k, kmax
      character(len=*), parameter :: &
           subname = "put_data_into_workspace"
      ! begin
      answer = RETURN_FAIL
      ! transfer dimensions and check; allocate and copy basic input
      ! data
      work%nrow_input_data = size(input_data, 1)
      work%ncol_input_data = size(input_data, 2)
      allocate( work%input_data( work%nrow_input_data, &
           work%ncol_input_data ), stat=status )
      if( status /= 0 ) goto 100
      work%input_data(:,:) = input_data(:,:)
      if( size(input_data_freq_int) /= work%nrow_input_data ) goto 150
      allocate( work%input_data_freq_int( work%nrow_input_data ), &
           work%input_data_freq( work%nrow_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%input_data_freq_int(:) = input_data_freq_int(:)
      do i = 1, work%nrow_input_data
         work%input_data_freq(i) = real( work%input_data_freq_int(i), &
              kind=our_dble)
      end do
      if( size(n_levels_matrix,1) /= work%ncol_input_data ) goto 160
      if( size(n_levels_matrix,2) /= 4 ) goto 160
      allocate( work%n_base_levels( work%ncol_input_data ), &
           work%n_coarse_levels( work%ncol_input_data ), &
           work%n_levels( work%ncol_input_data ), &
           work%fixed( work%ncol_input_data ), &
           work%fixed_tmp( work%ncol_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%n_base_levels(:) = n_levels_matrix(:,1)
      work%n_coarse_levels(:) = n_levels_matrix(:,2)
      work%n_levels(:) = n_levels_matrix(:,3)
      do i = 1, work%ncol_input_data
         if( n_levels_matrix(i,4) == 0 ) then
            work%fixed(i) = .false.
         else if( n_levels_matrix(i,4) == 1 ) then
            work%fixed(i) = .true.
         else
            goto 165
         end if
      end do
      ! create ragged 3d array
      allocate( work%mapping%vec( work%ncol_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      posn = 0
      do i = 1, work%ncol_input_data
         posn = posn + 1
         if( posn > size( packed_map ) ) goto 170 
         if( packed_map(posn) /= work%n_levels(i) ) goto 170
         allocate( work%mapping%vec(i)%vec( packed_map(posn) ), &
              stat=status )
         if( status /= 0 ) goto 180
         do j = 1, work%n_levels(i)
            posn = posn + 1
            if( posn > size( packed_map ) ) goto 170 
            kmax = packed_map(posn)
            allocate( work%mapping%vec(i)%vec(j)%vec( kmax ), &
                 stat=status )
            if( status /= 0 ) goto 180
            do k = 1, kmax
               posn = posn + 1
               if( posn > size( packed_map ) ) goto 170 
               work%mapping%vec(i)%vec(j)%vec(k) = packed_map(posn)
            end do
         end do
      end do
      ! additional objects in cvam workspace that pertain to the data and
      ! are the same for any model
      work%nvar = work%ncol_input_data
      allocate( work%nlev( work%nvar ), stat=status )
      if( status /= 0 ) goto 100
      work%nlev(:) = work%n_base_levels(:)
      allocate( work%cumprod( work%nvar ), stat=status )
      if( status /= 0 ) goto 100
      do j = 1, work%nvar
         if( j == 1 ) then
            work%cumprod(j) = 1
         else
            work%cumprod(j) = work%cumprod(j-1) * work%nlev(j-1)
         end if
      end do
      if( work%nvar == 0 ) then
         work%ncells = 1
      else
         work%ncells = work%cumprod( work%nvar ) * &
              work%nlev( work%nvar )
      end if
      allocate( work%mvcode( work%nvar ), stat=status )
      if( status /= 0 ) goto 100
      do i = 1, work%nvar
         work%mvcode(i) = size( work%mapping%vec(i)%vec )
      end do
      ! objects used in cycling through cells of table
      allocate( work%var( work%nvar ), work%var_index( work%nvar ), &
           work%var_max( work%nvar ), stat=status )
      if( status /= 0 ) goto 100
      work%var(:) = 0
      work%var_index(:) = 0
      work%var_max(:) = .false.
      ! objects used in cycling through cells of table
      allocate( work%var_2( work%nvar ), work%var_index_2( work%nvar ), &
           work%var_max_2( work%nvar ), stat=status )
      if( status /= 0 ) goto 100
      work%var_2(:) = 0
      work%var_index_2(:) = 0
      work%var_max_2(:) = .false.
      !###
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Sizes of input_data and input_data_freq_int mismatched" )
      goto 800
160   call err_handle(err, 1, &
            comment = "Size of n_levels_matrix not correct" )
      goto 800
165   call err_handle(err, 1, &
            comment = "Value in n_levels_matrix not correct" )
      goto 800
170   call err_handle(err, 1, &
            comment = "When creating mapping: array bounds exceeded" )
      goto 800
180   call err_handle(err, 1, &
            comment = "When creating mapping: allocation error " )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function put_data_into_workspace
   !##################################################################
   integer(kind=our_int) function put_model_into_workspace( &
        model_matrix, offset, str_zero_int, work, err) result(answer)
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, status
      character(len=*), parameter :: &
           subname = "put_model_into_workspace"
      ! begin
      answer = RETURN_FAIL
      ! model_matrix
      if( work%model_type == "saturated" ) then
         ! n = p = 0 and model_matrix is never allocated
         ! offset is never allocated
      else if( work%model_type == "log-linear" ) then
         if( size( model_matrix, 1 ) /= work%ncells ) goto 150
         work%n = work%ncells
         work%p = size( model_matrix, 2 ) 
         allocate( work%model_matrix( work%n, work%p ), stat=status )
         if( status /= 0 ) goto 100
         work%model_matrix(:,:) = model_matrix(:,:)
         if( size( offset ) /= work%ncells ) goto 160
         allocate( work%offset( work%n ), stat=status )
         if( status /= 0 ) goto 100
         work%offset(:) = offset(:)
      else
         goto 50
      end if
      ! structural zeros
      if( size( str_zero_int ) /= work%ncells ) goto 170
      allocate( work%str_zero( work%ncells ), stat=status )
      if( status /= 0 ) goto 100
      work%ncells_nonzero = 0
      do i = 1, work%ncells
         if( str_zero_int(i) == 1 ) then
            work%str_zero(i) = .true.
         else if( str_zero_int(i) == 0 ) then
            work%str_zero(i) = .false.
            work%ncells_nonzero = work%ncells_nonzero + 1
         else
            goto 200
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
50    call err_handle(err, 1, &
            comment = "Value of model_type not recognized" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Argument model_matrix has incorrect number of rows" )
      goto 800
160   call err_handle(err, 1, &
            comment = "Argument offset has incorrect size" )
      goto 800
170   call err_handle(err, 1, &
            comment = "Argument str_zero_int has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Element of str_zero_int is not 0 or 1" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function put_model_into_workspace
   !##################################################################
   integer(kind=our_int) function put_prior_into_workspace( &
        flatten, ridge, prior_data, prior_data_freq_int, &
        work, err ) result(answer)
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: flatten
      real(kind=our_dble), intent(in) :: ridge
      integer(kind=our_int), intent(in) :: prior_data(:,:)
      integer(kind=our_int), intent(in) :: prior_data_freq_int(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: status
      character(len=*), parameter :: &
           subname = "put_prior_into_workspace"
      ! begin
      answer = RETURN_FAIL
      if( flatten < 0.D0 ) goto 50
      work%flatten = flatten
      if( ridge < 0.D0 ) goto 60
      work%ridge = ridge
      if( ( work%ridge > 0.D0 ) .and. ( work%model_type == "saturated " ) ) &
           call err_handle(err, 1, &
           comment = "Prior ridge factor ignored; model is saturated" )
      work%nrow_prior_data = size(prior_data, 1)
      allocate( work%prior_data( work%nrow_prior_data, work%nvar ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%prior_data(:,:) = prior_data(:,:)
      allocate( work%prior_data_freq( work%nrow_prior_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      if( size( prior_data_freq_int ) /= work%nrow_prior_data ) goto 151
      allocate( work%prior_data_freq_int( work%nrow_prior_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%prior_data_freq_int(:) = prior_data_freq_int(:)
      work%prior_data_freq(:) = real( work%prior_data_freq_int(:), our_dble )
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
50    call err_handle(err, 1, &
            comment = "Negative value for flattening constant" )
      goto 800
60    call err_handle(err, 1, &
            comment = "Negative value for ridge factor" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
151   call err_handle(err, 1, &
            comment = "Sizes of prior_data and prior_data_freq_int mismatched" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function put_prior_into_workspace
   !##################################################################
   integer(kind=our_int) function put_ctrl_into_workspace( &
        start_val_jitter, crit_EM, crit_NR, crit_boundary, &
        iter_max_EM, iter_max_NR, start_val_use_int, &
        start_val_default_int, exclude_all_na_int, omit_data_int, &
        iter_mcmc, burn_mcmc, thin_mcmc, impute_every, &
        save_prob_series_int, type_mcmc_int, stuck_limit, &
        iter_approx_bayes, impute_approx_bayes_int, &
        df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
        work, err ) result(answer)
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: start_val_jitter
      real(kind=our_dble), intent(in) :: crit_EM
      real(kind=our_dble), intent(in) :: crit_NR
      real(kind=our_dble), intent(in) :: crit_boundary
      integer(kind=our_int), intent(in) :: iter_max_EM
      integer(kind=our_int), intent(in) :: iter_max_NR
      integer(kind=our_int), intent(in) :: start_val_use_int
      integer(kind=our_int), intent(in) :: start_val_default_int
      integer(kind=our_int), intent(in) :: exclude_all_na_int
      integer(kind=our_int), intent(in) :: omit_data_int
      integer(kind=our_int), intent(in) :: iter_mcmc
      integer(kind=our_int), intent(in) :: burn_mcmc
      integer(kind=our_int), intent(in) :: thin_mcmc
      integer(kind=our_int), intent(in) :: impute_every
      integer(kind=our_int), intent(in) :: save_prob_series_int
      integer(kind=our_int), intent(in) :: type_mcmc_int
      integer(kind=our_int), intent(in) :: stuck_limit
      integer(kind=our_int), intent(in) :: iter_approx_bayes
      integer(kind=our_int), intent(in) :: impute_approx_bayes_int
      real(kind=our_dble), intent(in) :: df_da
      real(kind=our_dble), intent(in) :: step_size_da
      real(kind=our_dble), intent(in) :: scale_fac_da
      real(kind=our_dble), intent(in) :: df_rwm
      real(kind=our_dble), intent(in) :: scale_fac_rwm
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      character(len=*), parameter :: &
           subname = "put_ctrl_into_workspace"
      ! begin
      answer = RETURN_FAIL
      if( start_val_jitter < 0.D0 ) goto 100
      work%start_val_jitter = start_val_jitter
      if( crit_EM <= 0.D0 ) goto 120
      work%crit_EM = crit_EM
      if( crit_NR <= 0.D0 ) goto 140
      work%crit_NR = crit_NR
      if( crit_boundary <= 0.D0 ) goto 150
      work%crit_boundary = crit_boundary
      if( iter_max_EM < 0 ) goto 160
      work%iter_max_EM = iter_max_EM
      if( iter_max_NR < 0 ) goto 180
      work%iter_max_NR = iter_max_NR
      if( start_val_use_int == 0 ) then
         work%start_val_use = .false.
      else if( start_val_use_int == 1 ) then
         work%start_val_use = .true.
      else
         goto 200 
      end if
      if( start_val_default_int == 1 ) then
         work%start_val_default = "center"
      else if( start_val_default_int == 2 ) then
         work%start_val_default = "uniform"
      else
         goto 220 
      end if
      if( exclude_all_na_int == 0 ) then
         work%exclude_all_na = .false.
      else if( exclude_all_na_int == 1 ) then
         work%exclude_all_na = .true.
      else
         goto 300
      end if
      if( omit_data_int == 0 ) then
         work%omit_data = .false.
      else if( omit_data_int == 1 ) then
         work%omit_data = .true.
      else
         goto 300
      end if
      !
      if( iter_mcmc < 0 ) goto 400
      work%iter_mcmc_nominal = iter_mcmc
      if( burn_mcmc < 0 ) goto 410
      work%burn_mcmc = burn_mcmc
      if( thin_mcmc < 0 ) goto 420
      work%thin_mcmc = thin_mcmc
      if( impute_every < 0 ) goto 430
      work%impute_every = impute_every
      work%iter_mcmc = floor( real(iter_mcmc, our_dble) / &
         real(thin_mcmc, our_dble) ) * thin_mcmc
      if( mod( iter_mcmc, thin_mcmc ) > 0 ) work%iter_mcmc = &
           work%iter_mcmc + work%thin_mcmc
      if( save_prob_series_int == 0 ) then
         work%save_prob_series = .false.
      else if( save_prob_series_int == 1 ) then
         work%save_prob_series = .true.
      else
         goto 440
      end if
      if( type_mcmc_int == 1 ) then
         work%type_mcmc = "DA"
      else if( type_mcmc_int == 2 ) then
         work%type_mcmc = "RWM"
      else
         goto 450
      end if
      if( stuck_limit <= 0 ) goto 455
      work%stuck_limit = stuck_limit
      if( iter_approx_bayes < 0 ) goto 456
      work%iter_approx_bayes = iter_approx_bayes
      if( impute_approx_bayes_int == 0 ) then
         work%impute_approx_bayes = .false.
      else if( impute_approx_bayes_int == 1 ) then
         work%impute_approx_bayes = .true.
      else
         goto 457
      end if
      if( df_da <= 0.D0 ) goto 460
      work%df_da = df_da
      work%step_size_da = step_size_da
      if( scale_fac_da <= 0.D0 ) goto 470
      work%scale_fac_da = scale_fac_da
      if( df_rwm <= 0.D0 ) goto 460
      work%df_rwm = df_rwm
      if( scale_fac_rwm <= 0.D0 ) goto 470
      work%scale_fac_rwm = scale_fac_rwm
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Negative value for start_val_jitter" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Non-positive value for crit_EM" )
      goto 800
140   call err_handle(err, 1, &
            comment = "Non-positive value for crit_NR" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Non-positive value for crit_boundary" )
      goto 800
160   call err_handle(err, 1, &
            comment = "Negative value for iter_max_EM" )
      goto 800
180   call err_handle(err, 1, &
            comment = "Negative value for iter_max_NR" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Value of start_val_use_int not recognized" )
      goto 800
220   call err_handle(err, 1, &
            comment = "Value of start_val_default_int not recognized" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Value of exclude_all_na_int not recognized" )
      goto 800
400   call err_handle(err, 1, &
            comment = "Negative value for iter_mcmc" )
      goto 800
410   call err_handle(err, 1, &
            comment = "Negative value for burn_mcmc" )
      goto 800
420   call err_handle(err, 1, &
            comment = "Negative value for thin_mcmc" )
      goto 800
430   call err_handle(err, 1, &
            comment = "Negative value for impute_every" )
      goto 800
440   call err_handle(err, 1, &
            comment = "Value for save_prob_series_int not recognized" )
      goto 800
450   call err_handle(err, 1, &
            comment = "Value for type_mcmc_int not recognized" )
      goto 800
455   call err_handle(err, 1, &
            comment = "Value for stuck_limit not positive" )
      goto 800
456   call err_handle(err, 1, &
            comment = "Negative value for iter_approx_bayes" )
      goto 800
457   call err_handle(err, 1, &
            comment = "Value for impute_approx_bayes_int not recognized" )
      goto 800
460   call err_handle(err, 1, &
            comment = "Non-positive value for mcmc proposal df" )
      goto 800
470   call err_handle(err, 1, &
            comment = "Non-positive value for mcmc proposal scale_fac" )
      goto 800

800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function put_ctrl_into_workspace
   !##################################################################
   integer(kind=our_int) function create_data_and_prior_use_objects( &
        work, err ) result(answer)
      implicit none
      ! declare inputs
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: status, isum, i, j
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "create_data_and_prior_use_objects"
      ! begin
      answer = RETURN_FAIL
      allocate( work%input_data_use( work%nrow_input_data ), &
           work%input_data_use_tmp( work%nrow_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%input_data_use(:) = .false.
      work%total_freq_use_data_int = 0
      if( work%omit_data ) then
         ! nothing to do
      else
         isum = 0
         if( work%exclude_all_na ) then
            do i = 1, work%nrow_input_data
               do j = 1, work%nvar
                  if( ( .not. work%fixed(j) ) .and. &
                       ( work%input_data(i,j) /= work%mvcode(j) ) ) then
                     work%input_data_use(i) = .true.
                     exit
                  end if
               end do
               if( work%input_data_use(i) ) isum = isum +  &
                    work%input_data_freq_int(i)
            end do
         else
            do i = 1, work%nrow_input_data
               work%input_data_use(i) = .true.
               isum = isum + work%input_data_freq_int(i)
            end do
         end if
         work%total_freq_use_data_int = isum
      end if
      !
      allocate( work%prior_data_use( work%nrow_prior_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%prior_data_use(:) = .false.
      sum = 0.D0
      if( work%exclude_all_na ) then
         do i = 1, work%nrow_prior_data
            do j = 1, work%nvar
               if( ( .not. work%fixed(j) ) .and. &
                    ( work%prior_data(i,j) /= work%mvcode(j) ) ) then
                  work%prior_data_use(i) = .true.
                  exit
               end if
            end do
            if( work%prior_data_use(i) ) sum = sum +  &
                 real( work%prior_data_freq_int(i), our_dble )
         end do
      else
         do i = 1, work%nrow_prior_data
            work%prior_data_use(i) = .true.
            sum = sum + real( work%prior_data_freq_int(i), our_dble )
         end do
      end if
      work%total_freq_use_prior = work%flatten + sum
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function create_data_and_prior_use_objects
   !##################################################################
   integer(kind=our_int) function create_model_fitting_objects( &
        work, err ) result(answer)
      ! creates workspace arrays for fitting the model
      implicit none
      ! declare inputs
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: status
      character(len=*), parameter :: &
           subname = "create_model_fitting_objects"
      ! begin
      answer = RETURN_FAIL
      allocate( work%freq( work%ncells ), work%freq_mean( work%ncells ), &
           work%freq_int( work%ncells ), work%prob( work%ncells ), &
           work%prob_new( work%ncells ), &
           work%freq_tmp( work%ncells ), &
           work%freq_int_tmp( work%ncells ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%freq(:) = 0.D0
      work%freq_mean(:) = 0.D0
      work%freq_int(:) = 0
      work%prob(:) = 0.D0
      work%prob_new(:) = 0.D0
      !
      if( work%model_type == "saturated" ) then
         ! n = p = 0 and beta, beta_new, mu, logmu are never allocated
         work%degrees_of_freedom = 0
      else if( work%model_type == "log-linear" ) then
         allocate( work%beta( work%p ), work%beta_new( work%p ), &
              work%beta_null( work%p ), work%beta_hat( work%p ), &
              work%mu( work%n ), work%mu_old( work%n ), &
              work%logmu( work%n ), stat=status )
         if( status /= 0 ) goto 100
         work%beta(:) = 0.D0
         work%beta_new(:) = 0.D0
         work%beta_null(:) = 0.D0
         work%beta_hat(:) = 0.D0
         work%mu(:) = 0.D0
         work%mu_old(:) = 0.D0
         work%logmu(:) = 0.D0
         work%degrees_of_freedom = work%ncells_nonzero - work%p
         allocate( work%wknA( work%n ), work%wkpA( work%p ), &
              work%wkpB( work%p ), work%beta_tmp( work%p ), &
              work%mu_tmp( work%n ), work%logmu_tmp( work%n ), &
              work%wkppA( work%p, work%p ), work%wkppB( work%p, work%p ), &
              work%lambda( work%p ), &
              work%score_mstep( work%p ), work%hess_mstep( work%p, work%p ), &
              work%score( work%p ), work%hess( work%p, work%p ), &
              work%Mrow( work%n ), work%fhat( work%n ), &
              work%bigFhat( work%n ), work%Mx( work%n, work%p ), &
              work%cfac_vhat( work%p, work%p ), &
              work%vhat_beta( work%p, work%p ), &
              work%vhat_beta_rwm( work%p, work%p ), &
              work%dvec( work%p ), &
              stat=status )
         if( status /= 0 ) goto 100
         work%wknA(:) = 0.D0
         work%wkpA(:) = 0.D0
         work%wkpB(:) = 0.D0
         work%beta_tmp(:) = 0.D0
         work%wkppA(:,:) = 0.D0
         work%wkppB(:,:) = 0.D0
         work%lambda(:) = 0.D0
         work%score_mstep(:) = 0.D0
         work%hess_mstep(:,:) = 0.D0
         work%score(:) = 0.D0
         work%hess(:,:) = 0.D0
         work%Mrow(:) = 0.D0
         work%fhat(:) = 0.D0
         work%bigFhat(:) = 0.D0
         work%Mx(:,:) = 0.D0
         work%cfac_vhat(:,:) = 0.D0
         work%vhat_beta(:,:) = 0.D0
         work%vhat_beta_rwm(:,:) = 0.D0
         work%dvec(:) = 0.D0
      else
         goto 50
      end if
      if( work%method == "EM" ) then
         allocate( work%loglik_vec( work%iter_max_EM ), &
              work%logP_vec( work%iter_max_EM ), stat=status )
         if( status /= 0 ) goto 100
         work%loglik_vec(:) = 0.D0
         work%logP_vec(:) = 0.D0
      else if( work%method == "MCMC" ) then
         allocate( work%loglik_vec( work%iter_mcmc + work%burn_mcmc ), &
              work%logP_vec( work%iter_mcmc + work%burn_mcmc ), &
              stat=status )
         if( status /= 0 ) goto 100
         work%loglik_vec(:) = 0.D0
         work%logP_vec(:) = 0.D0
      else if( work%method == "approxBayes" ) then
         allocate( work%loglik_vec( work%iter_approx_bayes ), &
              work%logP_vec( work%iter_approx_bayes ), &
              stat=status )
         if( status /= 0 ) goto 100
         work%loglik_vec(:) = 0.D0
         work%logP_vec(:) = 0.D0
      else
         goto 30
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
30    call err_handle(err, 1, &
            comment = "Value of method not recognized" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Value of model_type not recognized" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function create_model_fitting_objects
   !##################################################################
   integer(kind=our_int) function start_val_saturated( prob, &
        work, err ) result(answer)
      ! For model_type = "saturated"
      !   * If start_val_use is .true., begins with supplied values 
      !     in prob.  Then jitters them, and normalizes to sum to one.
      !   * If start_val_use is .false., creates starting
      !     values according to start_val_default. If "center",
      !     starts at the center of the parameter space 
      !     (equal across cells), then jitters them and normalizes.
      !     If "uniform", generates random probs from a uniform
      !     density over the simplex.
      ! puts the result into work%prob_new
      implicit none
      ! declare inouts
      real(kind=our_dble), intent(in) :: prob(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: rtmp, log_prob
      character(len=*), parameter :: &
           subname = "start_val_saturated"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "saturated" ) goto 10
      if( size(prob) /= work%ncells ) goto 20
      if( work%start_val_jitter < 0.D0 ) goto 90
      !
      if( work%start_val_use ) then
         work%prob(:) = prob(:)
         ! check for correct pattern of zeros and nonzeros
         do i = 1, work%ncells
            if( work%str_zero(i) .and. work%prob(i) /= 0.D0 ) goto 100
            if( .not. work%str_zero(i) .and. work%prob(i) <= 0.D0 ) goto 105
         end do
         if( work%start_val_jitter > 0.D0 ) then
            do i = 1, work%ncells
               if( work%str_zero(i) ) cycle
               log_prob = log( work%prob(i) )
               if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
               log_prob = log_prob + rtmp * work%start_val_jitter
               if( log_prob < log_tiny ) goto 110
               if( log_prob > log_huge ) goto 120
               work%prob(i) = exp( log_prob )
            end do
         end if
         if( normalize_prob( work%prob, work%prob_new, work, err ) &
              == RETURN_FAIL ) goto 800
      else
         if( work%start_val_default == "center" ) then
            ! fill prob with equal values
            do i = 1, work%ncells
               if( work%str_zero(i) ) then
                  work%prob(i) = 0.D0
               else
                  work%prob(i) = 1.D0
               end if
            end do
            if( work%start_val_jitter > 0.D0 ) then
               do i = 1, work%ncells
                  if( work%str_zero(i) ) cycle
                  log_prob = log( work%prob(i) )
                  if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
                  log_prob = log_prob + rtmp * work%start_val_jitter
                  if( log_prob < log_tiny ) goto 110
                  if( log_prob > log_huge ) goto 120
                  work%prob(i) = exp( log_prob )
               end do
            end if
            if( normalize_prob( work%prob, work%prob_new, &
                 work, err ) == RETURN_FAIL ) goto 800
         else if( work%start_val_default == "uniform" ) then
            do i = 1, work%ncells
               if( work%str_zero(i) ) cycle
               if( runif_R( rtmp, err ) == RETURN_FAIL ) goto 800
               work%prob(i) = rtmp
            end do
            if( normalize_prob( work%prob, work%prob_new, &
                 work, err ) == RETURN_FAIL ) goto 800
         else
            goto 50
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
10    call err_handle(err, 1, &
            comment = "Incorrect model type" )
      goto 800
20    call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Value of start_val_default not recognized" )
      goto 800
90    call err_handle(err, 1, &
           comment = "Negative value provided for start_val_jitter" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Non-zero starting value provided for structural zero" )
      call err_handle(err, 15, icell = i )
      goto 800
105   call err_handle(err, 1, &
            comment = "Starting value has non-positive probability" )
      call err_handle(err, 15, icell = i )
      goto 800
110   call err_handle(err, 1, &
            comment = "Underflow; jittered cell prob became zero" )
      call err_handle(err, 15, icell = i )
      goto 800
120   call err_handle(err, 1, &
            comment = "Overflow; jittered cell prob became too large" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function start_val_saturated
   !##################################################################
   integer(kind=our_int) function start_val_log_linear( beta, &
        work, err ) result(answer)
      ! For model_type "log-linear"
      !   * If start_val_use is .true., begins with supplied values 
      !     in beta, then jitters by adding random gaussian noise
      !     with standard deviation given by start_val_jitter
      !   * If start_val_use is .false., creates starting
      !     values according to start_val_default. If "center",
      !     starts at the center of the parameter space. The center
      !     means that prob is made proportional to offset. Then
      !     jitters them. If "uniform", starts with random probs that
      !     are uniformly distributed over the cells; performs an
      !     e-step, then fits the log-linear model to those frequencies.
      !   * Result is stored in beta_new, with corresponding probs
      !     in prob_new
      implicit none
      ! declare inouts
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: rtmp, sum
      character(len=*), parameter :: &
           subname = "start_val_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 10
      if( size(beta) /= work%p ) goto 40
      if( work%start_val_jitter < 0.D0 ) goto 90
      !
      ! regress vector of constants log(total N) - log(sum(exp(offset)))
      ! on model_matrix using ols, putting result in beta_null
      sum = 0.D0
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         if( work%offset(i) < log_tiny ) goto 110
         if( work%offset(i) > log_huge ) goto 120
         sum = sum + exp( work%offset(i) )
      end do
      rtmp = real(work%total_freq_use_data_int, our_dble) + &
           work%total_freq_use_prior
      if( sum <= 0.D0 ) goto 200
      if( rtmp <= 0.D0 ) goto 200
      rtmp = log(rtmp) - log(sum)
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         work%wknA(i) = rtmp
      end do
      if( compute_ls_fit( work%wknA, work%beta_null, work, err ) &
           == RETURN_FAIL ) goto 800
      !
      if( work%start_val_use ) then
         if( work%start_val_jitter > 0.D0 ) then
            work%beta(:) = beta(:)
            do i = 1, work%p
               if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
               work%beta_new(i) = work%beta(i) + &
                    rtmp * work%start_val_jitter
            end do
         else
            work%beta_new(:) = beta(:)
         end if
         if( compute_mu_from_beta( work%beta_new, work, err ) &
              == RETURN_FAIL ) goto 800
         if( normalize_prob( work%mu, work%prob_new, work, err ) &
              == RETURN_FAIL ) goto 800
      else
         if( work%start_val_default == "center" ) then
            if( work%start_val_jitter > 0.D0 ) then
               do i = 1, work%p
                  if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
                  work%beta_new(i) = work%beta_null(i) + &
                       rtmp * work%start_val_jitter
               end do
            else
               work%beta_new(:) = work%beta_null(:)   
            end if
            if( compute_mu_from_beta( work%beta_new, work, err ) &
                 == RETURN_FAIL ) goto 800
            if( normalize_prob( work%mu, work%prob_new, work, err ) &
                 == RETURN_FAIL ) goto 800
         else if( work%start_val_default == "uniform" ) then
            do i = 1, work%ncells
               if( work%str_zero(i) ) cycle
               if( runif_R( rtmp, err ) == RETURN_FAIL ) goto 800
               work%prob_new(i) = rtmp
            end do
            if( normalize_prob( work%prob_new, work%prob, &
                 work, err ) == RETURN_FAIL ) goto 800
            if( run_estep( work, err, &
                 use_flatten = .true., use_prior_data = .true., &
                 use_input_data = .true.) == RETURN_FAIL ) goto 800
            if( run_mstep_log_linear( work, err ) == RETURN_FAIL ) goto 800
            if( compute_mu_from_beta( work%beta_new, work, err ) &
                 == RETURN_FAIL ) goto 800
            if( normalize_prob( work%mu, work%prob_new, work, err ) &
                 == RETURN_FAIL ) goto 800
         else
            goto 50
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
10    call err_handle(err, 1, &
            comment = "Incorrect model type" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Value of start_val_default not recognized" )
      goto 800
90    call err_handle(err, 1, &
           comment = "Negative value provided for start_val_jitter" )
      goto 800

110   call err_handle(err, 1, &
            comment = "Underflow; offset too small" )
      call err_handle(err, 15, icell = i )
      goto 800
120   call err_handle(err, 1, &
            comment = "Overflow; offset too large" )
      call err_handle(err, 15, icell = i )
      goto 800

200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function start_val_log_linear
   !##################################################################
   integer(kind=our_int) function compute_mu_from_beta( beta, &
        work, err ) result(answer)
      ! only call this function if model_type = "log-linear"
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j
      real(kind=our_dble) :: sum
      character(len=*), parameter :: subname = "compute_mu_from_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( size(beta) /= work%p ) goto 50
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         sum = work%offset(i)
         do j = 1, work%p
            sum = sum + work%model_matrix(i,j) * beta(j)
         end do
         work%logmu(i) = sum
         if( sum < log_tiny ) goto 110
         if( sum > log_huge ) goto 120
         work%mu(i) = exp( sum )
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Underflow; cell mean became zero" )
      call err_handle(err, 15, icell = i )
      goto 800
120   call err_handle(err, 1, &
            comment = "Overflow; cell mean became too large" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_mu_from_beta
   !##################################################################
   integer(kind=our_int) function normalize_prob( mu, prob, &
        work, err ) result(answer)
      ! scale the cell means in mu to become cell probs, conditioning
      ! on the variables that are fixed in the model
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: mu(:)
      ! declare outputs
      real(kind=our_dble), intent(out) :: prob(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      logical :: non_zero_cells_present
      real(kind=our_dble) :: sum
      character(len=*), parameter :: subname = "normalize_prob"
      ! begin
      answer = RETURN_FAIL
      !
      if( size(mu) /= work%ncells ) goto 50
      if( size(prob) /= work%ncells ) goto 100
      !
      work%begin_cycle_fixed = .true.
      work%cycle_done_fixed = .false.
      do
         if( advance_cell_fixed_part( work, err ) &
              == RETURN_FAIL ) goto 800
         non_zero_cells_present = .false.
         sum = 0.D0
         work%begin_cycle_random = .true.
         work%cycle_done_random = .false.
         do
            if( advance_cell_random_part( work, err ) &
                 == RETURN_FAIL ) goto 800
            if( .not. work%str_zero( work%cell ) ) then
               non_zero_cells_present = .true.
               sum = sum + mu( work%cell )
            end if
            if( work%cycle_done_random ) exit
         end do
         if( non_zero_cells_present .and. ( sum == 0.D0 ) ) goto 150
         work%begin_cycle_random = .true.
         work%cycle_done_random = .false.
         do
            if( advance_cell_random_part( work, err ) &
                 == RETURN_FAIL ) goto 800
            if( .not. work%str_zero( work%cell ) ) then
               prob( work%cell ) = mu( work%cell ) / sum
            end if
            if( work%cycle_done_random ) exit
         end do
         if( work%cycle_done_fixed ) exit
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
50    call err_handle(err, 1, &
            comment = "Argument mu has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function normalize_prob
   !##################################################################
    integer(kind=our_int) function advance_cell_fixed_part( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, jj, var_old
      character(len=*), parameter :: &
           subname = "advance_cell_fixed_part"
      ! begin
      answer = RETURN_FAIL
      !
      if( work%begin_cycle_fixed ) then
         work%cell_fixed_part = 0
         do j = 1, work%nvar
            if( .not. work%fixed(j) ) cycle
            work%var(j) = 1
         end do
         work%begin_cycle_fixed = .false.
      else
         if( work%cycle_done_fixed ) goto 150
         do j = 1, work%nvar
            ! find the first fixed variable j that is not maxed out
            if( .not. work%fixed(j) ) cycle
            if( work%var(j) == work%nlev(j) ) cycle
            work%var(j) = work%var(j) + 1
            work%cell_fixed_part = work%cell_fixed_part + work%cumprod(j)
            ! go back and reset each previous fixed variable to 1
            do jj = 1, j - 1
               if( .not. work%fixed(jj) ) cycle
               var_old = work%var(jj)
               work%var(jj) = 1
               work%cell_fixed_part = work%cell_fixed_part + &
                    ( work%var(jj) - var_old ) * work%cumprod(jj)
            end do
            exit
         end do
      end if
      work%cell = 1 + work%cell_fixed_part + work%cell_random_part
      work%cycle_done_fixed = .true.
      do j = 1, work%nvar
         if( .not. work%fixed(j) ) cycle
         if( work%var(j) /= work%nlev(j) ) then
            work%cycle_done_fixed = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_cell_fixed_part
   !##################################################################
    integer(kind=our_int) function advance_cell_random_part( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, jj, var_old
      character(len=*), parameter :: &
           subname = "advance_cell_random_part"
      ! begin
      answer = RETURN_FAIL
      !
      if( work%begin_cycle_random ) then
         work%cell_random_part = 0
         do j = 1, work%nvar
            if( work%fixed(j) ) cycle
            work%var(j) = 1
         end do
         work%begin_cycle_random = .false.
      else
         if( work%cycle_done_random ) goto 150
         do j = 1, work%nvar
            ! find the first random variable j that is not maxed out
            if( work%fixed(j) ) cycle
            if( work%var(j) == work%nlev(j) ) cycle
            work%var(j) = work%var(j) + 1
            work%cell_random_part = work%cell_random_part + work%cumprod(j)
            ! go back and reset each previous random variable to 1
            do jj = 1, j - 1
               if( work%fixed(jj) ) cycle
               var_old = work%var(jj)
               work%var(jj) = 1
               work%cell_random_part = work%cell_random_part + &
                    ( work%var(jj) - var_old ) * work%cumprod(jj)
            end do
            exit
         end do
      end if
      work%cell = 1 + work%cell_fixed_part + work%cell_random_part
      work%cycle_done_random = .true.
      do j = 1, work%nvar
         if( work%fixed(j) ) cycle
         if( work%var(j) /= work%nlev(j) ) then
            work%cycle_done_random = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_cell_random_part
   !##################################################################
    integer(kind=our_int) function run_estep( work, err, &
         use_flatten, use_prior_data, use_input_data, use_cell_means ) &
         result(answer)
      ! Performs a single E-step
      ! Distributes the frequencies in prior_data_freq and input_data_freq
      ! across the cells of the complete-data table ( work%freq ),
      ! given the current estimate of the cell probs ( work%prob )
      ! or cell means (work%mu)
      ! Also accumulates the loglikelihood in work%loglik, and
      ! the log-prior in work%logprior
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare optional args
      logical, intent(in), optional :: use_flatten
      logical, intent(in), optional :: use_prior_data
      logical, intent(in), optional :: use_input_data
      logical, intent(in), optional :: use_cell_means
      ! declare locals
      integer(kind=our_int) :: row, i, j
      real(kind=our_dble) :: ntotal, nprior, ndata, log_ntotal, sum
      logical :: use_prior_data_local, use_input_data_local, &
           use_flatten_local, use_cell_means_local
      character(len=*), parameter :: &
           subname = "run_estep"
      ! begin
      answer = RETURN_FAIL
      ! apply defaults
      if( present( use_flatten ) ) then
         use_flatten_local = use_flatten
      else
         use_flatten_local = .true.
      end if
      if( present( use_prior_data ) ) then
         use_prior_data_local = use_prior_data
      else
         use_prior_data_local = .true.
      end if
      if( present( use_input_data ) ) then
         use_input_data_local = use_input_data
      else
         use_input_data_local = .true.
      end if
      if( present( use_cell_means ) ) then
         if( use_cell_means .and. (work%model_type /= "log-linear") ) goto 20
         use_cell_means_local = use_cell_means
      else
         use_cell_means_local = .false.
      end if
      ! initialize counts etc.
      work%freq(:) = 0.D0
      work%logprior = 0.D0
      work%loglik = 0.D0
      nprior = 0.D0
      ndata = 0.D0
      if( work%model_type == "log-linear" ) then
         ! poisson correction
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            work%loglik = work%loglik - work%mu(i)
         end do
         ! ridge contribution
         sum = 0.D0
         do j = 1, work%p
            sum = sum + work%beta(j)**2
         end do
         work%logprior = work%logprior - work%ridge * sum / 2.D0
      end if
      ! apply flattening constants
      if( use_flatten_local ) then
         if( flatten_table(work, err, use_cell_means_local ) &
              == RETURN_FAIL ) goto 800
         nprior = nprior + work%flatten
      end if
      ! distribute prior data
      if( use_prior_data_local ) then
         do row = 1, work%nrow_prior_data
            if( run_estep_single_row( row, work%prior_data, &
                 work%prior_data_freq, work%prior_data_use, &
                 work, err, prior = .true., &
                 use_cell_means = use_cell_means_local ) &
                 == RETURN_FAIL ) goto 800
            if( work%prior_data_use(row) ) nprior = nprior + &
                 work%prior_data_freq(row)
         end do
      end if
      ! distribute input data
      if( use_input_data_local ) then
         do row = 1, work%nrow_input_data
            if( run_estep_single_row( row, work%input_data, &
                 work%input_data_freq, work%input_data_use, &
                 work, err, prior = .false., &
                 use_cell_means = use_cell_means_local ) &
                 == RETURN_FAIL ) goto 800 
            if( work%input_data_use(row) ) ndata = ndata + &
                 work%input_data_freq(row)
         end do
      end if
      ! poisson corrections for saturated model
      if( work%model_type == "saturated" ) then
         ntotal = nprior + ndata
         if( ntotal > 0.D0 ) then
            log_ntotal = log(ntotal)
            work%loglik = work%loglik + ndata * log_ntotal - ntotal
            work%logprior = work%logprior + nprior * log_ntotal
         else
            ! nothing to do, because nprior and ndata are zero
         end if
      end if
      work%logP = work%logprior + work%loglik
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_estep
   !##################################################################
    integer(kind=our_int) function flatten_table( work, err, &
         use_cell_means ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: use_cell_means
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: flat_const, log_prob
      logical :: use_cell_means_local
      character(len=*), parameter :: &
           subname = "flatten_table"
      ! begin
      answer = RETURN_FAIL
      if( present( use_cell_means ) ) then
         if( use_cell_means .and. (work%model_type /= "log-linear" ) ) goto 20
         use_cell_means_local = use_cell_means
      else
         use_cell_means_local = .false.
      end if
      flat_const = work%flatten / &
           real( work%ncells_nonzero, our_dble )
      if( flat_const == 0.D0 ) goto 10 ! nothing to do
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         work%freq(i) = work%freq(i) + flat_const
         if( use_cell_means_local ) then
            if( work%mu(i) <= 0.D0 ) goto 200
            log_prob = log( work%mu(i) )
         else
            if( work%prob(i) <= 0.D0 ) goto 200
            log_prob = log( work%prob(i) )
         end if
         work%logprior = work%logprior + flat_const * &
              log_prob
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function flatten_table
   !##################################################################
    integer(kind=our_int) function run_estep_single_row( row, &
         data_matrix, data_freq, data_use, work, err, &
         prior, use_cell_means ) result( answer )
      ! performs the e-step for a single row of data_matrix
      ! (data_matrix is either work%input_data or work%prior_data)
      ! (data_freq is either work%input_data_freq or work%prior_data_freq)
      ! (data_use is either work%input_data_use or work%prior_data_use)
      ! Distributes the frequency in data_freq(row) 
      ! across the cells of the complete-data table ( work%counts ),
      ! given the current estimate of the cell probs ( work%prob ),
      ! and across the cells of each submodel table
      ! If prior = .true., accumulates work%logprior
      ! If prior = .false., accumulates work%loglik (default)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: row 
      integer(kind=our_int), intent(in) :: data_matrix(:,:)
      real(kind=our_dble), intent(in) :: data_freq(:)
      logical, intent(in) :: data_use(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: prior
      logical, intent(in), optional :: use_cell_means
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: sum, freq_part
      logical :: prior_local, use_cell_means_local, &
           non_zero_cells_present
      character(len=*), parameter :: &
           subname = "run_estep_single_row"
      ! begin
      answer = RETURN_FAIL
      if( present( prior ) ) then
         prior_local = prior
      else
         prior_local = .false.
      end if
      if( present( use_cell_means ) ) then
         if( use_cell_means .and. (work%model_type /= "log-linear" ) ) goto 20
         use_cell_means_local = use_cell_means
      else
         use_cell_means_local = .false.
      end if
      if( ( row < 0 ) .or. ( row > size(data_matrix, 1) ) ) goto 100
      if( .not. data_use(row) ) goto 10 ! nothing to do
      !
      ! accumulate total probability for all cells of the 
      ! complete-data table that contribute to the given row
      ! of coarsened data
      work%begin_cycle = .true.
      work%cycle_done = .false.
      non_zero_cells_present = .false.
      sum = 0.D0
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            if( use_cell_means_local ) then
               sum = sum + work%mu( i )
            else
               sum = sum + work%prob( i )
            end if
            non_zero_cells_present = .true.
         end if
         if( work%cycle_done ) exit
      end do
      if( non_zero_cells_present .and. ( sum == 0.D0 ) ) goto 200
      if( .not. non_zero_cells_present .and. &
           ( data_freq(row) > 0.D0 ) ) goto 205
      if( sum < 0.D0 ) goto 210
      if( prior_local ) then
         if( non_zero_cells_present ) &
              work%logprior = work%logprior + data_freq(row) * log(sum)
      else
         if( non_zero_cells_present ) &
              work%loglik = work%loglik + data_freq(row) * log(sum)
      end if
      !
      ! distribute the frequency proportionately across the cells
      ! of the complete-data table, and also the cells of the
      ! table for each sub-model
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            if( use_cell_means_local ) then
               freq_part = data_freq( row ) * work%mu( i ) / sum
            else
               freq_part = data_freq( row ) * work%prob( i ) / sum
            end if
            work%freq( i ) = work%freq( i ) + freq_part
         end if
         if( work%cycle_done ) exit
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument row out of bounds" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Bad row in model frame, positive freq for zero cells" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_estep_single_row
   !##################################################################
    integer(kind=our_int) function advance_to_next_cell( row, &
         data_matrix, work, err ) result(answer)
      ! advances to the next cell of the complete-data table
      ! that contributes to row of data_matrix
      ! (data_matrix is either work%input_data or work%prior_data)
      ! ( if begin_cycle is .true., finds the first cell )
      !  The cell number is stored in work%cell, and the corresponding
      !  levels of the variables (base levels) are stored in work%var
      implicit none
      ! declare workspaces
      integer(kind=our_int), intent(in) :: row
      integer(kind=our_int), intent(in) :: data_matrix(:,:) 
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k, ii, var_old
      character(len=*), parameter :: &
           subname = "advance_to_next_cell"
      ! begin
      answer = RETURN_FAIL
      if( work%begin_cycle ) then
         work%var_max(:) = .false.
         work%cell = 1
         do i = 1, work%nvar                     ! i=variable
            j = data_matrix( row, i )   ! j=observed level
            k = 1                                ! k=index of base_level
            work%var_index(i) = k
            if( size( work%mapping%vec(i)%vec(j)%vec ) == k ) &
                 work%var_max(i) = .true.
            work%var(i) = work%mapping%vec(i)%vec(j)%vec(k)
            work%cell = work%cell + ( work%var(i) - 1 ) * work%cumprod(i)
         end do
         work%begin_cycle = .false.
      else
         if( work%cycle_done ) goto 150
         do i = 1, work%nvar
            ! find the first variable that has not maxed out
            if( work%var_max(i) ) cycle
            ! advance variable i to next base_level in the set
            j = data_matrix( row, i ) 
            k = work%var_index(i)
            var_old = work%var(i)
            k = k + 1
            work%var_index(i) = k
            if( size( work%mapping%vec(i)%vec(j)%vec ) == k ) &
               work%var_max(i) = .true.
            work%var(i) = work%mapping%vec(i)%vec(j)%vec(k)
            work%cell = work%cell + ( work%var(i) - var_old ) * &
                 work%cumprod(i)
            ! go back and reset each of the previous variables
            ! to the first level in its mapping set, unless the
            ! mapping set has only one member
            do ii = 1, i - 1
               j = data_matrix( row, ii ) 
               if( size( work%mapping%vec(ii)%vec(j)%vec ) == 1 ) cycle
               var_old = work%var(ii)
               k = 1
               work%var_index(ii) = k
               work%var_max(ii) = .false.
               work%var(ii) = work%mapping%vec(ii)%vec(j)%vec(k)
               work%cell = work%cell + ( work%var(ii) - var_old ) * &
                    work%cumprod(ii)
            end do
            exit
         end do
      end if
      work%cycle_done = .true.
      do i = 1, work%nvar
         if( .not. work%var_max(i) ) then
            work%cycle_done = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_to_next_cell
   !##################################################################
   integer(kind=our_int) function compute_ls_fit( y, beta, &
        work, err ) result(answer)
      ! regress y on model_matrix, put result into beta
      ! do not call this function if model_type = "saturated"
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: y(:)
      ! declare outputs
      real(kind=our_dble), intent(out) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k
      real(kind=our_dble) :: sum
      character(len=*), parameter :: subname = "compute_ls_fit"
      ! begin
      answer = RETURN_FAIL
      !
      if( work%model_type /= "log-linear" ) goto 20 
      if( size(y) /= work%n ) goto 50
      if( size(beta) /= work%p ) goto 100
      !
      do j = 1, work%p
         sum = 0.D0
         do i = 1, work%n
            if( work%str_zero(i) ) cycle
            sum = sum + work%model_matrix(i,j) * y(i)
         end do
         work%wkpA(j) = sum
         do k = 1, j
            sum = 0.D0
            do i = 1, work%n
               if( work%str_zero(i) ) cycle
               sum = sum + work%model_matrix(i,j) * &
                    work%model_matrix(i,k)
            end do
            work%wkppA(j,k) = sum
            work%wkppA(k,j) = work%wkppA(j,k)
         end do
      end do
      if( cholesky_in_place(work%wkppA, err ) &
           == RETURN_FAIL ) goto 750
      if( invert_lower(work%wkppA, err ) == RETURN_FAIL ) goto 750
      if( premult_lower_by_transpose( work%wkppA, &
           work%wkppB, err) == RETURN_FAIL ) goto 800
      do j = 1, work%p
         sum = 0.D0
         do k = 1, work%p
            sum = sum + work%wkppB(j,k) * work%wkpA(k)
         end do
         beta(j) = sum
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Argument y has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
750   call err_handle(err, 1, &
           comment = "Model matrix not full rank" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_ls_fit
   !##################################################################
   integer(kind=our_int) function run_mstep_log_linear( work, err ) &
        result(answer)
      ! fit a loglinear model to frequencies in work%freq, using
      ! work%beta as a starting value, and storing the result in 
      ! work%beta_new
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k, i
      real(kind=our_dble) :: sum
      character(len=*), parameter :: subname = "run_mstep_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      ! store beta in beta_tmp, then restore before exiting 
      work%beta_tmp(:) = work%beta(:) 
      work%iter_mstep = 0
      work%converged_mstep = .false.
      do
         if( work%converged_mstep ) exit
         if( work%iter_mstep >= work%iter_max_NR ) exit
         work%iter_mstep = work%iter_mstep + 1
         if(work%iter_mstep == 1 ) then
            if( compute_mu_from_beta( work%beta, work, err ) &
                 == RETURN_FAIL ) goto 800
         else
            work%beta(:) = work%beta_new(:)
            ! mu is already up-to-date
         end if
         if( compute_score_hessian_mstep( work, err ) &
              == RETURN_FAIL ) goto 800
         ! put negative hessian into wkppA and inverse into wkppB
         work%wkppA(:,:) = - work%hess_mstep(:,:)
         if( cholesky_in_place(work%wkppA, err ) &
              == RETURN_FAIL ) goto 200
         if( invert_lower(work%wkppA, err ) == RETURN_FAIL ) goto 200
         if( premult_lower_by_transpose( work%wkppA, &
              work%wkppB, err) == RETURN_FAIL ) goto 700
         do j = 1, work%p
            sum = 0.D0
            do k = 1, work%p
               sum = sum + work%wkppB(j,k) * work%score_mstep(k)
            end do
            work%wkpA(j) = work%beta(j) + sum
         end do
         work%beta_new(:) = work%wkpA(:)
         work%wknA(:) = work%mu(:)
         if( compute_mu_from_beta( work%beta_new, work, err ) &
                 == RETURN_FAIL ) goto 800
         ! detect convergence 
         work%max_diff_mstep = 0.D0
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            work%max_diff_mstep = max( work%max_diff_mstep, &
                 abs( work%mu(i) - work%wknA(i) ) )
         end do
         if( work%max_diff_mstep <= work%crit_NR ) &
              work%converged_mstep = .true.
      end do
      if( .not. work%converged_mstep ) goto 750 ! EM will stop
      work%beta(:) = work%beta_tmp(:)           ! EM will continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Hessian matrix not neg-def" )
      goto 700
700   continue
      call err_handle(err, 1, &
           comment = "Newton-Raphson aborted" )
      call err_handle(err, 5, iiter = work%iter_mstep )
      goto 800
750   call err_handle(err, 1, &
            comment = "Newton-Raphson failed to converge by" )
      call err_handle(err, 5, iiter = work%iter_mstep )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_mstep_log_linear
   !##################################################################
   integer(kind=our_int) function compute_score_hessian_mstep( &
        work, err ) result(answer)
      ! ordinary score and hessian for log-linear model, using the 
      ! frequencies stored in work%freq(:)
      ! assumes that work%mu and work%logmu are up-to-date
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_score_hessian_mstep"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      work%logP_mstep = 0.D0
      work%score_mstep(:) = 0.D0
      work%hess_mstep(:,:) = 0.D0
      ! ridge contribution
      sum = 0.D0
      do j = 1, work%p
         sum = sum + work%beta(j)**2
      end do
      work%logP_mstep = work%logP_mstep - work%ridge * sum / 2.D0
      do j = 1, work%p
         work%score_mstep(j) = work%score_mstep(j) &
              - work%ridge * work%beta(j)
         work%hess_mstep(j,j) = work%hess_mstep(j,j) - work%ridge
      end do
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         work%logP_mstep = work%logP_mstep - work%mu(i) + &
              work%freq(i) * work%logmu(i)
         do j = 1, work%p
            work%score_mstep(j) = work%score_mstep(j) + ( work%freq(i) &
                 - work%mu(i) ) * work%model_matrix(i,j)
            do k = 1, j
               work%hess_mstep(j,k) = work%hess_mstep(j,k) - work%mu(i) * &
                    work%model_matrix(i,j) * work%model_matrix(i,k)
            end do
         end do
      end do
      do j = 1, work%p
         do k = 1, j
            work%hess_mstep(k,j) = work%hess_mstep(j,k)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_score_hessian_mstep
   !##################################################################
   integer(kind=our_int) function structural_zero_check( &
        work, err, check_prob ) result(answer)
      ! detects observed or prior observations that are impossible
      ! given the supplied structural zeros.
      ! if check_probs = .true., checks work%probs to make sure there
      ! is zero prob for every structural zero (default is false)
      ! 
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: check_prob
      ! declare locals
      integer(kind=our_int) :: row, i
      logical :: check_prob_local
      character(len=*), parameter :: &
           subname = "structural_zero_check"
      ! begin
      answer = RETURN_FAIL
      if( present( check_prob ) ) then
         check_prob_local = check_prob
      else
         check_prob_local = .false.
      end if
      !
      do row = 1, work%nrow_prior_data
         if( str_zero_check_single_row( row, work%prior_data, &
              work%prior_data_freq, work, err ) &
              == RETURN_FAIL ) goto 100
      end do
      do row = 1, work%nrow_input_data
         if( str_zero_check_single_row( row, work%input_data, &
              work%input_data_freq, work, err ) &
              == RETURN_FAIL ) goto 200
      end do
      !
      if( check_prob_local ) then
         do i = 1, work%ncells
            if( work%str_zero(i) .and. work%prob(i) /= 0.D0 ) goto 300
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
100   call err_handle(err, 1, &
            comment = "Bad row in model frame mfPrior" )
      call err_handle(err, 3, iobs=row)
      goto 800
200   call err_handle(err, 1, &
            comment = "Bad row in model frame mfSeen" )
      call err_handle(err, 3, iobs=row)
      goto 800
300   call err_handle(err, 1, &
            comment = "Non-zero prob provided for structural zero" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function structural_zero_check
   !##################################################################
   integer(kind=our_int) function str_zero_check_single_row( row, &
        data_matrix, data_freq, work, err ) result(answer)
      implicit none
      integer(kind=our_int), intent(in) :: row 
      integer(kind=our_int), intent(in) :: data_matrix(:,:)
      real(kind=our_dble), intent(in) :: data_freq(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer :: i
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "str_zero_check_single_row"
      ! begin
      answer = RETURN_FAIL
      if( ( row < 0 ) .or. ( row > size(data_matrix, 1) ) ) goto 100
      if( row > size( data_freq ) ) goto 100
      if( data_freq(row) == 0.D0 )  goto 10 ! nothing to do
      work%begin_cycle = .true.
      work%cycle_done = .false.
      non_zero_cells_present = .false.
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            non_zero_cells_present = .true.
            exit
         end if
         if( work%cycle_done ) exit
      end do
      if( .not. non_zero_cells_present ) goto 200
      ! normal exit
10    continue
      answer = RETURN_SUCCESS
      goto 999
100   call err_handle(err, 1, &
            comment = "Argument row out of bounds" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Data inconsistent with structural zeros")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function str_zero_check_single_row
   !##################################################################
   integer(kind=our_int) function run_em_log_linear( work, err ) &
        result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      logical :: aborted
      character(len=*), parameter :: subname = "run_em_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      work%iter = 0
      work%converged = .false.
      aborted = .false.
      do
         if( work%converged ) exit
         if( work%iter >= work%iter_max_EM ) exit
         work%iter = work%iter + 1
         if( work%iter == 1 ) then
            if( compute_mu_from_beta( work%beta, work, err ) &
                 == RETURN_FAIL ) then
               aborted = .true.
               goto 600
            else
               work%beta(:) = work%beta_new(:)
               ! mu is already up-to-date
            end if
         end if
         if( normalize_prob( work%mu, work%prob, work, err ) &
              == RETURN_FAIL ) then
            aborted = .true.
            goto 600
         end if
         if( run_estep( work, err, &
              use_flatten = .true., use_prior_data = .true., &
              use_input_data = .true., use_cell_means = .true. ) &
              == RETURN_FAIL ) then
            aborted = .true.
            goto 600
         end if
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         if( work%iter == 1 ) work%start_logP = work%logP
         work%mu_old(:) = work%mu(:)
         if( run_mstep_log_linear( work, err ) == RETURN_FAIL ) then
            aborted = .true.
            goto 600
         end if
         ! At this point, mu is consistent with beta_new 
         ! detect convergence 
         work%max_diff = 0.D0
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            work%max_diff = max( work%max_diff, &
                 abs( work%mu(i) - work%mu_old(i) ) )
         end do
         if( work%max_diff <= work%crit_EM ) work%converged = .true.
      end do
600   continue
      if( aborted ) then
         call err_handle(err, 1, &
              comment = "EM algorithm aborted" )
         call err_handle(err, 5, iiter = work%iter )
         work%beta_new(:) = work%beta(:)
         work%prob_new(:) = work%prob(:)
      else
         ! put latest estimates implied by work%beta_new into work%prob_new
         if( compute_mu_from_beta( work%beta_new, work, err ) &
              == RETURN_FAIL ) goto 800
         if( normalize_prob( work%mu, work%prob_new, work, err ) &
              == RETURN_FAIL ) goto 800
         work%prob(:) = work%prob_new(:)
         if( run_estep_with_derivs( work%beta_new, work, err ) &
              == RETURN_FAIL ) goto 800
         if( compute_vhat_beta( work, err ) == RETURN_FAIL ) goto 800
      end if
      ! compute freq from observed data for putting into mfTrue
      work%input_data_use(:) = .true.
      if( run_estep( work, err, &
           use_flatten = .false., use_prior_data = .false., &
           use_input_data = .true., use_cell_means = .true. ) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_em_log_linear
   !##################################################################
   integer(kind=our_int) function run_em_saturated( work, err ) &
        result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      logical :: aborted, boundary
      character(len=*), parameter :: subname = "run_em_saturated"
      ! begin
      answer = RETURN_FAIL
      work%iter = 0
      work%converged = .false.
      aborted = .false.
      do
         if( work%converged ) exit
         if( work%iter >= work%iter_max_EM ) exit
         work%iter = work%iter + 1
         if( work%iter > 1 ) work%prob(:) = work%prob_new(:)
         if( run_estep( work, err, &
              use_flatten = .true., use_prior_data = .true., &
              use_input_data = .true.) == RETURN_FAIL ) then
            aborted = .true.
            goto 600
         end if
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         if( work%iter == 1 ) work%start_logP = work%logP
         if( normalize_prob( work%freq, work%prob_new, &
              work, err ) == RETURN_FAIL ) then
            aborted = .true.
            goto 600
         end if
         ! detect convergence 
         work%max_diff = 0.D0
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            work%max_diff = max( work%max_diff, &
                 abs( work%prob_new(i) - work%prob(i) ) )
         end do
         if( work%max_diff <= work%crit_EM ) work%converged = .true.
      end do
      ! normal exit
600   continue
      if( aborted ) then
         call err_handle(err, 1, &
              comment = "EM algorithm aborted" )
         call err_handle(err, 5, iiter = work%iter )
         work%prob_new(:) = work%prob(:)
      else
         work%prob(:) = work%prob_new(:)
      end if
      ! compute freq
      work%input_data_use(:) = .true.
      if( run_estep( work, err, &
           use_flatten = .false., use_prior_data = .false., &
           use_input_data = .true.) == RETURN_FAIL ) goto 800
      ! check boundary
      boundary = .false.
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         if( work%prob_new(i) <= work%crit_boundary ) then
            boundary = .true.
            exit
         end if
      end do
      if( boundary ) then
         call err_handle(err, 1, &
              comment = "Estimate at or near boundary" )
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_em_saturated
   !##################################################################
   integer(kind=our_int) function compute_lambda_from_beta( beta, &
        lambda, work, err ) result(answer)
      ! only call this function if model_type = "log-linear"
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare outputs
      real(kind=our_dble), intent(out) :: lambda(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: rtmp
      character(len=*), parameter :: subname = "compute_lambda_from_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( size(beta) /= work%p ) goto 50
      if( size(lambda) /= work%p ) goto 60
      if( compute_mu_from_beta( beta, work, err ) == RETURN_FAIL ) goto 800
      if( normalize_prob( work%mu, work%wknA, work, err ) &
           == RETURN_FAIL ) goto 800
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         rtmp = work%wknA(i)
         if( rtmp <= 0.D0 ) goto 100
         work%wknA(i) = log(rtmp) - work%offset(i)
      end do
      if( compute_ls_fit( work%wknA, lambda, work, err ) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
60    call err_handle(err, 1, &
            comment = "Argument lambda has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_lambda_from_beta
   !##################################################################
    integer(kind=our_int) function run_estep_with_derivs( beta, &
         work, err ) result(answer)
      ! Performs a single E-step, computing first and second derivatives
      ! of the observed-data loglikelihood wrt beta uder a Poisson model
      ! Distributes the frequencies in prior_data_freq and input_data_freq
      ! across the cells of the complete-data table ( work%freq ),
      ! given the current estimate of the cell means ( work%mu )
      ! Important: work%mu must be consistent with current beta
      ! Note that after this is run, work%loglik contains the 
      ! value of the loglikelihood under the Poisson surrogate model,
      ! not the multinomial model
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare optional args
      ! declare locals
      integer(kind=our_int) :: row, i, j, k
      real(kind=our_dble) :: sum, rtmp
      character(len=*), parameter :: &
           subname = "run_estep_with_derivs"
      ! begin
      answer = RETURN_FAIL
      ! apply defaults
      if( work%model_type /= "log-linear" ) goto 20
      ! initialize counts etc.
      work%logprior = 0.D0
      work%loglik = 0.D0
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         work%loglik = work%loglik - work%mu(i)
      end do
      work%fhat(:) = 0.D0
      work%Mx(:,:) = 0.D0
      ! ridge contribution to prior
      sum = 0.D0
      do j = 1, work%p
         sum = sum + beta(j)**2
      end do
      work%logprior = work%logprior - work%ridge * sum / 2.D0
      ! apply flattening constants
      if( flatten_table_with_derivs(work, err ) &
           == RETURN_FAIL ) goto 800
      do row = 1, work%nrow_prior_data
         if( run_estep_single_row_with_derivs( row, work%prior_data, &
              work%prior_data_freq, work%prior_data_use, &
              work, err, prior = .true. ) &
              == RETURN_FAIL ) goto 800
      end do
      do row = 1, work%nrow_input_data
         if( run_estep_single_row_with_derivs( row, work%input_data, &
              work%input_data_freq, work%input_data_use, &
              work, err, prior = .false. ) &
              == RETURN_FAIL ) goto 800
      end do
      work%logP = work%logprior + work%loglik
      !
      work%score(:) = 0.D0
      work%hess(:,:) = 0.D0
      ! ridge contributions
      do j = 1, work%p
         work%score(j) = work%score(j) - work%ridge * beta(j)
         work%hess(j,j) = work%hess(j,j) - work%ridge
      end do
      ! add score and Hessian term 1
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         rtmp = work%fhat(i) - work%mu(i)
         do j = 1, work%p
            work%score(j) = work%score(j) + rtmp * work%model_matrix(i,j)
            do k = 1, j
               work%hess(j,k) = work%hess(j,k) + &
                    rtmp * work%model_matrix(i,j) * work%model_matrix(i,k)
            end do
         end do
      end do
      ! add Hessian term 2
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         do j = 1, work%p
            do k = 1, j
               work%hess(j,k) = work%hess(j,k) + &
                    work%model_matrix(i,j) * work%Mx(i,k)
            end do
         end do
      end do
      do j = 1, work%p
         do k = 1, j
            work%hess(k,j) = work%hess(j,k)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_estep_with_derivs
   !##################################################################
    integer(kind=our_int) function flatten_table_with_derivs( work, err ) &
         result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      ! declare locals
      integer(kind=our_int) :: i, j
      real(kind=our_dble) :: flat_const, log_mean
      character(len=*), parameter :: &
           subname = "flatten_table_with_derivs"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      flat_const = work%flatten / &
           real( work%ncells_nonzero, our_dble )
      if( flat_const == 0.D0 ) goto 10 ! nothing to do
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         ! increment fhat and logprior
         work%fhat(i) = work%fhat(i) + flat_const
         if( work%mu(i) <= 0.D0 ) goto 200
         log_mean = log( work%mu(i) )
         work%logprior = work%logprior + flat_const * log_mean
         ! increment Mx
         do j = 1, work%p
            work%Mx(i,j) = work%Mx(i,j) - flat_const * &
                 work%model_matrix(i,j)
         end do
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      call err_handle(err, 15, icell = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function flatten_table_with_derivs
   !##################################################################
    integer(kind=our_int) function run_estep_single_row_with_derivs( &
         row, data_matrix, data_freq, data_use, work, err, &
         prior ) result( answer )
      ! performs the e-step for a single row of data_matrix
      ! (data_matrix is either work%input_data or work%prior_data)
      ! (data_freq is either work%input_data_freq or work%prior_data_freq)
      ! (data_use is either work%input_data_use or work%prior_data_use)
      ! Distributes the frequency in data_freq(row) 
      ! across the cells of the complete-data table ( work%counts ),
      ! given the current estimate of the cell means ( work%mu ),
      ! and across the cells of each submodel table
      ! If prior = .true., accumulates work%logprior
      ! If prior = .false., accumulates work%loglik (default)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: row 
      integer(kind=our_int), intent(in) :: data_matrix(:,:)
      real(kind=our_dble), intent(in) :: data_freq(:)
      logical, intent(in) :: data_use(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: prior
      ! declare locals
      integer(kind=our_int) :: i, ii, j
      real(kind=our_dble) :: tau, freq_part, sum
      logical :: prior_local, non_zero_cells_present
      character(len=*), parameter :: &
           subname = "run_estep_single_row_with_derivs"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      if( present( prior ) ) then
         prior_local = prior
      else
         prior_local = .false.
      end if
      if( ( row < 0 ) .or. ( row > size(data_matrix, 1) ) ) goto 100
      if( .not. data_use(row) ) goto 10 ! nothing to do
      !
      ! accumulate total mean for all cells of the 
      ! complete-data table that contribute to the given row
      ! of coarsened data
      work%begin_cycle = .true.
      work%cycle_done = .false.
      non_zero_cells_present = .false.
      tau = 0.D0
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            tau = tau + work%mu( i )
            non_zero_cells_present = .true.
         end if
         if( work%cycle_done ) exit
      end do
      if( non_zero_cells_present .and. ( tau == 0.D0 ) ) goto 200
      if( .not. non_zero_cells_present .and. &
           ( data_freq(row) > 0.D0 ) ) goto 205
      if( tau < 0.D0 ) goto 210
      if( prior_local ) then
         if( non_zero_cells_present ) &
              work%logprior = work%logprior + data_freq(row) * log(tau)
      else
         if( non_zero_cells_present ) &
              work%loglik = work%loglik + data_freq(row) * log(tau)
      end if
      ! compute bigFhat and increment fhat
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            freq_part = data_freq( row ) * work%mu( i ) / tau
            work%bigFhat(i) = freq_part
            work%fhat(i) = work%fhat(i) + freq_part
         end if
         if( work%cycle_done ) exit
      end do
      ! increment Mx
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) then
            work%begin_cycle_2 = .true.
            work%cycle_done_2 = .false.
            do
               if( advance_to_next_cell_2( row, data_matrix, work, err ) &
                    == RETURN_FAIL ) goto 800
               ii = work%cell_2
               if( .not. work%str_zero(ii) ) then
                  work%Mrow(ii) = - work%bigFhat(i) * work%bigFhat(ii) &
                       / data_freq(row)
               end if
               if( work%cycle_done_2 ) exit
            end do
            !
            do j = 1, work%p
               sum = 0.D0
               work%begin_cycle_2 = .true.
               work%cycle_done_2 = .false.
               do
                  if( advance_to_next_cell_2( row, data_matrix, work, err ) &
                       == RETURN_FAIL ) goto 800
                  ii = work%cell_2
                  if( .not. work%str_zero(ii) ) then
                     sum = sum + work%Mrow(ii) * work%model_matrix(ii,j)
                  end if
                  if( work%cycle_done_2 ) exit
               end do
               work%Mx(i,j) = work%Mx(i,j) + sum
            end do
            !
            if( work%cycle_done ) exit
         end if
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument row out of bounds" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Bad row in model frame, positive freq for zero cells" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_estep_single_row_with_derivs
   !##################################################################
    integer(kind=our_int) function advance_to_next_cell_2( row, &
         data_matrix, work, err ) result(answer)
      ! advances to the next cell of the complete-data table
      ! that contributes to row of data_matrix
      ! (data_matrix is either work%input_data or work%prior_data)
      ! ( if begin_cycle is .true., finds the first cell )
      !  The cell number is stored in work%cell, and the corresponding
      !  levels of the variables (base levels) are stored in work%var
      implicit none
      ! declare workspaces
      integer(kind=our_int), intent(in) :: row
      integer(kind=our_int), intent(in) :: data_matrix(:,:) 
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k, ii, var_old
      character(len=*), parameter :: &
           subname = "advance_to_next_cell_2"
      ! begin
      answer = RETURN_FAIL
      if( work%begin_cycle_2 ) then
         work%var_max_2(:) = .false.
         work%cell_2 = 1
         do i = 1, work%nvar                     ! i=variable
            j = data_matrix( row, i )   ! j=observed level
            k = 1                                ! k=index of base_level
            work%var_index_2(i) = k
            if( size( work%mapping%vec(i)%vec(j)%vec ) == k ) &
                 work%var_max_2(i) = .true.
            work%var_2(i) = work%mapping%vec(i)%vec(j)%vec(k)
            work%cell_2 = work%cell_2 + ( work%var_2(i) - 1 ) * work%cumprod(i)
         end do
         work%begin_cycle_2 = .false.
      else
         if( work%cycle_done_2 ) goto 150
         do i = 1, work%nvar
            ! find the first variable that has not maxed out
            if( work%var_max_2(i) ) cycle
            ! advance variable i to next base_level in the set
            j = data_matrix( row, i ) 
            k = work%var_index_2(i)
            var_old = work%var_2(i)
            k = k + 1
            work%var_index_2(i) = k
            if( size( work%mapping%vec(i)%vec(j)%vec ) == k ) &
               work%var_max_2(i) = .true.
            work%var_2(i) = work%mapping%vec(i)%vec(j)%vec(k)
            work%cell_2 = work%cell_2 + ( work%var_2(i) - var_old ) * &
                 work%cumprod(i)
            ! go back and reset each of the previous variables
            ! to the first level in its mapping set, unless the
            ! mapping set has only one member
            do ii = 1, i - 1
               j = data_matrix( row, ii ) 
               if( size( work%mapping%vec(ii)%vec(j)%vec ) == 1 ) cycle
               var_old = work%var_2(ii)
               k = 1
               work%var_index_2(ii) = k
               work%var_max_2(ii) = .false.
               work%var_2(ii) = work%mapping%vec(ii)%vec(j)%vec(k)
               work%cell_2 = work%cell_2 + ( work%var_2(ii) - var_old ) * &
                    work%cumprod(ii)
            end do
            exit
         end do
      end if
      work%cycle_done_2 = .true.
      do i = 1, work%nvar
         if( .not. work%var_max_2(i) ) then
            work%cycle_done_2 = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_to_next_cell_2
   !##################################################################
    integer(kind=our_int) function compute_vhat_beta( work, err ) &
         result(answer)
      ! run this after run_estep_with_derive
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare optional args
      ! declare locals
      logical :: failed, boundary
      integer(kind=our_int) :: i
      character(len=*), parameter :: &
           subname = "compute_vhat_beta"
      ! begin
      answer = RETURN_FAIL
      ! apply defaults
      if( work%model_type /= "log-linear" ) goto 20
      work%wkppA(:,:) = - work%hess(:,:)
      failed = .false.
      if( cholesky_in_place(work%wkppA, err ) &
           == RETURN_FAIL ) failed = .true.
      if( .not. failed ) then
         if( invert_lower(work%wkppA, err ) &
              == RETURN_FAIL ) failed = .true.
      end if
      if( failed ) then
         call err_handle(err, 1, &
              comment = "logP is not concave" )
         call err_handle(err, 1, &
              comment = "Estimated cov. matrix for beta not available" )
      end if
      if( .not. failed ) then
         if( premult_lower_by_transpose( work%wkppA, &
              work%vhat_beta, err) == RETURN_FAIL ) goto 800
      end if
      boundary = .false.
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         if( work%mu(i) <= work%crit_boundary ) then
            boundary = .true.
            exit
         end if
      end do
      if( boundary ) then
         call err_handle(err, 1, &
              comment = "Estimate at or near boundary" )
         if( .not. failed ) call err_handle(err, 1, &
              comment = "Estimated variances may be unreliable" )
      end if
      work%vhat_ok = .not. failed
      if( work%vhat_ok ) then
         work%cfac_vhat(:,:) = work%vhat_beta(:,:)
         if( cholesky_in_place(work%cfac_vhat, err ) &
              == RETURN_FAIL ) then
            work%vhat_ok = .false.
            call err_handle(err, 1, &
                 comment = "Matrix vhat_beta not positive definite" )
            call err_handle(err, 1, &
                 comment = "Standard errors cannot be computed" )
            work%vhat_ok = .false.
         end if
      end if
      if( .not. work%vhat_ok ) work%vhat_beta(:,:) = 0.D0
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_vhat_beta
   !##################################################################
   integer(kind=our_int) function run_cvam_estimate_em( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        estimate_info, estimate_var_info, skip_SEs_int, &
        ! workspaces
        work, err, &
        ! outputs
        packed_estimates, packed_SEs &
        ) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      integer(kind=our_int), intent(in) :: estimate_info(:,:)
      integer(kind=our_int), intent(in) :: estimate_var_info(:,:)
      integer(kind=our_int), intent(in) :: skip_SEs_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! outputs
      real(kind=our_dble), optional, intent(out) :: packed_estimates(:)
      real(kind=our_dble), optional, intent(out) :: packed_SEs(:)
      ! declare locals
      integer(kind=our_int) :: ijunk
      logical :: skip_SEs
      character(len=*), parameter :: &
           subname = "run_cvam_estimate_em"
      ! begin
      answer = RETURN_FAIL
      if( skip_SEs_int == 1 ) then
         skip_SEs = .true.
      else if( skip_SEs_int == 0 ) then
         skip_SEs = .false.
      else
         goto 20
      end if
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( input_data, input_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_params_into_workspace( prob, beta, vhat_beta, work, err ) &
           == RETURN_FAIL ) goto 800
      work%vhat_ok = .true.
      if( ( work%model_type == "log-linear" ) .and. ( .not. skip_SEs ) ) then
         work%cfac_vhat(:,:) = work%vhat_beta(:,:)
         if( cholesky_in_place(work%cfac_vhat, err ) &
              == RETURN_FAIL ) then
            work%vhat_ok = .false.
            call err_handle(err, 1, &
                 comment = "Matrix vhat_beta not positive definite" )
            call err_handle(err, 1, &
                 comment = "Standard errors cannot be computed" )
         end if
      end if
      if( create_estimate_objects( estimate_info, &
           estimate_var_info, work, err ) == RETURN_FAIL ) goto 800
      if( compute_estimates( work%prob, work, err ) &
           == RETURN_FAIL ) goto 800
      if( size(packed_estimates) /= work%n_packed_estimates ) goto 100
      packed_estimates(:) = work%packed_estimates(:)
      if( present( packed_SEs ) .and. &
           ( work%model_type == "log-linear" ) .and. ( .not. skip_SEs ) ) then
         if( size( packed_SEs ) /= work%n_packed_estimates ) goto 110
         if( compute_estimate_SEs( work%prob, work, err ) &
              == RETURN_FAIL ) goto 800
         packed_SEs(:) = work%packed_SEs(:)
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "Value of skip_SEs_int not recognized" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument packed_estimates has incorrect size" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Argument packed_SEs has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_estimate_em
   !##################################################################
   integer(kind=our_int) function put_params_into_workspace( &
        prob, beta, vhat_beta, work, err) result(answer)
      implicit none
      ! declare inputs
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "put_params_into_workspace"
      ! begin
      answer = RETURN_FAIL
      ! model_matrix
      if( size(prob) /= work%ncells ) goto 10
      if( work%model_type == "saturated" ) then
         if( normalize_prob( prob, work%prob, work, err ) &
              == RETURN_FAIL ) goto 800
      else if( work%model_type == "log-linear" ) then
         if( size( beta ) /= work%p ) goto 20
         work%beta(:) = beta(:)
         if( compute_mu_from_beta( work%beta, work, err ) &
              == RETURN_FAIL ) goto 800
         if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 800
         if( size( vhat_beta, 1 ) /= work%p ) goto 30
         if( size( vhat_beta, 2 ) /= work%p ) goto 30
         work%vhat_beta(:,:) = vhat_beta(:,:)
      else
         goto 50
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
10    call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
20    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument vhat_beta has incorrect size" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Value of model_type not recognized" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function put_params_into_workspace
   !##################################################################
    integer(kind=our_int) function create_estimate_objects( &
         estimate_info, estimate_var_info, work, err ) result(answer)
      ! additional workspace objects
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: estimate_info(:,:)
      integer(kind=our_int), intent(in) :: estimate_var_info(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, status, posn, nvar_tot, &
           n_estimates, type_int, nvar, st, fin, n_packed_estimates, &
           cell, ncells
      integer(kind=our_int), allocatable :: data_matrix(:,:)
      character(len=*), parameter :: &
           subname = "create_estimate_objects"
      ! begin
      answer = RETURN_FAIL
      ! set up estimates
      if( size( estimate_info, 2 ) /= 4 ) goto 200
      work%n_estimates = size( estimate_info, 1 )
      allocate( work%estimates(work%n_estimates), stat=status )
      if( status /= 0 ) goto 100
      nvar_tot = 0
      n_estimates = 0
      n_packed_estimates = 0
      do i = 1, work%n_estimates
         type_int = estimate_info(i,1)
         work%estimates(i)%estimate_type_int = type_int
         if( type_int == 1 ) then
            work%estimates(i)%estimate_type = "default"
         else
            goto 205
         end if
         work%estimates(i)%estimate_number = i
         nvar = estimate_info(i,2)
         if( nvar == 0) goto 210
         work%estimates(i)%nvar = nvar
         nvar_tot = nvar_tot + nvar
         st = estimate_info(i,3)
         fin = estimate_info(i,4)
         if( fin < st ) goto 220
         work%estimates(i)%first_in_estimates = st
         work%estimates(i)%last_in_estimates = fin
         n_packed_estimates = n_packed_estimates + ( fin - st + 1 )
      end do
      work%n_packed_estimates = n_packed_estimates
      allocate( work%packed_estimates(work%n_packed_estimates), &
           work%packed_estimates_mean(work%n_packed_estimates), &
           stat=status )
      if( status /= 0 ) goto 100
      work%packed_estimates(:) = 0.D0
      work%packed_estimates_mean(:) = 0.D0
      if( work%model_type == "log-linear" ) then
         allocate( work%packed_SEs(work%n_packed_estimates), &
              stat=status )
         if( status /= 0 ) goto 100
         work%packed_SEs(:) = 0.D0
      end if
      !
      if( size( estimate_var_info, 1 ) /= nvar_tot ) goto 300
      if( size( estimate_var_info, 2 ) /= 3 ) goto 300
      !
      posn = 0
      do i = 1, work%n_estimates
         nvar = work%estimates(i)%nvar
         allocate( work%estimates(i)%nlev(nvar), &
              work%estimates(i)%model_posn(nvar), &
              work%estimates(i)%fixed(nvar), &
              stat = status )
         if( status /= 0 ) goto 100
         do j = 1, nvar
            posn = posn + 1
            work%estimates(i)%nlev(j) = estimate_var_info(posn,1)
            work%estimates(i)%model_posn(j) = estimate_var_info(posn,2)
            if( estimate_var_info(posn,3) == 1 ) then
               work%estimates(i)%fixed(j) = .false.
            else if( estimate_var_info(posn,3) == 2 ) then
               work%estimates(i)%fixed(j) = .true.
            else
               goto 320
            end if
         end do
      end do
      !
      do i = 1, work%n_estimates
         nvar = work%nvar ! number of vars in full model
         allocate( work%estimates(i)%var_present(nvar), &
              work%estimates(i)%estimate_posn(nvar), &
              stat = status )
         if( status /= 0 ) goto 100
         work%estimates(i)%var_present(:) = .false.
         work%estimates(i)%estimate_posn(:) = 0
         work%estimates(i)%nvar_present = 0
         work%estimates(i)%nvar_absent = 0
         do j = 1, work%estimates(i)%nvar
            posn = work%estimates(i)%model_posn(j)
            work%estimates(i)%var_present(posn) = .true.
            work%estimates(i)%estimate_posn(posn) = j
         end do
         work%estimates(i)%nvar_present = work%estimates(i)%nvar
         work%estimates(i)%nvar_absent = work%nvar - work%estimates(i)%nvar
      end do
      !
      do i = 1, work%n_estimates
         nvar = work%estimates(i)%nvar
         allocate( work%estimates(i)%cumprod(nvar), stat=status )
         if( status /= 0 ) goto 100
         do j = 1, nvar
            if( j == 1 ) then
               work%estimates(i)%cumprod(j) = 1
            else
               work%estimates(i)%cumprod(j) = &
                    work%estimates(i)%cumprod(j-1) * &
                    work%estimates(i)%nlev(j-1)
            end if
         end do
         if( nvar == 0 ) then
            work%estimates(i)%ncells = 1
         else
            work%estimates(i)%ncells = &
                 work%estimates(i)%cumprod( nvar ) * &
                 work%estimates(i)%nlev( nvar )
         end if
         ncells = work%estimates(i)%ncells
         st = work%estimates(i)%first_in_estimates
         fin = work%estimates(i)%last_in_estimates
         if( ncells /= ( fin - st + 1 ) ) goto 350
         allocate( work%estimates(i)%prob(ncells), stat=status )
         if( status /= 0 ) goto 100
         work%estimates(i)%prob(:) = 0.D0
         allocate( work%estimates(i)%str_zero(ncells), stat=status )
         if( status /= 0 ) goto 100
         work%estimates(i)%str_zero(:) = .false.
         if( work%model_type == "log-linear" ) then
            allocate( work%estimates(i)%SEs(ncells), stat=status )
            if( status /= 0 ) goto 100
            work%estimates(i)%SEs(:) = 0.D0
            allocate( work%estimates(i)%xpi(ncells, work%p), stat=status )
            if( status /= 0 ) goto 100
            work%estimates(i)%xpi(:,:) = 0.D0
         end if
         allocate( work%estimates(i)%var(nvar), stat=status )
         if( status /= 0 ) goto 100
         work%estimates(i)%var(:) = 0
      end do
      ! 
      ! create a dummy data_matrix with a single row and missing
      ! values for all variables
      allocate( data_matrix(1, work%nvar), stat=status )
      if( status /= 0 ) goto 100
      data_matrix(1,:) = work%mvcode(:)
      ! create cells(:)
      do i = 1, work%n_estimates
         allocate( work%estimates(i)%cells( work%ncells ), stat=status )
         if( status /= 0 ) goto 100
         work%estimates(i)%cells(:) = 0
         work%begin_cycle = .true.
         work%cycle_done = .false.
         do
            if( advance_to_next_cell( 1, data_matrix, work, err ) &
                 == RETURN_FAIL ) goto 800
            cell = 1
            do j = 1, work%nvar
               if( .not. work%estimates(i)%var_present(j) ) cycle
               posn = work%estimates(i)%estimate_posn(j)
               cell = cell + ( work%var(j) - 1 ) * &
                    work%estimates(i)%cumprod(posn) 
            end do
            work%estimates(i)%cells( work%cell ) = cell
            if( work%cycle_done ) exit
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Argument estimate_info has incorrect size" )
      goto 800
205   call err_handle(err, 1, &
            comment = "estimate_type not recognized" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Estimate has no variables" )
      goto 500
220   call err_handle(err, 1, &
            comment = "Start, finish positions for estimate incorrect")
      goto 500
300   call err_handle(err, 1, &
            comment = "Argument estimate_var_info has incorrect size" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Invalid value in estimate_var_info, column 3" )
      goto 800
350   call err_handle(err, 1, &
            comment = "Incorrect number of cells for estimate" )
      goto 500
500   call err_handle(err, 12, iest = i)
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(data_matrix) ) deallocate( data_matrix ) 
    end function create_estimate_objects
   !##################################################################
   integer(kind=our_int) function compute_estimates( &
         prob, work, err ) result( answer )
      ! computes the estimates
      ! given the probabilities for the 
      ! complete-data table in work%prob, which may have already been 
      ! conditioned on the variables fixed in the model, but not
      ! necessarily all variables that are fixed in the estimate
      implicit none
      real(kind=our_dble), intent(in) :: prob(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, cfull, cest, st, fin
      real(kind=our_dble) :: sum
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "compute_estimates"
      ! begin
      answer = RETURN_FAIL
      if( size(prob) /= work%ncells ) goto 50
      do i = 1, work%n_estimates
         ! marginalize over the variables not present
         work%estimates(i)%prob(:) = 0.D0
         work%estimates(i)%str_zero(:) = .true.
         do cfull = 1, work%ncells
            if( work%str_zero(cfull) ) cycle
            cest = work%estimates(i)%cells( cfull )
            work%estimates(i)%str_zero( cest ) = .false.
            work%estimates(i)%prob( cest ) = &
                 work%estimates(i)%prob( cest ) + prob( cfull )
         end do
         ! condition on the variables that are fixed in the estimate,
         ! which must include any variables fixed in the model, and
         ! possibly some additional ones
         work%estimates(i)%begin_cycle_fixed = .true.
         work%estimates(i)%cycle_done_fixed = .false.
         do
            if( advance_cell_fixed_part_estimate( &
                 work%estimates(i), err ) &
                 == RETURN_FAIL ) goto 800
            non_zero_cells_present = .false.
            sum = 0.D0
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  non_zero_cells_present = .true.
                  sum = sum + work%estimates(i)%prob( cest )
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            if( non_zero_cells_present .and. ( sum == 0.D0 ) ) goto 100
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  work%estimates(i)%prob( cest ) = &
                       work%estimates(i)%prob( cest ) / sum
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            if( work%estimates(i)%cycle_done_fixed ) exit
         end do
         st = work%estimates(i)%first_in_estimates
         fin = work%estimates(i)%last_in_estimates 
         work%packed_estimates( st : fin ) = work%estimates(i)%prob(:)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
50    call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Non-positive probability in denominator" )
      call err_handle(err, 12, iest = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_estimates
    !##################################################################
   integer(kind=our_int) function compute_estimate_SEs( prob, &
         work, err ) result( answer )
      ! computes the SEs, make sure that compute_estimates was already
      ! done. The cell probabilities in prob may have been conditioned
      ! on all variables fixed in the model.
      implicit none
      real(kind=our_dble), intent(in) :: prob(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, cfull, cest, st, fin, j, jj, k
      real(kind=our_dble) :: sum
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "compute_estimate_SEs"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      if( .not. work%vhat_ok ) then
         goto 700 ! nothing to do
      end if
      ! copy fixed into fixed_tmp, and restore at the end
      work%fixed_tmp(:) = work%fixed(:)
      do i = 1, work%n_estimates
         ! for each variable fixed in the estimate, set element
         ! of fised to .true.
         work%fixed(:) = .false.
         do j = 1, work%estimates(i)%nvar
            jj = work%estimates(i)%model_posn(j)
            if( work%estimates(i)%fixed(j) ) work%fixed(jj) = .true.
         end do
         ! conditioh prob on all variables fixed in the estimate,
         ! putting result into wkNA
         work%wknA(:) = 0.D0
         work%begin_cycle_fixed = .true.
         work%cycle_done_fixed = .false.
         do
            if( advance_cell_fixed_part( work, err ) &
                 == RETURN_FAIL ) goto 800
            non_zero_cells_present = .false.
            sum = 0.D0
            work%begin_cycle_random = .true.
            work%cycle_done_random = .false.
            do
               if( advance_cell_random_part( work, err ) &
                    == RETURN_FAIL ) goto 800
               if( .not. work%str_zero( work%cell ) ) then
                  non_zero_cells_present = .true.
                  sum = sum + prob( work%cell )
               end if
               if( work%cycle_done_random ) exit
            end do
            if( non_zero_cells_present .and. ( sum == 0.D0 ) ) goto 150
            work%begin_cycle_random = .true.
            work%cycle_done_random = .false.
            do
               if( advance_cell_random_part( work, err ) &
                    == RETURN_FAIL ) goto 800
               if( .not. work%str_zero( work%cell ) ) then
                  work%wknA( work%cell ) = prob( work%cell ) / sum
               end if
               if( work%cycle_done_random ) exit
            end do
            if( work%cycle_done_fixed ) exit
         end do
         ! accumulate xpi
         work%estimates(i)%xpi(:,:) = 0.D0
         work%estimates(i)%str_zero(:) = .true.
         do cfull = 1, work%ncells
            if( work%str_zero(cfull) ) cycle
            cest = work%estimates(i)%cells( cfull )
            work%estimates(i)%str_zero( cest ) = .false.
            do j = 1, work%p
               work%estimates(i)%xpi( cest, j ) = &
                    work%estimates(i)%xpi( cest, j ) + &
                    work%model_matrix( cfull, j ) * work%wknA( cfull )
            end do
         end do
         ! compute dvec and SE
         work%estimates(i)%begin_cycle_fixed = .true.
         work%estimates(i)%cycle_done_fixed = .false.
         do
            if( advance_cell_fixed_part_estimate( &
                 work%estimates(i), err ) &
                 == RETURN_FAIL ) goto 800
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            work%wkpA(:) = 0.D0
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  do j = 1, work%p
                     work%wkpA(j) = work%wkpA(j) + &
                          work%estimates(i)%xpi( cest, j )
                  end do
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  do j = 1, work%p
                     work%dvec(j) = work%estimates(i)%xpi( cest, j ) - &
                          work%wkpA(j) * work%estimates(i)%prob( cest )
                  end do
                  ! premultiply dvec by t(cfac_vhat)
                  do j = 1, work%p
                     sum = 0.D0
                     do k = j, work%p
                        sum = sum + work%cfac_vhat(k,j) * work%dvec(k)
                     end do
                     work%wkpB(j) = sum
                  end do
                  ! take norm
                  sum = 0.D0
                  do j = 1, work%p
                     sum = sum + work%wkpB(j)**2
                  end do
                  work%estimates(i)%SEs( cest ) = sqrt( sum )
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            if( work%estimates(i)%cycle_done_fixed ) exit
         end do
         st = work%estimates(i)%first_in_estimates
         fin = work%estimates(i)%last_in_estimates 
         work%packed_SEs( st : fin ) = work%estimates(i)%SEs(:)
      end do
      work%fixed = work%fixed_tmp(:)
700   continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_estimate_SEs
    !##################################################################
    integer(kind=our_int) function advance_cell_fixed_part_estimate( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam_estimate), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, jj, var_old
      character(len=*), parameter :: &
           subname = "advance_cell_fixed_part_estimate"
      ! begin
      answer = RETURN_FAIL
      !
      if( work%begin_cycle_fixed ) then
         work%cell_fixed_part = 0
         do j = 1, work%nvar
            if( .not. work%fixed(j) ) cycle
            work%var(j) = 1
         end do
         work%begin_cycle_fixed = .false.
      else
         if( work%cycle_done_fixed ) goto 150
         do j = 1, work%nvar
            ! find the first fixed variable j that is not maxed out
            if( .not. work%fixed(j) ) cycle
            if( work%var(j) == work%nlev(j) ) cycle
            work%var(j) = work%var(j) + 1
            work%cell_fixed_part = work%cell_fixed_part + work%cumprod(j)
            ! go back and reset each previous fixed variable to 1
            do jj = 1, j - 1
               if( .not. work%fixed(jj) ) cycle
               var_old = work%var(jj)
               work%var(jj) = 1
               work%cell_fixed_part = work%cell_fixed_part + &
                    ( work%var(jj) - var_old ) * work%cumprod(jj)
            end do
            exit
         end do
      end if
      work%cell = 1 + work%cell_fixed_part + work%cell_random_part
      work%cycle_done_fixed = .true.
      do j = 1, work%nvar
         if( .not. work%fixed(j) ) cycle
         if( work%var(j) /= work%nlev(j) ) then
            work%cycle_done_fixed = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_cell_fixed_part_estimate
   !##################################################################
    integer(kind=our_int) function advance_cell_random_part_estimate( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam_estimate), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, jj, var_old
      character(len=*), parameter :: &
           subname = "advance_cell_random_part_estimate"
      ! begin
      answer = RETURN_FAIL
      if( work%begin_cycle_random ) then
         work%cell_random_part = 0
         do j = 1, work%nvar
            if( work%fixed(j) ) cycle
            work%var(j) = 1
         end do
         work%begin_cycle_random = .false.
      else
         if( work%cycle_done_random ) goto 150
         do j = 1, work%nvar
            ! find the first random variable j that is not maxed out
            if( work%fixed(j) ) cycle
            if( work%var(j) == work%nlev(j) ) cycle
            work%var(j) = work%var(j) + 1
            work%cell_random_part = work%cell_random_part + work%cumprod(j)
            ! go back and reset each previous random variable to 1
            do jj = 1, j - 1
               if( work%fixed(jj) ) cycle
               var_old = work%var(jj)
               work%var(jj) = 1
               work%cell_random_part = work%cell_random_part + &
                    ( work%var(jj) - var_old ) * work%cumprod(jj)
            end do
            exit
         end do
      end if
      work%cell = 1 + work%cell_fixed_part + work%cell_random_part
      work%cycle_done_random = .true.
      do j = 1, work%nvar
         if( work%fixed(j) ) cycle
         if( work%var(j) /= work%nlev(j) ) then
            work%cycle_done_random = .false.
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
150   call err_handle(err, 1, &
            comment = "Cycle is already done" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function advance_cell_random_part_estimate
   !##################################################################
   integer(kind=our_int) function run_cvam_predict_em( &
        model_type_int, method_int, &
        pred_data, pred_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        predict_var_info, &
        ! workspaces
        work, err, &
        ! outputs
        pred_mat &
        ) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: pred_data(:,:)
      integer(kind=our_int), intent(in) :: pred_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      integer(kind=our_int), intent(in) :: predict_var_info(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! outputs
      real(kind=our_dble), intent(out) :: pred_mat(:,:)
      ! declare locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_cvam_predict_em"
      ! begin
      answer = RETURN_FAIL
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( pred_data, pred_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_params_into_workspace( prob, beta, vhat_beta, work, err ) &
           == RETURN_FAIL ) goto 800
      work%nrow_prior_data = 0
      if( structural_zero_check( work, err, check_prob = .true. ) &
           == RETURN_FAIL ) goto 800
      if( create_predict_objects( &
           predict_var_info, work, err ) == RETURN_FAIL ) goto 800
      if( compute_predictions( work%prob, pred_mat, work, err ) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_predict_em
   !##################################################################
    integer(kind=our_int) function create_predict_objects( &
         predict_var_info, work, err ) result(answer)
      ! additional workspace objects
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: predict_var_info(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, status, posn, &
           nvar, cell, ncells
      integer(kind=our_int), allocatable :: data_matrix(:,:)
      character(len=*), parameter :: &
           subname = "create_predict_objects"
      ! begin
      answer = RETURN_FAIL
      ! set up predicts
      if( size( predict_var_info, 2 ) /= 2 ) goto 200
      nvar = size( predict_var_info, 1 )
      if( nvar == 0 ) goto 205
      work%predict%nvar = nvar
      !
      posn = 0
      allocate( work%predict%nlev(nvar), &
           work%predict%model_posn(nvar), &
           work%predict%fixed(nvar), &
           stat = status )
      if( status /= 0 ) goto 100
      do j = 1, work%predict%nvar
         work%predict%nlev(j) = predict_var_info(j,1)
         work%predict%model_posn(j) = predict_var_info(j,2)
         work%predict%fixed(j) = .false.
      end do
      !
      nvar = work%nvar
      allocate( work%predict%var_present(nvar), &
           work%predict%predict_posn(nvar), &
           stat = status )
      if( status /= 0 ) goto 100
      work%predict%var_present(:) = .false.
      work%predict%predict_posn(:) = 0
      work%predict%nvar_present = 0
      work%predict%nvar_absent = 0
      do j = 1, work%predict%nvar
         posn = work%predict%model_posn(j)
         work%predict%var_present(posn) = .true.
         work%predict%predict_posn(posn) = j
      end do
      work%predict%nvar_present = work%predict%nvar
      work%predict%nvar_absent = work%nvar - work%predict%nvar
      !
      nvar = work%predict%nvar
      allocate( work%predict%cumprod(nvar), stat=status )
      if( status /= 0 ) goto 100
      do j = 1, nvar
         if( j == 1 ) then
            work%predict%cumprod(j) = 1
         else
            work%predict%cumprod(j) = &
                 work%predict%cumprod(j-1) * &
                 work%predict%nlev(j-1)
         end if
      end do
      if( nvar == 0 ) then
         work%predict%ncells = 1
      else
         work%predict%ncells = &
              work%predict%cumprod( nvar ) * &
              work%predict%nlev( nvar )
      end if
      ncells = work%predict%ncells
      allocate( work%predict%prob(ncells), stat=status )
      if( status /= 0 ) goto 100
      work%predict%prob(:) = 0.D0
      allocate( work%predict%str_zero(ncells), stat=status )
      if( status /= 0 ) goto 100
      work%predict%str_zero(:) = .false.
      allocate( work%predict%var(nvar), stat=status )
      if( status /= 0 ) goto 100
      work%predict%var(:) = 0
      ! create a dummy data_matrix with a single row and missing
      ! values for all variables
      allocate( data_matrix(1, work%nvar), stat=status )
      if( status /= 0 ) goto 100
      data_matrix(1,:) = work%mvcode(:)
      ! create cells(:)
      allocate( work%predict%cells( work%ncells ), stat=status )
      if( status /= 0 ) goto 100
      work%predict%cells(:) = 0
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( 1, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         cell = 1
         do j = 1, work%nvar
            if( .not. work%predict%var_present(j) ) cycle
            posn = work%predict%predict_posn(j)
            cell = cell + ( work%var(j) - 1 ) * &
                 work%predict%cumprod(posn) 
         end do
         work%predict%cells( work%cell ) = cell
         if( work%cycle_done ) exit
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Argument predict_var_info has incorrect size" )
      goto 800
205   call err_handle(err, 1, &
            comment = "No variables to predict" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(data_matrix) ) deallocate( data_matrix ) 
    end function create_predict_objects
   !##################################################################
   integer(kind=our_int) function compute_predictions( &
         prob, pred_mat, work, err ) result( answer )
      ! computes the predictions
      ! given the probabilities for the 
      ! complete-data table in work%prob, which may have already been 
      ! conditioned on the variables fixed in the model
      implicit none
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(out) :: pred_mat(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: cfull, cpred, i
      real(kind=our_dble) :: sum, rtmp
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "compute_predictions"
      ! begin
      answer = RETURN_FAIL
      if( size(prob) /= work%ncells ) goto 50
      if( size( pred_mat, 1 ) /= work%nrow_input_data ) goto 60
      if( size( pred_mat, 2 ) /= work%predict%ncells ) goto 60
      !
      do i = 1, work%nrow_input_data
         if( work%input_data_freq_int(i) < 0 ) goto 70
         if( work%input_data_freq_int(i) == 0 ) then
            pred_mat(i,:) = 0.D0
            cycle
         end if
         work%predict%prob(:) = 0.D0
         work%begin_cycle = .true.
         work%cycle_done = .false.
         non_zero_cells_present = .false.
         sum = 0.D0
         do
            if( advance_to_next_cell( i, work%input_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            cfull = work%cell
            if( .not. work%str_zero(cfull) ) then
               non_zero_cells_present = .true.
               cpred = work%predict%cells( cfull )
               work%predict%prob(cpred) = work%predict%prob(cpred) &
                    + prob(cfull)
               sum = sum + prob(cfull)
            end if
            if( work%cycle_done ) exit
         end do
         if( non_zero_cells_present ) then
            if ( sum <= 0.D0 ) goto 100
         else
            goto 250
         end if
         rtmp = real( work%input_data_freq_int(i), our_dble)
         do cpred = 1, work%predict%ncells
            pred_mat(i,cpred) = rtmp * work%predict%prob(cpred) / sum
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
50    call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
60    call err_handle(err, 1, &
            comment = "Argument pred_mat has incorrect size" )
      goto 800
70    call err_handle(err, 1, &
            comment = "Negative frequency for row of pred_data" )
      call err_handle(err, 3, iobs=i)
      goto 800
100   call err_handle(err, 1, &
            comment = "Non-positive probability in denominator" )
      call err_handle(err, 1, &
            comment = "Bad row in pred_data" )
      call err_handle(err, 3, iobs=i)
      goto 800
250   call err_handle(err, 1, &
            comment = "Bad row in pred_data, with positive frequency" )
      call err_handle(err, 1, &
            comment = "assigned to structural zero cells" )
      call err_handle(err, 3, iobs=i)
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_predictions
    !##################################################################
   integer(kind=our_int) function run_cvam_impute_freq( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, synthetic_int, &
        ! workspaces
        work, err, &
        ! outputs
        result_mat, result_freq_int &
        ) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      integer(kind=our_int), intent(in) :: synthetic_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! outputs
      integer(kind=our_int), intent(out) :: result_mat(:,:)
      integer(kind=our_int), intent(out) :: result_freq_int(:)
      ! declare locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_cvam_impute_freq"
      ! begin
      answer = RETURN_FAIL
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( input_data, input_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_params_into_workspace( prob, beta, vhat_beta, work, err ) &
           == RETURN_FAIL ) goto 800
      if( fill_impute_result_mat( &
           result_mat, work, err ) == RETURN_FAIL ) goto 800
      work%nrow_prior_data = 0
      if( structural_zero_check( work, err, check_prob = .true. ) &
           == RETURN_FAIL ) goto 800
      work%exclude_all_na = .false.
      if( create_data_and_prior_use_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( fill_input_data_with_missing( synthetic_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( run_istep( work, err, use_flatten=.false., use_prior_data=.false., &
           use_input_data=.true. ) == RETURN_FAIL ) goto 800
      if( size( result_freq_int ) /= work%ncells ) goto 100
      result_freq_int(:) = work%freq_int(:)
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Argument result_freq_int has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_impute_freq
   !##################################################################
   integer(kind=our_int) function fill_impute_result_mat( &
         result_mat, work, err ) result( answer )
      implicit none
      integer(kind=our_int), intent(out) :: result_mat(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: status
      integer(kind=our_int), allocatable :: data_matrix(:,:)
      character(len=*), parameter :: &
           subname = "fill_impute_result_mat"
      ! begin
      answer = RETURN_FAIL
      if( size( result_mat, 1 ) /= work%ncells ) goto 60
      if( size( result_mat, 2 ) /= work%ncol_input_data ) goto 60
      ! create a dummy data_matrix with a single row and missing
      ! values for all variables
      allocate( data_matrix(1, work%nvar), stat=status )
      if( status /= 0 ) goto 100
      data_matrix(1,:) = work%mvcode(:)
      ! cycle through all cells
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( 1, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         result_mat( work%cell, : ) = work%var(:)
         if( work%cycle_done ) exit
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
60    call err_handle(err, 1, &
            comment = "Argument result_mat has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(data_matrix) ) deallocate( data_matrix ) 
    end function fill_impute_result_mat
    !##################################################################
    integer(kind=our_int) function run_istep( work, err, &
         use_flatten, use_prior_data, use_input_data, use_cell_means ) &
         result(answer)
      ! Performs a single I-step
      ! Distributes the frequencies in prior_data_freq and input_data_freq
      ! across the cells of the complete-data table ( work%freq_int ),
      ! given the current estimate of the cell probs ( work%prob )
      ! Final result is put into work%freq, which also includes flattening
      ! constant
      ! When finished, work%loglik contains the loglikelohood under
      ! the poisson model
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare optional args
      logical, intent(in), optional :: use_flatten
      logical, intent(in), optional :: use_prior_data
      logical, intent(in), optional :: use_input_data
      logical, intent(in), optional :: use_cell_means
      ! declare locals
      integer(kind=our_int) :: row, i, j
      real(kind=our_dble) :: ntotal, nprior, ndata, log_ntotal, sum
      logical :: use_flatten_local, use_prior_data_local, &
           use_input_data_local, use_cell_means_local
      character(len=*), parameter :: &
           subname = "run_istep"
      ! begin
      answer = RETURN_FAIL
      ! apply defaults
      if( present( use_flatten ) ) then
         use_flatten_local = use_flatten
      else
         use_flatten_local = .true.
      end if
      if( present( use_prior_data ) ) then
         use_prior_data_local = use_prior_data
      else
         use_prior_data_local = .true.
      end if
      if( present( use_input_data ) ) then
         use_input_data_local = use_input_data
      else
         use_input_data_local = .true.
      end if
      if( present( use_cell_means ) ) then
         if( use_cell_means .and. ( work%model_type /= "log-linear" ) ) goto 20
         use_cell_means_local = use_cell_means
      else
         use_cell_means_local = .false.
      end if
      ! initialize counts etc.
      work%freq(:) = 0.D0
      work%freq_int(:) = 0
      work%logprior = 0.D0
      work%loglik = 0.D0
      nprior = 0.D0
      ndata = 0.D0
      if( work%model_type == "log-linear" ) then
         ! poisson correction
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            work%loglik = work%loglik - work%mu(i)
         ! ridge contribution
         sum = 0.D0
         do j = 1, work%p
            sum = sum + work%beta(j)**2
         end do
         work%logprior = work%logprior - work%ridge * sum / 2.D0
         end do
      end if
      ! call flatten_table to accumulate logprior
      if( use_flatten_local ) then
         if( flatten_table( work, err, &
              use_cell_means = use_cell_means_local ) &
              == RETURN_FAIL ) goto 800
         nprior = nprior + work%flatten
      end if
      ! distribute prior data into freq_int
      if( use_prior_data_local ) then
         do row = 1, work%nrow_prior_data
            if( run_istep_single_row( row, work%prior_data, &
                 work%prior_data_freq_int, work%prior_data_use, &
                 work, err, prior = .true., &
                 use_cell_means = use_cell_means_local ) &
                 == RETURN_FAIL ) goto 800
            if( work%prior_data_use(row) ) nprior = nprior + &
                 work%prior_data_freq(row)
         end do
      end if
      ! distribute input data into freq_int
      if( use_input_data_local ) then
         do row = 1, work%nrow_input_data
            if( run_istep_single_row( row, work%input_data, &
                 work%input_data_freq_int, work%input_data_use, &
                 work, err, prior = .false., &
                 use_cell_means = use_cell_means_local ) &
                 == RETURN_FAIL ) goto 800
            if( work%input_data_use(row) ) ndata = ndata + &
                 work%input_data_freq(row)
         end do
      end if
      ! poisson corrections for saturated model
      if( work%model_type == "saturated" ) then
         ntotal = nprior + ndata
         if( ntotal > 0.D0 ) then
            log_ntotal = log(ntotal)
            work%loglik = work%loglik + ndata * log_ntotal - ntotal
            work%logprior = work%logprior + nprior * log_ntotal
         else
            ! nothing to do, because nprior and ndata are zero
         end if
      end if
      work%logP = work%logprior + work%loglik
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         work%freq(i) = work%freq(i) + real( work%freq_int(i), our_dble )
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_istep
    !##################################################################
    integer(kind=our_int) function run_istep_single_row( row, &
         data_matrix, data_freq_int, data_use, work, err, &
         prior, use_cell_means ) result( answer )
      ! performs the i-step for a single row of data_matrix
      ! (data_matrix is either work%input_data or work%prior_data)
      ! (data_freq_int is either work%input_data_freq_int
      ! or work%prior_data_freq_int)
      ! (data_use is either work%input_data_use or work%prior_data_use)
      ! Distributes the frequency in data_freq(row) 
      ! across the cells of the complete-data table ( work%counts ),
      ! given the current estimate of the cell probs ( work%prob ),
      ! and across the cells of each submodel table
      ! If prior = .true., accumulates work%logprior
      ! If prior = .false., accumulates work%loglik (default)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: row 
      integer(kind=our_int), intent(in) :: data_matrix(:,:)
      integer(kind=our_int), intent(in) :: data_freq_int(:)
      logical, intent(in) :: data_use(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: prior
      logical, intent(in), optional :: use_cell_means
      ! declare locals
      integer(kind=our_int) :: i, n_non_zero_cells, n_cells, f_remaining, &
           n_non_zero_cells_remaining, rb
      real(kind=our_dble) :: sum, ntmp, ptmp, rtmp
      logical :: prior_local, use_cell_means_local
      character(len=*), parameter :: &
           subname = "run_estep_single_row"
      ! begin
      answer = RETURN_FAIL
      if( present( prior ) ) then
         prior_local = prior
      else
         prior_local = .false.
      end if
      if( present( use_cell_means ) ) then
         if( use_cell_means .and. ( work%model_type /= "log-linear" ) ) goto 20
         use_cell_means_local = use_cell_means
      else
         use_cell_means_local = .false.
      end if
      if( ( row < 0 ) .or. ( row > size(data_matrix, 1) ) ) goto 100
      if( .not. data_use(row) ) goto 10 ! nothing to do
      if( data_freq_int(row) == 0 ) goto 10 ! nothing to do
      ! accumulate total probability for all cells of the 
      ! complete-data table that contribute to the given row
      ! of coarsened data
      work%begin_cycle = .true.
      work%cycle_done = .false.
      n_non_zero_cells = 0
      n_cells = 0
      sum = 0.D0
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         n_cells = n_cells + 1
         if( .not. work%str_zero(i) ) then
            n_non_zero_cells = n_non_zero_cells + 1
            if( use_cell_means_local ) then
               sum = sum + work%mu( i )
            else
               sum = sum + work%prob( i )
            end if
         end if
         if( work%cycle_done ) exit
      end do
      if( ( n_non_zero_cells > 0 ) .and. ( sum == 0.D0 ) ) goto 200
      if( ( n_non_zero_cells == 0 ) .and. &
           ( data_freq_int(row) > 0 ) ) goto 205
      if( sum < 0.D0 ) goto 210
      if( prior_local ) then
         if( n_non_zero_cells > 0 ) &
              work%logprior = work%logprior + &
              real(data_freq_int(row), our_dble) * log(sum)
      else
         if( n_non_zero_cells > 0 ) &
              work%loglik = work%loglik + &
              real(data_freq_int(row), our_dble) * log(sum)
      end if
      ! distribute the frequency proportionately across the cells
      ! of the complete-data table, and also the cells of the
      ! table for each sub-model
      work%begin_cycle = .true.
      work%cycle_done = .false.
      f_remaining = data_freq_int(row)
      n_non_zero_cells_remaining = n_non_zero_cells
      do
         if( advance_to_next_cell( row, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         i = work%cell
         if( .not. work%str_zero(i) ) n_non_zero_cells_remaining &
              = n_non_zero_cells_remaining - 1
         if( f_remaining < 0 ) goto 300
         if( work%cycle_done ) then
            if( n_non_zero_cells_remaining > 0 ) goto 310
            if( work%str_zero(i) .and. ( f_remaining > 0 ) ) goto 320
            work%freq_int(i) = work%freq_int(i) + f_remaining
            exit
         else
            if( n_non_zero_cells_remaining == 0 ) then
               ! this is the last non_zero cell
               work%freq_int(i) = work%freq_int(i) + f_remaining
               exit
            end if
         end if
         if( .not. work%str_zero(i) ) then
            if( use_cell_means_local ) then
               ptmp = work%mu(i) / sum
            else
               ptmp = work%prob(i) / sum
            end if
            ntmp = real( f_remaining, our_dble )
            if( rbinom_R( ntmp, ptmp, rtmp, err ) == RETURN_FAIL ) goto 800
            rb = int(rtmp, our_int)
            work%freq_int(i) = work%freq_int(i) + rb
            f_remaining = f_remaining - rb
            if( use_cell_means_local ) then
               sum = sum - work%mu(i)
            else
               sum = sum - work%prob(i)
            end if
            if( f_remaining == 0 ) exit
         end if
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument row out of bounds" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Bad row in input data, positive freq for zero cells" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Negative value for f_remaining" )
      goto 800
310   call err_handle(err, 1, &
            comment = "Something bad happened" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Something really bad happened" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_istep_single_row
    !##################################################################
   integer(kind=our_int) function run_cvam_impute_microdata( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, synthetic_int, &
        ! workspaces
        work, err, &
        ! outputs
        result_mat &
        ) result(answer)
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      integer(kind=our_int), intent(in) :: synthetic_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! outputs
      integer(kind=our_int), intent(out) :: result_mat(:,:)
      ! declare locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_cvam_impute_microdata"
      ! begin
      answer = RETURN_FAIL
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( input_data, input_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_params_into_workspace( prob, beta, vhat_beta, work, err ) &
           == RETURN_FAIL ) goto 800
      work%nrow_prior_data = 0
      if( structural_zero_check( work, err, check_prob = .true. ) &
           == RETURN_FAIL ) goto 800
      work%exclude_all_na = .false.
      if( create_data_and_prior_use_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( fill_input_data_with_missing( synthetic_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( run_istep_microdata( result_mat, work, err ) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_impute_microdata
   !##################################################################
    integer(kind=our_int) function fill_input_data_with_missing( &
         synthetic_int, work, err) result(answer)
      implicit none
      integer(kind=our_int), intent(in) :: synthetic_int
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: row, j
      character(len=*), parameter :: &
           subname = "fill_input_data_with_missing"
      ! begin
      answer = RETURN_FAIL
      if( synthetic_int == 1 ) then
         do j = 1, work%nvar
            if( work%fixed(j) ) cycle
            do row = 1, work%nrow_input_data
               work%input_data(row,j) = work%mvcode(j)
            end do
         end do
      else if( synthetic_int == 0 ) then
         ! do nothing
      else
         goto 100
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
100   call err_handle(err, 1, &
            comment = "Value of synthetic_int not recognized" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function fill_input_data_with_missing
    !##################################################################
    integer(kind=our_int) function run_istep_microdata( output_data, &
         work, err) result(answer)
      ! Performs a single I-step
      ! Distributes the frequencies in prior_data_freq and input_data_freq
      ! across the cells of the complete-data table ( work%freq_int ),
      ! given the current estimate of the cell probs ( work%prob )
      implicit none
      ! output
      integer(kind=our_int), intent(out) :: output_data(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: row
      integer(kind=our_int) :: i, n_non_zero_cells, n_cells, &
           n_non_zero_cells_remaining
      real(kind=our_dble) :: sum, ptmp, rtmp
      character(len=*), parameter :: &
           subname = "run_istep_microdata"
      ! begin
      answer = RETURN_FAIL
      if( size(output_data, 1) /= work%nrow_input_data ) goto 100
      if( size(output_data, 2) /= work%nvar ) goto 100
      !
      do row = 1, work%nrow_input_data
         if( .not. work%input_data_use(row) ) cycle
         work%begin_cycle = .true.
         work%cycle_done = .false.
         n_non_zero_cells = 0
         n_cells = 0
         sum = 0.D0
         do
            if( advance_to_next_cell( row, work%input_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            i = work%cell
            n_cells = n_cells + 1
            if( .not. work%str_zero(i) ) then
               n_non_zero_cells = n_non_zero_cells + 1
               sum = sum + work%prob( i )
            end if
            if( work%cycle_done ) exit
         end do
         if( ( n_non_zero_cells > 0 ) .and. ( sum == 0.D0 ) ) goto 200
         if( n_non_zero_cells == 0 ) goto 205
         if( sum < 0.D0 ) goto 210
         ! draw a single cell from available non-zero cells
         work%begin_cycle = .true.
         work%cycle_done = .false.
         n_non_zero_cells_remaining = n_non_zero_cells
         do
            if( advance_to_next_cell( row, work%input_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            i = work%cell
            if( .not. work%str_zero(i) ) n_non_zero_cells_remaining &
                 = n_non_zero_cells_remaining - 1
            if( work%cycle_done ) then
               if( n_non_zero_cells_remaining > 0 ) goto 310
               output_data(row,:) = work%var(:)
               exit
            else
               if( n_non_zero_cells_remaining == 0 ) then
                  ! this is the last non_zero cell
                  output_data(row,:) = work%var(:)
                  exit
               end if
            end if
            if( .not. work%str_zero(i) ) then
               ptmp = work%prob(i) / sum
               if( runif_R( rtmp, err ) == RETURN_FAIL ) goto 800
               if( rtmp <= ptmp ) then
                  output_data(row,:) = work%var(:)
                  exit
               end if
               sum = sum - work%prob(i)
            end if
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
100   call err_handle(err, 1, &
            comment = "Argument output_data has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Bad row in input data" )
      call err_handle(err, 1, &
            comment = "Observation found in structural zero cell" )
      call err_handle(err, 3, iobs=row)
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
310   call err_handle(err, 1, &
            comment = "Something bad happened" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_istep_microdata
    !##################################################################
   integer(kind=our_int) function run_cvam_lik( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        lik_var_info, &
        ! workspaces
        work, err, &
        ! outputs
        lik_values &
        ) result(answer)
      ! When calling this function, all values of input_data_freq_int
      ! should be zero. Otherwise it could mess up the structural
      ! zero check
      implicit none
      ! declare inputs
      integer(kind=our_int), intent(in) :: model_type_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: input_data(:,:)
      integer(kind=our_int), intent(in) :: input_data_freq_int(:)
      integer(kind=our_int), intent(in) :: n_levels_matrix(:,:)
      integer(kind=our_int), intent(in) :: packed_map(:)
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: offset(:)
      integer(kind=our_int), intent(in) :: str_zero_int(:)
      real(kind=our_dble), intent(in) :: prob(:)
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: vhat_beta(:,:)
      integer(kind=our_int), intent(in) :: lik_var_info(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! outputs
      real(kind=our_dble), optional, intent(out) :: lik_values(:)
      ! declare locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_cvam_lik"
      ! begin
      answer = RETURN_FAIL
      if( put_model_type_into_workspace( model_type_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_method_into_workspace( method_int, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_data_into_workspace( input_data, input_data_freq_int, &
           n_levels_matrix, packed_map, work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_model_into_workspace( model_matrix, offset, str_zero_int, &
           work, err ) == RETURN_FAIL ) goto 800
      if( create_model_fitting_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( put_params_into_workspace( prob, beta, vhat_beta, work, err ) &
           == RETURN_FAIL ) goto 800
      work%nrow_prior_data = 0
      if( structural_zero_check( work, err, check_prob = .true. ) &
           == RETURN_FAIL ) goto 800
      work%exclude_all_na = .false.
      if( create_data_and_prior_use_objects( work, err ) &
           == RETURN_FAIL ) goto 800
      if( create_lik_objects( &
           lik_var_info, work, err ) == RETURN_FAIL ) goto 800
      if( compute_lik_table( work%prob, work, err ) &
           == RETURN_FAIL ) goto 800
      if( compute_lik_values( lik_values, work, err ) &
           == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_cvam( work, err )
    end function run_cvam_lik
    !##################################################################
    integer(kind=our_int) function create_lik_objects( &
         lik_var_info, work, err ) result(answer)
      ! This routine creates a single estimate object for use in a
      ! likelihood calculation, but many of its components are never
      ! used or even allocated
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: lik_var_info(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, status, posn, &
           nvar, cell, ncells
      integer(kind=our_int), allocatable :: data_matrix(:,:)
      character(len=*), parameter :: &
           subname = "create_lik_objects"
      ! begin
      answer = RETURN_FAIL
      ! set up estimates
      if( size( lik_var_info, 2 ) /= 3 ) goto 200
      work%n_estimates = 1
      allocate( work%estimates(work%n_estimates), stat=status )
      if( status /= 0 ) goto 100
      nvar = size( lik_var_info, 1 )
      if( nvar == 0) goto 210
      i = 1
      !
      work%estimates(i)%estimate_number = i
      work%estimates(i)%nvar = nvar
      allocate( work%estimates(i)%nlev(nvar), &
           work%estimates(i)%model_posn(nvar), &
           work%estimates(i)%fixed(nvar), &
           stat = status )
      if( status /= 0 ) goto 100
      do j = 1, nvar
         work%estimates(i)%nlev(j) = lik_var_info(j,1)
         work%estimates(i)%model_posn(j) = lik_var_info(j,2)
         if( lik_var_info(j,3) == 1 ) then
            work%estimates(i)%fixed(j) = .false.
         else if( lik_var_info(j,3) == 2 ) then
            work%estimates(i)%fixed(j) = .true.
         else
            goto 320
         end if
      end do
      !
      nvar = work%nvar ! number of vars in full model
      allocate( work%estimates(i)%var_present(nvar), &
           work%estimates(i)%estimate_posn(nvar), &
           stat = status )
      if( status /= 0 ) goto 100
      work%estimates(i)%var_present(:) = .false.
      work%estimates(i)%estimate_posn(:) = 0
      work%estimates(i)%nvar_present = 0
      work%estimates(i)%nvar_absent = 0
      do j = 1, work%estimates(i)%nvar
         posn = work%estimates(i)%model_posn(j)
         work%estimates(i)%var_present(posn) = .true.
         work%estimates(i)%estimate_posn(posn) = j
      end do
      work%estimates(i)%nvar_present = work%estimates(i)%nvar
      work%estimates(i)%nvar_absent = work%nvar - work%estimates(i)%nvar
      !
      nvar = work%estimates(i)%nvar
      allocate( work%estimates(i)%cumprod(nvar), stat=status )
      if( status /= 0 ) goto 100
      do j = 1, nvar
         if( j == 1 ) then
            work%estimates(i)%cumprod(j) = 1
         else
            work%estimates(i)%cumprod(j) = &
                 work%estimates(i)%cumprod(j-1) * &
                 work%estimates(i)%nlev(j-1)
         end if
      end do
      if( nvar == 0 ) then
         work%estimates(i)%ncells = 1
      else
         work%estimates(i)%ncells = &
              work%estimates(i)%cumprod( nvar ) * &
              work%estimates(i)%nlev( nvar )
      end if
      ncells = work%estimates(i)%ncells
      allocate( work%estimates(i)%prob(ncells), stat=status )
      if( status /= 0 ) goto 100
      work%estimates(i)%prob(:) = 0.D0
      allocate( work%estimates(i)%str_zero(ncells), stat=status )
      if( status /= 0 ) goto 100
      work%estimates(i)%str_zero(:) = .false.
      allocate( work%estimates(i)%var(nvar), stat=status )
      if( status /= 0 ) goto 100
      work%estimates(i)%var(:) = 0
      ! 
      ! create a dummy data_matrix with a single row and missing
      ! values for all variables
      allocate( data_matrix(1, work%nvar), stat=status )
      if( status /= 0 ) goto 100
      data_matrix(1,:) = work%mvcode(:)
      ! create cells(:)
      allocate( work%estimates(i)%cells( work%ncells ), stat=status )
      if( status /= 0 ) goto 100
      work%estimates(i)%cells(:) = 0
      work%begin_cycle = .true.
      work%cycle_done = .false.
      do
         if( advance_to_next_cell( 1, data_matrix, work, err ) &
              == RETURN_FAIL ) goto 800
         cell = 1
         do j = 1, work%nvar
            if( .not. work%estimates(i)%var_present(j) ) cycle
            posn = work%estimates(i)%estimate_posn(j)
            cell = cell + ( work%var(j) - 1 ) * &
                 work%estimates(i)%cumprod(posn) 
         end do
         work%estimates(i)%cells( work%cell ) = cell
         if( work%cycle_done ) exit
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Argument lik_var_info has incorrect size" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Estimate has no variables" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Invalid value in lik_var_info, column 3" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(data_matrix) ) deallocate( data_matrix ) 
    end function create_lik_objects
    !##################################################################
   integer(kind=our_int) function compute_lik_table( &
         prob, work, err ) result( answer )
      ! computes the estimates
      ! given the probabilities for the 
      ! complete-data table in work%prob, which may have already been 
      ! conditioned on the variables fixed in the model, but not
      ! necessarily all variables that are fixed in the estimate
      implicit none
      real(kind=our_dble), intent(in) :: prob(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, cfull, cest
      real(kind=our_dble) :: sum
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "compute_lik_table"
      ! begin
      answer = RETURN_FAIL
      if( size(prob) /= work%ncells ) goto 50
      do i = 1, work%n_estimates
         ! marginalize over the variables not present
         work%estimates(i)%prob(:) = 0.D0
         work%estimates(i)%str_zero(:) = .true.
         do cfull = 1, work%ncells
            if( work%str_zero(cfull) ) cycle
            cest = work%estimates(i)%cells( cfull )
            work%estimates(i)%str_zero( cest ) = .false.
            work%estimates(i)%prob( cest ) = &
                 work%estimates(i)%prob( cest ) + prob( cfull )
         end do
         ! condition on the variables that are fixed in the estimate,
         ! which must include any variables fixed in the model, and
         ! possibly some additional ones
         work%estimates(i)%begin_cycle_fixed = .true.
         work%estimates(i)%cycle_done_fixed = .false.
         do
            if( advance_cell_fixed_part_estimate( &
                 work%estimates(i), err ) &
                 == RETURN_FAIL ) goto 800
            non_zero_cells_present = .false.
            sum = 0.D0
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  non_zero_cells_present = .true.
                  sum = sum + work%estimates(i)%prob( cest )
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            if( non_zero_cells_present .and. ( sum == 0.D0 ) ) goto 100
            work%estimates(i)%begin_cycle_random = .true.
            work%estimates(i)%cycle_done_random = .false.
            do
               if( advance_cell_random_part_estimate( &
                    work%estimates(i), err ) &
                    == RETURN_FAIL ) goto 800
               cest = work%estimates(i)%cell
               if( .not. work%estimates(i)%str_zero( cest ) ) then
                  work%estimates(i)%prob( cest ) = &
                       work%estimates(i)%prob( cest ) / sum
               end if
               if( work%estimates(i)%cycle_done_random ) exit
            end do
            if( work%estimates(i)%cycle_done_fixed ) exit
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
50    call err_handle(err, 1, &
            comment = "Argument prob has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Non-positive probability in denominator" )
      call err_handle(err, 12, iest = i )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_lik_table
    !##################################################################
   integer(kind=our_int) function compute_lik_values( &
         lik_values, work, err ) result( answer )
      ! this computation assumed there are no missing or 
      ! coarsened values in any of the variables being conditioned on
      implicit none
      real(kind=our_dble), intent(out) :: lik_values(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, row, cfull, cest, n_non_zero_cells, n_cells
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_lik_values"
      ! begin
      answer = RETURN_FAIL
      if( size(lik_values) /= work%nrow_input_data ) goto 50
      i = 1
      do row = 1, work%nrow_input_data
         if( .not. work%input_data_use(row) ) cycle ! nothing to do
         work%begin_cycle = .true.
         work%cycle_done = .false.
         n_non_zero_cells = 0
         n_cells = 0
         sum = 0.D0
         do
            if( advance_to_next_cell( row, work%input_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            cfull = work%cell
            cest = work%estimates(i)%cells( cfull )
            n_cells = n_cells + 1
            if( .not. work%str_zero(cfull) ) then
               n_non_zero_cells = n_non_zero_cells + 1
               sum = sum + work%estimates(i)%prob( cest )
            end if
            if( work%cycle_done ) exit
         end do
         lik_values(row) = sum
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
50    call err_handle(err, 1, &
            comment = "Argument lik_values has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_lik_values
    !##################################################################
    integer(kind=our_int) function create_mcmc_objects( &
         beta_series, prob_series, logp_series, imputed_freq_int, &
         packed_estimates_series, work, err ) result(answer)
      ! checks dimensions of mcmc output objects and puts relevant
      ! dimensions into workspace
      implicit none
      ! outputs (we only query the sizes in this function)
      real(kind=our_dble), intent(out) :: beta_series(:,:)
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: status
      character(len=*), parameter :: subname = "create_mcmc_objects"
      ! begin
      answer = RETURN_FAIL
      if( .not. ( work%method == "MCMC" .or. &
          work%method == "approxBayes" ) ) goto 10    ! nothing to do
      if( work%method == "MCMC" ) then      
         if( mod( work%iter_mcmc, work%thin_mcmc ) /= 0  ) goto 50
         work%series_length = work%iter_mcmc / work%thin_mcmc
         if( work%impute_every == 0 ) then
            work%n_impute = 0
         else
            work%n_impute = floor( real(work%iter_mcmc, our_dble) / &
                 real(work%impute_every, our_dble) )
         end if
      else
         work%series_length = work%iter_approx_bayes
         if( work%impute_approx_bayes ) then
            work%n_impute = work%iter_approx_bayes
         else
            work%n_impute = 0
         end if
      end if
      if( work%save_prob_series ) then
         work%series_length_prob = work%series_length
      else
         work%series_length_prob = 0
      end if
      if( work%model_type == "log-linear" ) then
         if( size( beta_series, 1 ) /= work%series_length ) goto 200
         if( size( beta_series, 2 ) /= work%p ) goto 200
      end if
      if( work%save_prob_series ) then
         if( size( prob_series, 1 ) /= work%series_length ) goto 210
      else
         if( size( prob_series, 1 ) /= 0 ) goto 210
      end if
      if( size( prob_series, 2 ) /= work%ncells ) goto 210
      if( size( logp_series ) /= work%series_length ) goto 220
      if( size( packed_estimates_series, 1 ) /= work%series_length ) &
           goto 225
      if( size( packed_estimates_series, 2 ) /= work%n_packed_estimates ) &
           goto 225
      if( size( imputed_freq_int, 1 ) /= work%ncells ) goto 230
      if( size( imputed_freq_int, 2 ) /= work%n_impute ) goto 230
      ! allocate workspace items
      if( work%model_type == "log-linear" ) then
         allocate( work%mh_ratios_beta( &
              work%iter_mcmc + work%burn_mcmc ), &
              work%mh_accept_beta( &
              work%iter_mcmc + work%burn_mcmc ), stat=status )
         if( status /= 0 ) goto 100
         allocate( work%beta_can( work%p ), work%beta_center( work%p ), &
              work%beta_scale( work%p, work%p ), &
              work%beta_scale_inv( work%p, work%p ), &
              work%beta_scale_inv_sqrt( work%p, work%p ), &
              work%beta_scale_sqrt( work%p, work%p ), &
              work%beta_mean( work%p ), &
              work%beta_cov_mat( work%p, work%p ), stat=status )
         if( status /= 0 ) goto 100
      end if
      allocate( work%prob_mean( work%ncells ), stat=status )
      if( status /= 0 ) goto 100
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
50    call err_handle(err, 1, &
            comment = "iter_mcmc is not a multiple of thin_mcmc" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate workspace array component" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Argument beta_series has incorrect size" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Argument prob_series has incorrect size" )
      goto 800
220   call err_handle(err, 1, &
            comment = "Argument prob_series has incorrect size" )
      goto 800
225   call err_handle(err, 1, &
            comment = "Argument packed_estimates_series has incorrect size" )
      goto 800
230   call err_handle(err, 1, &
            comment = "Argument imputed_freq_int has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function create_mcmc_objects
    !##################################################################
    integer(kind=our_int) function run_da_log_linear( &
         beta_series, prob_series, logp_series, imputed_freq_int, &
         packed_estimates_series, &
         n_iter_actual, n_sample_actual, n_imp_actual, &
         work, err ) result(answer)
      implicit none
      ! outputs 
      real(kind=our_dble), intent(out) :: beta_series(:,:)
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter, thin_step, imp_step
      real(kind=our_dble) :: rtmp, logP_tmp, loglik_tmp, logprior_tmp
      logical :: aborted
      character(len=*), parameter :: subname = "run_da_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( work%method /= "MCMC" ) goto 30 
      ! initialize iteration counters
      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%prob_mean(:) = 0.D0
      work%freq_mean(:) = 0.D0
      if( work%n_estimates > 0 ) work%packed_estimates_mean = 0.D0
      if( size(beta_series) > 0 ) beta_series(:,:) = 0.D0
      if( size(prob_series) > 0 ) prob_series(:,:) = 0.D0
      if( size(logp_series) > 0 ) logp_series(:) = 0.D0
      if( size( imputed_freq_int) > 0 ) imputed_freq_int(:,:) = 0
      if( size(packed_estimates_series) > 0 ) &
           packed_estimates_series(:,:) = 0.D0
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      aborted = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         if( iter == 1 ) then
            if( compute_mu_from_beta( work%beta, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( compute_estimates( work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
         else
            ! mu, logmu and prob are already consistent with beta
         end if
         if( run_istep( work, err, use_flatten = .true., &
              use_prior_data = .true., &
              use_input_data = .true., use_cell_means = .true.) &
              == RETURN_FAIL ) goto 10
         if( work%iter == 1 ) work%start_logP = work%logP
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         if( run_mh_step_beta_da( work, err ) == RETURN_FAIL ) goto 10
         if( work%mh_accept_beta( work%iter ) ) then
            ! mu, logmu are already up-to-date, but prob is not
            if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( compute_estimates( work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
         end if
         ! update counters and running estimates
         if( work%mh_accept_beta( work%iter ) ) then
            ! candidate was accepted
            work%beta_accept_count = work%beta_accept_count + 1
            work%beta_current_reject_run = 0
            rtmp = 1.D0
         else
            ! candidate was rejected
            work%beta_current_reject_run = &
                 work%beta_current_reject_run + 1
            rtmp = 0.D0
         end if
         work%beta_accept_rate = work%beta_accept_rate + &
              ( rtmp - work%beta_accept_rate ) / &
              real( work%iter, our_int )
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 200
            if( work%store_count > size(beta_series, 1) ) goto 200
            if( work%store_count > size(logp_series) ) goto 200
            beta_series( work%store_count, : ) = work%beta(:)
            logp_series( work%store_count ) = work%logP
            if( work%save_prob_series ) then
               if( work%store_count > size(prob_series, 1) ) goto 200
               prob_series( work%store_count, : ) = work%prob(:)
            end if
            if( work%n_estimates > 0 ) then
               if( work%store_count > size(packed_estimates_series, 1) ) &
                    goto 200
               packed_estimates_series( work%store_count, : ) = &
                    work%packed_estimates(:)
            end if
         end if
         ! generate imputation, and put the result in imputed_freq_int
         work%freq_int_tmp(:) = work%freq_int(:)
         work%freq_tmp(:) = work%freq(:)
         logP_tmp = work%logP
         loglik_tmp = work%loglik
         logprior_tmp = work%logprior
         work%input_data_use_tmp(:) = work%input_data_use(:)
         work%input_data_use(:) = .true.
         if( run_istep( work, err, use_flatten = .false., &
              use_prior_data = .false., &
              use_input_data = .true. ) == RETURN_FAIL ) goto 10
         if( work%imp_this_iter ) then
            if( ( work%imp_count < 1 ) .or. &
                 ( work%imp_count > size(imputed_freq_int, 2) ) ) goto 200
            imputed_freq_int(:, work%imp_count) = work%freq_int(:)
         end if
         if( update_freq_mean( work, err ) == RETURN_FAIL ) goto 800
         work%freq_int(:) = work%freq_int_tmp(:)
         work%freq(:) = work%freq_tmp(:)
         work%logP = logP_tmp
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%input_data_use(:) = work%input_data_use_tmp(:)
         ! update the running means
         if( update_running_means( work, err ) == RETURN_FAIL ) goto 800
         if( work%beta_current_reject_run > work%stuck_limit ) then
            call err_handle(err, 1, &
                 comment = "Metropolis-Hastings got stuck" )
            goto 10
         end if
         !#######################
         aborted = .false.
      end do
      !### end main iteration
10    continue
      if( aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = work%iter )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_da_log_linear
    !##################################################################
    integer(kind=our_int) function run_mh_step_beta_da( work, err ) &
         result( answer )
      ! performs one step of metropolis-hastings on beta, conditioning
      ! on the frequencies in work%freq
      ! uses the augmented-data poisson likelihood
      implicit none
      ! args
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: rtmp, log_mh_ratio, mh_ratio, runif
      logical :: accept
      character(len=*), parameter :: &
           subname = "run_mh_step_beta_da"
      ! begin
      answer = RETURN_FAIL
      ! store mu and logmu to restore later if step is rejected
      work%mu_tmp(:) = work%mu(:)
      work%logmu_tmp(:) = work%logmu(:)
      !
      if( compute_logP_score_hessian_beta_da( work%beta, work, err ) &
           == RETURN_FAIL ) goto 800
      log_mh_ratio = - work%logP_mstep
      if( compute_center_and_scale_beta_da( work%beta, &
           work%df_da, work%step_size_da, work%scale_fac_da, &
           work, err ) == RETURN_FAIL ) goto 800
      if( draw_candidate_beta( work%df_da, &
           work, err ) == RETURN_FAIL ) goto 800
      if( compute_log_proposal_beta( work%beta_can, work%df_da, &
           rtmp, work, err ) == RETURN_FAIL ) goto 800
      log_mh_ratio = log_mh_ratio - rtmp
      if( compute_mu_from_beta( work%beta_can, work, err ) &
           == RETURN_FAIL ) goto 800
      if( compute_logP_score_hessian_beta_da( work%beta_can, work, err ) &
           == RETURN_FAIL ) goto 800
      log_mh_ratio = log_mh_ratio + work%logP_mstep
      if( compute_center_and_scale_beta_da( work%beta_can, &
           work%df_da, work%step_size_da, work%scale_fac_da, &
           work, err ) == RETURN_FAIL ) goto 800
      if( compute_log_proposal_beta( work%beta, work%df_da, &
           rtmp, work, err ) == RETURN_FAIL ) goto 800
      log_mh_ratio = log_mh_ratio + rtmp
      ! to prevent over/underflow
      if( log_mh_ratio > log_huge ) then
         mh_ratio = huge(0.D0)
      else if( log_mh_ratio < log_tiny ) then
         mh_ratio = 0.D0
      else
         mh_ratio = exp( log_mh_ratio )
      end if
      ! compare to uniform variate
      if( runif_R( runif, err ) == RETURN_FAIL ) goto 800
      if( runif <= mh_ratio ) then
         accept = .true.
         work%beta(:) = work%beta_can(:)
      else
         accept = .false.
         work%mu(:) = work%mu_tmp(:)
         work%logmu(:) = work%logmu_tmp(:)
      end if
      work%mh_ratios_beta( work%iter ) = mh_ratio
      work%mh_accept_beta( work%iter ) = accept
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_mh_step_beta_da
    !##################################################################
    integer(kind=our_int) function compute_logP_score_hessian_beta_da( &
        beta, work, err ) result(answer)
      ! computes logP and first two derivatives at beta
      ! under the Poisson model
      ! given the augmented-data frequencies stored in work%freq
      ! assumes that work%mu and work%logmu are consistent with beta
      ! results are stored in logP_mstep, score_mstep, hess_mstep
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_logP_score_hessian_beta_da"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      work%logP_mstep = 0.D0
      work%score_mstep(:) = 0.D0
      work%hess_mstep(:,:) = 0.D0
      ! ridge contributions
      sum = 0.D0
      do j = 1, work%p
         sum = sum + beta(j)**2
      end do
      work%logP_mstep = work%logP_mstep - work%ridge * sum / 2.D0
      do j = 1, work%p
         work%score_mstep(j) = work%score_mstep(j) &
              - work%ridge * beta(j)
         work%hess_mstep(j,j) = work%hess_mstep(j,j) - work%ridge
      end do
      ! apply Poisson correction
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         work%logP_mstep = work%logP_mstep - work%mu(i)
      end do
      do i = 1, work%n
         if( work%str_zero(i) ) cycle
         work%logP_mstep = work%logP_mstep + work%freq(i) * work%logmu(i)
         do j = 1, work%p
            work%score_mstep(j) = work%score_mstep(j) + ( work%freq(i) &
                 - work%mu(i) ) * work%model_matrix(i,j)
            do k = 1, j
               work%hess_mstep(j,k) = work%hess_mstep(j,k) - work%mu(i) * &
                    work%model_matrix(i,j) * work%model_matrix(i,k)
            end do
         end do
      end do
      do j = 1, work%p
         do k = 1, j
            work%hess_mstep(k,j) = work%hess_mstep(j,k)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_logP_score_hessian_beta_da
    !##################################################################
    integer(kind=our_int) function compute_center_and_scale_beta_da( &
        beta, df, step_size, scale_fac, work, err ) result(answer)
      ! computes beta_center and beta_scale for MH proposal given the first
      ! two derivatives in score_mstep and hess_mstep
      ! also computes beta_scale_inv,beta_ scale_inv_sqrt, and
      ! beta_scale_sqrt
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(in) :: step_size
      real(kind=our_dble), intent(in) :: scale_fac
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum
      character(len=*), parameter :: &
           subname = "compute_center_and_scale_beta_da"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( size(beta) /= work%p ) goto 30
      rtmp = - ( df / &
           ( df + real( work%p, our_dble) ) ) &
                 * ( 1.D0 / ( scale_fac**2 ) )
      work%beta_scale_inv(:,:) = rtmp * work%hess_mstep(:,:)
      work%beta_scale_inv_sqrt(:,:) = work%beta_scale_inv(:,:)
      if( cholesky_in_place( work%beta_scale_inv_sqrt, err ) &
           == RETURN_FAIL ) goto 200
      work%beta_scale_sqrt(:,:) = work%beta_scale_inv_sqrt(:,:)
      if( invert_lower(work%beta_scale_sqrt, err ) == RETURN_FAIL ) goto 200
      if( premult_lower_by_transpose( work%beta_scale_sqrt, &
           work%beta_scale, err) == RETURN_FAIL ) goto 800
      ! put inverse of negative hessian into wkppA
      work%wkppA(:,:) = - rtmp * work%beta_scale(:,:)
      ! find center of proposal
      do j = 1, work%p
         sum = 0.D0
         do k = 1, work%p
            sum = sum + work%wkppA(j,k) * work%score_mstep(k)
         end do
         work%beta_center(j) = beta(j) + step_size * sum
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Hessian not neg-def" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_center_and_scale_beta_da
    !##################################################################
    integer(kind=our_int) function draw_candidate_beta( df, &
         work, err ) result(answer)
      ! Draws candidate beta from multivariate t proposal, storing
      ! the result in mmw%beta_can 
      ! Depends on beta_center and beta_scale, so should be run after
      ! compute_center_and_scale_beta_da
      implicit none
      real(kind=our_dble), intent(in) :: df
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum, rnorm, rchisq
      character(len=*), parameter :: &
           subname = "draw_candidate_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      if( df <= 0.D0 ) goto 30
      !####
      ! draw independent t variates, put into wkpA
      if( rchisq_R( df, rchisq, err ) == RETURN_FAIL ) goto 800
      rtmp = sqrt( df / rchisq )
      do j = 1, work%p
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkpA(j) = rnorm * rtmp
      end do
      ! premultiply by t(scale_sqrt), put into wkpB
      do j = 1, work%p
         sum = 0.D0
         do k = j, work%p
            sum = sum + work%beta_scale_sqrt(k,j) * work%wkpA(k)
         end do
         work%wkpB(j) = sum
      end do
      ! add center to obtain candidate
      do j = 1, work%p
         work%beta_can(j) = work%wkpB(j) + work%beta_center(j)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_candidate_beta
    !##################################################################
    integer(kind=our_int) function compute_log_proposal_beta( beta, df, &
         ans, work, err ) result(answer)
      ! computes log-proposal density at beta, assuming that the
      ! proposal center and scale are in beta_center and beta_scale
      ! the result in mmw%beta_can 
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_log_proposal_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      if( df <= 0.D0 ) goto 30
      if( size(beta) /= work%p ) goto 40
      !####
      ! wkqB = t( beta - center) %*% (lower tri scale_inv_sqrt)
      work%wkpA(:) = beta(:) - work%beta_center(:)
      do j = 1, work%p
         sum = 0.D0
         do k = j, work%p
            sum = sum + work%wkpA(k) * work%beta_scale_inv_sqrt(k,j)
         end do
         work%wkpB(j) = sum
      end do
      ! sum = t( beta - center) %*% scale_inv %*% ( beta - center )
      sum = 0.D0
      do j = 1, work%p
         sum = sum + work%wkpB(j)**2
      end do
      ans = - (( df + real(work%p,our_dble) ) &
           / 2.D0 ) * log( 1.D0 + sum / df )
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_log_proposal_beta
    !##################################################################
    integer(kind=our_int) function update_running_means( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j, k 
      real(kind=our_dble) :: rtmp, nfloat
      character(len=*), parameter :: &
           subname = "update_running_means"
      ! begin
      answer = RETURN_FAIL
      if( work%iter_past_burn_in <= 0 ) goto 10 ! nothing to do
      nfloat = real( work%iter_past_burn_in, our_dble )
      if( nfloat == 0.D0 ) goto 200
      ! update prob_mean
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         rtmp = work%prob(i) - work%prob_mean(i)
         work%prob_mean(i) = work%prob_mean(i) + rtmp / nfloat
      end do
      !
      if( work%model_type == "log-linear" ) then
         ! update beta_mesn and beta_cov_mat
         do j = 1, work%p
            work%wkpA(j) = work%beta(j) - work%beta_mean(j)
         end do
         do j = 1, work%p
            do k = j, work%p
               work%wkppA(j,k) = work%wkpA(j) * work%wkpA(k)
            end do
         end do
         do j = 1, work%p
            work%beta_mean(j) = work%beta_mean(j) + work%wkpA(j) / nfloat
            do k = j, work%p
               work%beta_cov_mat(j,k) = &
                    ( nfloat - 1.D0 ) * work%beta_cov_mat(j,k) + &
                    ( ( nfloat - 1.D0 ) / nfloat ) * work%wkppA(j,k)
               work%beta_cov_mat(j,k) = work%beta_cov_mat(j,k) / nfloat
               work%beta_cov_mat(k,j) = work%beta_cov_mat(j,k)
            end do
         end do
      end if
      ! update packed_estimates_mean
      if( work%n_estimates > 0 ) then
         do i = 1, work%n_packed_estimates
            rtmp = work%packed_estimates(i) - &
                 work%packed_estimates_mean(i)
            work%packed_estimates_mean(i) = &
                 work%packed_estimates_mean(i) + rtmp / nfloat
         end do
      end if
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Attempted division by zero")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function update_running_means
    !##################################################################
    integer(kind=our_int) function run_da_saturated( &
         prob_series, logp_series, imputed_freq_int, &
         packed_estimates_series, &
         n_iter_actual, n_sample_actual, n_imp_actual, &
         work, err ) result(answer)
      implicit none
      ! outputs 
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter, thin_step, imp_step, i
      real(kind=our_dble) :: rtmp, logP_tmp, loglik_tmp, logprior_tmp
      logical :: aborted
      character(len=*), parameter :: subname = "run_da_saturated"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "saturated" ) goto 20 
      if( work%method /= "MCMC" ) goto 30 
      ! initialize iteration counters      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      work%prob_mean(:) = 0.D0
      work%freq_mean(:) = 0.D0
      if( size(prob_series) > 0 ) prob_series(:,:) = 0.D0
      if( size(logp_series) > 0 ) logp_series(:) = 0.D0
      if( size( imputed_freq_int) > 0 ) imputed_freq_int(:,:) = 0
      if( work%n_estimates > 0 ) work%packed_estimates_mean = 0.D0
      aborted = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         if( run_istep( work, err, use_flatten = .true., &
              use_prior_data = .true., &
              use_input_data = .true., use_cell_means = .false.) &
              == RETURN_FAIL ) goto 10
         if( work%iter == 1 ) work%start_logP = work%logP
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            if( work%freq(i) + 1.D0 <= 0.D0 ) then
               call err_handle(err, 1, &
                    comment = "Non-positive frequency encountered" )
               call err_handle(err, 15, icell = i )
               call err_handle(err, 1, &
                    comment = "Posterior distribution is not proper;" )
               call err_handle(err, 1, &
                    comment = "consider using a flattening prior" )
               goto 10
            end if
            if( rgamma_R( work%freq(i) + 1.D0, 1.D0, rtmp, err ) &
                 == RETURN_FAIL ) goto 800
            work%prob_new(i) = rtmp
         end do
         if( normalize_prob( work%prob_new, work%prob, work, err ) &
              == RETURN_FAIL ) goto 10
         if( compute_estimates( work%prob, work, err ) &
              == RETURN_FAIL ) goto 10
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 200
            if( work%store_count > size(logp_series) ) goto 200
            logp_series( work%store_count ) = work%logP
            if( work%save_prob_series ) then
               if( work%store_count > size(prob_series, 1) ) goto 200
               prob_series( work%store_count, : ) = work%prob(:)
            end if
            if( work%n_estimates > 0 ) then
               if( work%store_count > size(packed_estimates_series, 1) ) &
                    goto 200
               packed_estimates_series( work%store_count, : ) = &
                    work%packed_estimates(:)
            end if
         end if
         ! generate imputation, and put the result in imputed_freq_int
         work%freq_int_tmp(:) = work%freq_int(:)
         work%freq_tmp(:) = work%freq(:)
         logP_tmp = work%logP
         loglik_tmp = work%loglik
         logprior_tmp = work%logprior
         work%input_data_use_tmp(:) = work%input_data_use(:)
         work%input_data_use(:) = .true.
         if( run_istep( work, err, use_flatten = .false., &
              use_prior_data = .false., &
              use_input_data = .true. ) == RETURN_FAIL ) goto 10
         if( work%imp_this_iter ) then
            if( ( work%imp_count < 1 ) .or. &
                 ( work%imp_count > size(imputed_freq_int, 2) ) ) goto 200
            imputed_freq_int(:, work%imp_count) = work%freq_int(:)
         end if
         if( update_freq_mean( work, err ) == RETURN_FAIL ) goto 800
         work%freq_int(:) = work%freq_int_tmp(:)
         work%freq(:) = work%freq_tmp(:)
         work%logP = logP_tmp
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%input_data_use(:) = work%input_data_use_tmp(:)
         ! update the running means
         if( update_running_means( work, err ) == RETURN_FAIL ) goto 800
         aborted = .false.
      end do
      !### end main iteration
10    continue
      if( aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = work%iter )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error trap
20    call err_handle(err, 1, &
            comment = "The model is not saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC")
      goto 800
200   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_da_saturated
    !##################################################################
    integer(kind=our_int) function update_freq_mean( &
         work, err ) result(answer)
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: rtmp, nfloat
      character(len=*), parameter :: &
           subname = "update_freq_mean"
      ! begin
      answer = RETURN_FAIL
      if( work%iter_past_burn_in <= 0 ) goto 10 ! nothing to do
      nfloat = real( work%iter_past_burn_in, our_dble )
      if( nfloat == 0.D0 ) goto 200
      ! update freq_mean
      do i = 1, work%ncells
         if( work%str_zero(i) ) cycle
         rtmp = real( work%freq_int(i), our_dble ) - work%freq_mean(i)
         work%freq_mean(i) = work%freq_mean(i) + rtmp / nfloat
      end do
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
           comment = "Attempted division by zero")
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function update_freq_mean
    !##################################################################
    integer(kind=our_int) function run_rwm_log_linear( &
         beta_series, prob_series, logp_series, imputed_freq_int, &
         packed_estimates_series, &
         n_iter_actual, n_sample_actual, n_imp_actual, &
         work, err ) result(answer)
      implicit none
      ! outputs 
      real(kind=our_dble), intent(out) :: beta_series(:,:)
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter, thin_step, imp_step
      real(kind=our_dble) :: rtmp, logP_tmp, loglik_tmp, logprior_tmp
      logical :: aborted
      character(len=*), parameter :: subname = "run_rwm_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( work%method /= "MCMC" ) goto 30 
      ! initialize iteration counters
      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%prob_mean(:) = 0.D0
      work%freq_mean(:) = 0.D0
      if( work%n_estimates > 0 ) work%packed_estimates_mean = 0.D0
      if( size(beta_series) > 0 ) beta_series(:,:) = 0.D0
      if( size(prob_series) > 0 ) prob_series(:,:) = 0.D0
      if( size(logp_series) > 0 ) logp_series(:) = 0.D0
      if( size( imputed_freq_int) > 0 ) imputed_freq_int(:,:) = 0
      if( size(packed_estimates_series) > 0 ) &
           packed_estimates_series(:,:) = 0.D0
      if( compute_scale_rwm( work%df_rwm, work%scale_fac_rwm, work, err ) &
           == RETURN_FAIL ) goto 800
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      aborted = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         if( work%iter == 1 ) then
            if( compute_mu_from_beta( work%beta, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( compute_estimates( work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
         else
            ! mu, logmu and prob are already consistent with beta
         end if
         if( work%iter == 1 ) then
            if( compute_loglik_logprior( work%beta, work, err, &
                 use_cell_means = .true., logprior = logprior_tmp, &
                 loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
            work%loglik = loglik_tmp
            work%logprior = logprior_tmp
            work%logP = logprior_tmp + loglik_tmp
            work%start_logP = work%logP
         end if
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         if( run_rwm_step_beta( work, err ) == RETURN_FAIL ) goto 10
         if( work%mh_accept_beta( work%iter ) ) then
            ! mu, logmu are already up-to-date, but prob is not
            if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( compute_estimates( work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
         end if
         ! update counters and running estimates
         if( work%mh_accept_beta( work%iter ) ) then
            ! candidate was accepted
            work%beta_accept_count = work%beta_accept_count + 1
            work%beta_current_reject_run = 0
            rtmp = 1.D0
         else
            ! candidate was rejected
            work%beta_current_reject_run = &
                 work%beta_current_reject_run + 1
            rtmp = 0.D0
         end if
         work%beta_accept_rate = work%beta_accept_rate + &
              ( rtmp - work%beta_accept_rate ) / &
              real( work%iter, our_int )
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 200
            if( work%store_count > size(beta_series, 1) ) goto 200
            if( work%store_count > size(logp_series) ) goto 200
            beta_series( work%store_count, : ) = work%beta(:)
            logp_series( work%store_count ) = work%logP
            if( work%save_prob_series ) then
               if( work%store_count > size(prob_series, 1) ) goto 200
               prob_series( work%store_count, : ) = work%prob(:)
            end if
            if( work%n_estimates > 0 ) then
               if( work%store_count > size(packed_estimates_series, 1) ) &
                    goto 200
               packed_estimates_series( work%store_count, : ) = &
                    work%packed_estimates(:)
            end if
         end if
         ! generate imputation, and put the result in imputed_freq_int
         work%freq_int_tmp(:) = work%freq_int(:)
         work%freq_tmp(:) = work%freq(:)
         logP_tmp = work%logP
         loglik_tmp = work%loglik
         logprior_tmp = work%logprior
         work%input_data_use_tmp(:) = work%input_data_use(:)
         work%input_data_use(:) = .true.
         if( run_istep( work, err, use_flatten = .false., &
              use_prior_data = .false., &
              use_input_data = .true. ) == RETURN_FAIL ) goto 10
         if( work%imp_this_iter ) then
            if( ( work%imp_count < 1 ) .or. &
                 ( work%imp_count > size(imputed_freq_int, 2) ) ) goto 200
            imputed_freq_int(:, work%imp_count) = work%freq_int(:)
         end if
         if( update_freq_mean( work, err ) == RETURN_FAIL ) goto 800
         work%freq_int(:) = work%freq_int_tmp(:)
         work%freq(:) = work%freq_tmp(:)
         work%logP = logP_tmp
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%input_data_use(:) = work%input_data_use_tmp(:)
         ! update the running means
         if( update_running_means( work, err ) == RETURN_FAIL ) goto 800
         if( work%beta_current_reject_run > work%stuck_limit ) then
            call err_handle(err, 1, &
                 comment = "Random-walk Metropolis got stuck" )
            goto 10
         end if
         !#######################
         aborted = .false.
      end do
      !### end main iteration
10    continue
      if( aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = work%iter )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_rwm_log_linear
    !##################################################################
    integer(kind=our_int) function compute_scale_rwm( &
        df, scale_fac, work, err ) result(answer)
      ! beta_scale for RWM proposal from work%vhat_beta_rwm
      ! also computes beta_scale_sqrt
      implicit none
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(in) :: scale_fac
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      real(kind=our_dble) :: rtmp
      character(len=*), parameter :: &
           subname = "compute_scale_rwm"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      rtmp = ( df / &
           ( df + real( work%p, our_dble) ) ) &
                 * ( 1.D0 / ( scale_fac**2 ) )
      rtmp = 1.D0 / rtmp
      work%beta_scale(:,:) = rtmp * work%vhat_beta_rwm
      work%beta_scale_sqrt(:,:) = work%beta_scale(:,:)
      if( cholesky_in_place( work%beta_scale_sqrt, err ) &
           == RETURN_FAIL ) goto 200
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Matrix vhat_beta_rm not positive definite" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_scale_rwm
    !##################################################################
    integer(kind=our_int) function compute_loglik_logprior( beta, &
         work, err, &
         use_cell_means, logprior, loglik ) &
         result(answer)
      ! Given the current estimate of the cell probs ( work%prob )
      ! or cell means (work%mu),
      ! accumulates the loglik and logprior
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare inputs
      logical, intent(in) :: use_cell_means
      ! declare outputs
      real(kind=our_dble), intent(out) :: loglik, logprior
      ! declare locals
      integer(kind=our_int) :: row, i, j
      real(kind=our_dble) :: flat_const, log_prob, sum
      real(kind=our_dble) :: ntotal, nprior, ndata, log_ntotal
      logical :: non_zero_cells_present
      character(len=*), parameter :: &
           subname = "compute_loglik_logprior"
      ! begin
      answer = RETURN_FAIL
      if( use_cell_means .and. (work%model_type /= "log-linear") ) goto 20
      ! initialize counts etc.
      logprior = 0.D0
      loglik = 0.D0
      nprior = 0.D0
      ndata = 0.D0
      if( work%model_type == "log-linear" ) then
         ! poisson correction
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            loglik = loglik - work%mu(i)
         end do
         ! ridge contribution
         sum = 0.D0
         do j = 1, work%p
            sum = sum + beta(j)**2
         end do
         logprior = logprior - work%ridge * sum / 2.D0
      end if
      ! apply flattening constant
      flat_const = work%flatten / &
           real( work%ncells_nonzero, our_dble )
      if( flat_const > 0.D0 ) then 
         do i = 1, work%ncells
            if( work%str_zero(i) ) cycle
            if( use_cell_means ) then
               if( work%mu(i) <= 0.D0 ) goto 200
               log_prob = log( work%mu(i) )
            else
               if( work%prob(i) <= 0.D0 ) goto 200
               log_prob = log( work%prob(i) )
            end if
            logprior = logprior + flat_const * &
                 log_prob
         end do
      end if
      nprior = nprior + work%flatten
      ! accumulate rest of logprior
      do row = 1, work%nrow_prior_data
         if( .not. work%prior_data_use(row) ) cycle
         work%begin_cycle = .true.
         work%cycle_done = .false.
         non_zero_cells_present = .false.
         sum = 0.D0
         do
            if( advance_to_next_cell( row, work%prior_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            i = work%cell
            if( .not. work%str_zero(i) ) then
               if( use_cell_means ) then
                  sum = sum + work%mu( i )
               else
                  sum = sum + work%prob( i )
               end if
               non_zero_cells_present = .true.
            end if
            if( work%cycle_done ) exit
         end do
         if( non_zero_cells_present .and. ( sum <= 0.D0 ) ) goto 200
         if( .not. non_zero_cells_present .and. &
               ( work%prior_data_freq(row) > 0.D0 ) ) goto 205
         if( sum < 0.D0 ) goto 210
         if( non_zero_cells_present ) &
              logprior = logprior + work%prior_data_freq(row) * log(sum)
         nprior = nprior + work%prior_data_freq(row)
      end do
      ! accumulate rest of loglik
      do row = 1, work%nrow_input_data
         if( .not. work%input_data_use(row) ) cycle
         work%begin_cycle = .true.
         work%cycle_done = .false.
         non_zero_cells_present = .false.
         sum = 0.D0
         do
            if( advance_to_next_cell( row, work%input_data, work, err ) &
                 == RETURN_FAIL ) goto 800
            i = work%cell
            if( .not. work%str_zero(i) ) then
               if( use_cell_means ) then
                  sum = sum + work%mu( i )
               else
                  sum = sum + work%prob( i )
               end if
               non_zero_cells_present = .true.
            end if
            if( work%cycle_done ) exit
         end do
         if( non_zero_cells_present .and. ( sum <= 0.D0 ) ) goto 200
         if( .not. non_zero_cells_present .and. &
               ( work%input_data_freq(row) > 0.D0 ) ) goto 206
         if( sum < 0.D0 ) goto 210
         if( non_zero_cells_present ) &
              loglik = loglik + work%input_data_freq(row) * log(sum)
         ndata = ndata + work%input_data_freq(row)
      end do
      ! poisson corrections for saturated model
      if( work%model_type == "saturated" ) then
         ntotal = nprior + ndata
         if( ntotal > 0.D0 ) then
            log_ntotal = log(ntotal)
            work%loglik = work%loglik + ndata * log_ntotal - ntotal
            work%logprior = work%logprior + nprior * log_ntotal
         else
            ! nothing to do, because nprior and ndata are zero
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Bad row in prior frame, positive freq for zero cells" )
      call err_handle(err, 3, iobs=row)
      goto 800
206   call err_handle(err, 1, &
            comment = "Bad row in data frame, positive freq for zero cells" )
      call err_handle(err, 3, iobs=row)
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_loglik_logprior
   !##################################################################
    integer(kind=our_int) function draw_candidate_beta_rwm( df, &
         work, err ) result(answer)
      ! Draws candidate beta from multivariate t proposal, storing
      ! the result in mmw%beta_can 
      ! Depends on beta_scale_sqrt, so should be run after 
      ! compute_scale_rwm
      implicit none
      real(kind=our_dble), intent(in) :: df
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum, rnorm, rchisq
      character(len=*), parameter :: &
           subname = "draw_candidate_beta_rwm"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      if( df <= 0.D0 ) goto 30
      !####
      ! draw independent t variates, put into wkpA
      if( rchisq_R( df, rchisq, err ) == RETURN_FAIL ) goto 800
      rtmp = sqrt( df / rchisq )
      do j = 1, work%p
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkpA(j) = rnorm * rtmp
      end do
      ! premultiply by beta_scale_sqrt, put into wkpB
      do j = 1, work%p
         sum = 0.D0
         do k = 1, j
            sum = sum + work%beta_scale_sqrt(j,k) * work%wkpA(k)
         end do
         work%wkpB(j) = sum
      end do
      ! add center to obtain candidate
      do j = 1, work%p
         work%beta_can(j) = work%wkpB(j) + work%beta(j)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_candidate_beta_rwm
    !##################################################################
    integer(kind=our_int) function run_rwm_step_beta( work, err ) &
         result( answer )
      ! performs one step of metropolis-hastings on beta, conditioning
      ! on the frequencies in work%freq
      ! uses the augmented-data poisson likelihood
      implicit none
      ! args
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: log_mh_ratio, mh_ratio, runif, &
         logprior_tmp, loglik_tmp
      logical :: accept
      character(len=*), parameter :: &
           subname = "run_rwm_step_beta"
      ! begin
      answer = RETURN_FAIL
      ! store mu and logmu to restore later if step is rejected
      work%mu_tmp(:) = work%mu(:)
      work%logmu_tmp(:) = work%logmu(:)
      !
      log_mh_ratio = - work%logP
      if( draw_candidate_beta_rwm( work%df_da, &
           work, err ) == RETURN_FAIL ) goto 800
      if( compute_mu_from_beta( work%beta_can, work, err ) &
           == RETURN_FAIL ) goto 800
      if( compute_loglik_logprior( work%beta_can, work, err, &
           use_cell_means = .true., logprior = logprior_tmp, &
           loglik = loglik_tmp ) == RETURN_FAIL ) goto 800
      log_mh_ratio = log_mh_ratio + logprior_tmp + loglik_tmp
      ! to prevent over/underflow
      if( log_mh_ratio > log_huge ) then
         mh_ratio = huge(0.D0)
      else if( log_mh_ratio < log_tiny ) then
         mh_ratio = 0.D0
      else
         mh_ratio = exp( log_mh_ratio )
      end if
      ! compare to uniform variate
      if( runif_R( runif, err ) == RETURN_FAIL ) goto 800
      if( runif <= mh_ratio ) then
         accept = .true.
         work%beta(:) = work%beta_can(:)
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = logprior_tmp + loglik_tmp
      else
         accept = .false.
         work%mu(:) = work%mu_tmp(:)
         work%logmu(:) = work%logmu_tmp(:)
      end if
      work%mh_ratios_beta( work%iter ) = mh_ratio
      work%mh_accept_beta( work%iter ) = accept
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function run_rwm_step_beta
    !##################################################################
    integer(kind=our_int) function run_approx_bayes_log_linear( &
         beta_series, prob_series, logp_series, imputed_freq_int, &
         packed_estimates_series, &
         n_iter_actual, n_sample_actual, n_imp_actual, &
         work, err ) result(answer)
      implicit none
      ! outputs 
      real(kind=our_dble), intent(out) :: beta_series(:,:)
      real(kind=our_dble), intent(out) :: prob_series(:,:)
      real(kind=our_dble), intent(out) :: logp_series(:)
      integer(kind=our_int), intent(out) :: imputed_freq_int(:,:)
      real(kind=our_dble), intent(out) :: packed_estimates_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter
      real(kind=our_dble) :: rtmp, logP_tmp, loglik_tmp, logprior_tmp
      logical :: aborted
      character(len=*), parameter :: subname = "run_approx_bayes_log_linear"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20 
      if( work%method /= "approxBayes" ) goto 30 
      ! initialize iteration counters
      work%iter = 0
      work%store_count = 0
      work%imp_count = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%prob_mean(:) = 0.D0
      work%freq_mean(:) = 0.D0
      if( work%n_estimates > 0 ) work%packed_estimates_mean = 0.D0
      if( size(beta_series) > 0 ) beta_series(:,:) = 0.D0
      if( size(prob_series) > 0 ) prob_series(:,:) = 0.D0
      if( size(logp_series) > 0 ) logp_series(:) = 0.D0
      if( size( imputed_freq_int) > 0 ) imputed_freq_int(:,:) = 0
      if( size(packed_estimates_series) > 0 ) &
           packed_estimates_series(:,:) = 0.D0
      work%beta_scale(:,:) = work%vhat_beta_rwm(:,:)
      work%beta_scale_sqrt(:,:) = work%beta_scale(:,:)
      if( cholesky_in_place( work%beta_scale_sqrt, err ) &
           == RETURN_FAIL ) goto 100
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_accept_rate = 0.D0
      work%store_this_iter = .true.
      if( work%impute_approx_bayes ) then
         work%imp_this_iter = .true.
      else
         work%imp_this_iter = .false.
      end if
      aborted = .false.
      ! main iteration
      do iter = 1, work%iter_approx_bayes
         aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = iter
         work%store_count = work%store_count + 1
         if( work%imp_this_iter ) work%imp_count = work%imp_count + 1
         !#######################
         if( work%iter == 1 ) then
            if( compute_mu_from_beta( work%beta, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( normalize_prob( work%mu, work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
            if( compute_estimates( work%prob, work, err ) &
                 == RETURN_FAIL ) goto 10
         else
            ! mu, logmu and prob are already consistent with beta
         end if
         if( compute_loglik_logprior( work%beta, work, err, &
              use_cell_means = .true., logprior = logprior_tmp, &
              loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = logprior_tmp + loglik_tmp
         if( work%iter == 1 ) work%start_logP = work%logP
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%logP
         if( draw_approx_bayes_beta( work, err ) == RETURN_FAIL ) goto 10
         if( compute_mu_from_beta( work%beta, work, err ) &
              == RETURN_FAIL ) goto 10
         if( normalize_prob( work%mu, work%prob, work, err ) &
              == RETURN_FAIL ) goto 10
         if( compute_estimates( work%prob, work, err ) &
              == RETURN_FAIL ) goto 10
         work%beta_accept_count = work%beta_accept_count + 1
         ! update counters and running estimates
         rtmp = 1.D0
         work%beta_accept_rate = work%beta_accept_rate + &
              ( rtmp - work%beta_accept_rate ) / &
              real( work%iter, our_int )
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 200
            if( work%store_count > size(beta_series, 1) ) goto 200
            if( work%store_count > size(logp_series) ) goto 200
            beta_series( work%store_count, : ) = work%beta(:)
            logp_series( work%store_count ) = work%logP
            if( work%save_prob_series ) then
               if( work%store_count > size(prob_series, 1) ) goto 200
               prob_series( work%store_count, : ) = work%prob(:)
            end if
            if( work%n_estimates > 0 ) then
               if( work%store_count > size(packed_estimates_series, 1) ) &
                    goto 200
               packed_estimates_series( work%store_count, : ) = &
                    work%packed_estimates(:)
            end if
         end if
         ! generate imputation, and put the result in imputed_freq_int
         work%freq_int_tmp(:) = work%freq_int(:)
         work%freq_tmp(:) = work%freq(:)
         logP_tmp = work%logP
         loglik_tmp = work%loglik
         logprior_tmp = work%logprior
         work%input_data_use_tmp(:) = work%input_data_use(:)
         work%input_data_use(:) = .true.
         if( run_istep( work, err, use_flatten = .false., &
              use_prior_data = .false., &
              use_input_data = .true. ) == RETURN_FAIL ) goto 10
         if( work%imp_this_iter ) then
            if( ( work%imp_count < 1 ) .or. &
                 ( work%imp_count > size(imputed_freq_int, 2) ) ) goto 200
            imputed_freq_int(:, work%imp_count) = work%freq_int(:)
         end if
         if( update_freq_mean( work, err ) == RETURN_FAIL ) goto 800
         work%freq_int(:) = work%freq_int_tmp(:)
         work%freq(:) = work%freq_tmp(:)
         work%logP = logP_tmp
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%input_data_use(:) = work%input_data_use_tmp(:)
         ! update the running means
         if( update_running_means( work, err ) == RETURN_FAIL ) goto 800
         !#######################
         aborted = .false.
      end do
      !### end main iteration
10    continue
      if( aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "Approximate Bayes procedure aborted" )
         call err_handle(err, 5, iiter = work%iter )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error trap
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for approximate Bayes" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Matrix vhat_beta_rwm not positive definite" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_approx_bayes_log_linear
    !##################################################################
    integer(kind=our_int) function draw_approx_bayes_beta( &
         work, err ) result(answer)
      ! Draws beta from multivariate normal
      ! puts result in mmw%beta
      ! Depends on beta_scale_sqrt
      implicit none
      ! declare workspaces
      type(workspace_type_cvam), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: sum, rnorm
      character(len=*), parameter :: &
           subname = "draw_approx_bayes_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%model_type /= "log-linear" ) goto 20
      !####
      do j = 1, work%p
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkpA(j) = rnorm
      end do
      ! premultiply by beta_scale_sqrt, put into wkpB
      do j = 1, work%p
         sum = 0.D0
         do k = 1, j
            sum = sum + work%beta_scale_sqrt(j,k) * work%wkpA(k)
         end do
         work%wkpB(j) = sum
      end do
      ! add center to obtain beta
      do j = 1, work%p
         work%beta(j) = work%wkpB(j) + work%beta_hat(j)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "There is no log-linear model" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_approx_bayes_beta
    !##################################################################
    integer(kind=our_int) function run_mlogit( n, p, r, x, y, &
         baseline, iter_max, criterion, &
         iter, converged_int, loglik, score, hess, &
         beta, beta_vec, vhat_beta_vec, pi_mat, err ) result(answer)
      ! fits a baseline-category logit model using Newton-Raphson
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: y(:,:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      integer(kind=our_int), intent(in) :: iter_max
      real(kind=our_dble), intent(in) :: criterion
      ! outputs
      integer(kind=our_int), intent(out) :: iter
      integer(kind=our_int), intent(out) :: converged_int
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
      real(kind=our_dble), intent(out) :: beta(:,:)
      real(kind=our_dble), intent(out) :: beta_vec(:)
      real(kind=our_dble), intent(out) :: vhat_beta_vec(:,:)
      real(kind=our_dble), intent(out) :: pi_mat(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      real(kind=our_dble), allocatable :: wkrA(:), wkrB(:), &
           wkprA(:), wkprprA(:,:), nvec(:), beta_vec_new(:)
      real(kind=our_dble) :: max_diff, sum
      integer(kind=our_int) :: status, i, j, k, jj
      logical :: converged, aborted
      character(len=*), parameter :: &
           subname = "run_mlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      if( check_mlogit_args( n, p, r, x, y, &
         baseline, iter_max, criterion, &
         score, hess, &
         beta, beta_vec, vhat_beta_vec, pi_mat, &
         err ) == RETURN_FAIL ) goto 800
      allocate( wkrA(r), wkrB(r), &
           wkprA( p*(r-1) ), wkprprA( p*(r-1), p*(r-1) ), &
           nvec(n), beta_vec_new( p*(r-1) ), &
           stat=status )
      if( status /= 0 ) goto 100
      if( create_nvec( y, nvec, err ) == RETURN_FAIL ) goto 800
      beta_vec_new(:) = 0.D0
      iter = 0
      converged = .false.
      converged_int = 0
      do
         if( converged ) exit
         if( iter >= iter_max ) exit 
         aborted = .true.  ! set to .false. at end of iteration
         iter = iter + 1
         beta_vec(:) = beta_vec_new(:)
         if( compute_score_hess_mlogit(x, y, baseline, nvec, &
              beta_vec, loglik, score, wkprprA, pi_mat, &
              err, wkrA, wkrB ) == RETURN_FAIL ) goto 20
         hess(:,:) = - wkprprA(:,:)
         if( cholesky_in_place(wkprprA, err ) == RETURN_FAIL ) then
            call err_handle(err, 1, &
                 comment = "Hessian not neg-definite" )
            goto 20
         end if
         if( invert_lower(wkprprA, err ) == RETURN_FAIL ) then
            call err_handle(err, 1, &
                 comment = "Hessian apparently singular" )
            goto 20
         end if
         if( premult_lower_by_transpose( wkprprA, &
              vhat_beta_vec, err) == RETURN_FAIL ) goto 800
         do j = 1, size(beta_vec)
            sum = 0.D0
            do k = 1, size(beta_vec)
               sum = sum + vhat_beta_vec(j,k) * score(k)
            end do
            beta_vec_new(j) = beta_vec(j) + sum
         end do
         !
         max_diff = 0.D0
         do i = 1, size(beta_vec)
            max_diff = max( max_diff, &
                 abs( beta_vec_new(i) - beta_vec(i) ) )
         end do
         if( max_diff <= criterion ) then
            converged = .true.
            converged_int = 1
         end if
         aborted = .false.
      end do
20    continue
      if( aborted ) then
         call err_handle(err, 1, &
              comment = "Newton-Raphson aborted" )
         call err_handle(err, 5, iiter = iter )
      end if
      if( .not. converged ) then
         call err_handle(err, 1, &
              comment = "Newton-Raphson failed to converge by" )
         call err_handle(err, 5, iiter = iter )
      end if
      if( aborted .or. ( .not. converged ) ) then
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      !
      beta_vec(:) = beta_vec_new(:)
      beta(:,baseline) = 0.D0 
      jj = 0
      do j = 1, r
         if(j == baseline) cycle
         do k = 1, p
            jj = jj + 1
            beta(k,j) = beta_vec_new(jj)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(wkrA) ) deallocate(wkrA)
      if( allocated(wkrB) ) deallocate(wkrB)
      if( allocated(wkprA) ) deallocate(wkprA)
      if( allocated(wkprprA) ) deallocate(wkprprA)
      if( allocated(nvec) ) deallocate(nvec)
      if( allocated(beta_vec_new) ) deallocate(beta_vec_new)
    end function run_mlogit
    !##################################################################
    integer(kind=our_int) function check_mlogit_args( n, p, r, x, y, &
         baseline, iter_max, criterion, score, hess, &
         beta, beta_vec, vhat_beta_vec, pi_mat, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: y(:,:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      integer(kind=our_int), intent(in) :: iter_max
      real(kind=our_dble), intent(in) :: criterion
      ! outputs
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
      real(kind=our_dble), intent(out) :: beta(:,:)
      real(kind=our_dble), intent(out) :: beta_vec(:)
      real(kind=our_dble), intent(out) :: vhat_beta_vec(:,:)
      real(kind=our_dble), intent(out) :: pi_mat(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      character(len=*), parameter :: &
           subname = "check_mlogit_args"
      ! begin
      answer = RETURN_FAIL
      !####
      if( ( size(x,1) /= n ) .or. ( size(x,2) /= p ) ) goto 20
      if( ( size(y,1) /= n ) .or. ( size(y,2) /= r ) ) goto 30
      if( ( baseline < 1 ) .or. ( baseline > r ) ) goto 40
      if( iter_max < 0 ) goto 50
      if( criterion <= 0.D0 ) goto 60
      if( size(score) /= p*(r-1) ) goto 70
      if( size(hess,1) /= p*(r-1) ) goto 80
      if( size(hess,2) /= p*(r-1) ) goto 80
      if( ( size(beta,1) /= p ) .or. ( size(beta,2) /= r ) )  goto 90
      if( size(beta_vec) /= p*(r-1) ) goto 100
      if( size(vhat_beta_vec,1) /= p*(r-1) ) goto 110
      if( size(vhat_beta_vec,2) /= p*(r-1) ) goto 110
      if( ( size(pi_mat,1) /= n ) .or. ( size(pi_mat,2) /= r ) ) goto 120
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "Argument x has incorrect size" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument y has incorrect size" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument baseline is out of bounds" )
      goto 800
50    call err_handle(err, 1, &
            comment = "Argument iter_max is negative" )
      goto 800
60    call err_handle(err, 1, &
            comment = "Argument criterion is not positive" )
      goto 800
70    call err_handle(err, 1, &
            comment = "Argument score has incorrect size" )
      goto 800
80    call err_handle(err, 1, &
            comment = "Argument hess has incorrect size" )
      goto 800
90    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument beta_vec has incorrect size" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Argument vhat_beta_vec has incorrect size" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Argument pi_mat has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function check_mlogit_args
    !##################################################################
    integer(kind=our_int) function create_nvec( y, nvec, err ) &
         result(answer)
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: y(:,:)
      ! outputs
      real(kind=our_dble), intent(out) :: nvec(:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: i, j
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "create_nvec"
      ! begin
      answer = RETURN_FAIL
      !####
      do i = 1, size(y,1)
         sum = 0.D0
         do j = 1, size(y,2)
            if( y(i,j) < 0.D0 ) goto 100
            sum = sum + y(i,j)
         end do
         nvec(i) = sum
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Negative value in response matrix y" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function create_nvec
    !##################################################################
    integer(kind=our_int) function compute_score_hess_mlogit( &
         x, y, baseline, nvec, beta_vec, loglik, &
         score, neg_hess, pi_mat, &
         err, wkrA, wkrB ) &
         result(answer)
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: y(:,:)
      integer(kind=our_int), intent(in) :: baseline
      real(kind=our_dble), intent(in) :: nvec(:)
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: neg_hess(:,:)
      real(kind=our_dble), intent(out) :: pi_mat(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      real(kind=our_dble), intent(inout) :: wkrA(:)
      real(kind=our_dble), intent(inout) :: wkrB(:)
      ! declare locals
      integer(kind=our_int) :: n, p, r, i, j, jj, k, kk, el, kprime, m
      real(kind=our_dble) :: q1, q2, q3
      character(len=*), parameter :: &
           subname = "compute_score_hess_mlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      n = size(x,1)
      p = size(x,2)
      r = size(y,2)
      loglik = 0.D0
      score(:) = 0.D0
      neg_hess(:,:) = 0.D0
      if( compute_pi_mat( x, baseline, beta_vec, pi_mat, err, wkrA, wkrB ) &
           == RETURN_FAIL ) goto 800
      do i = 1, n
         ! increment loglik
         do j = 1, r
            if( y(i,j) > 0.D0 ) then
               if( pi_mat(i,j) <= 0.D0 ) goto 200
               loglik = loglik + y(i,j) * log( pi_mat(i,j) )
            end if
         end do
         !### accumulate derivatives
         jj = 0
         do k = 1, r
            if(k == baseline) cycle
            !### first derivatives wrt beta, and minus second
            !### derivatives wrt beta, diagonal blocks
            q1 = y(i,k) - nvec(i) * pi_mat(i,k)
            q2 = nvec(i) * pi_mat(i,k) * ( 1.D0 - pi_mat(i,k) )
            do el = 1, p
               jj = jj + 1
               score(jj) = score(jj) + q1 * x(i,el)
               kk = jj - 1
               do m = el, p
                  kk = kk + 1
                  neg_hess(jj,kk)  =  neg_hess(jj,kk) + q2 * x(i,el) * x(i,m)
                  neg_hess(kk,jj)  =  neg_hess(jj,kk)
               end do
               !### minus 2nd derivatives, off-diagonal blocks
               do kprime = (k+1), r
                  if( kprime /=  baseline ) then
                     q3 = - nvec(i) * pi_mat(i,k) * pi_mat(i,kprime)
                     do m = 1, p
                        kk = kk + 1
                        neg_hess(jj,kk) = neg_hess(jj,kk) + &
                             q3 * x(i,el) * x(i,m)
                        neg_hess(kk,jj) = neg_hess(jj,kk)
                     end do
                  end if
               end do
            end do
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted logarithn of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_score_hess_mlogit
    !##################################################################
    integer(kind=our_int) function compute_pi_mat( x, baseline, &
         beta_vec, pi_mat, err, wkrA, wkrB ) result(answer)
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: x(:,:)
      integer(kind=our_int), intent(in) :: baseline
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: pi_mat(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      real(kind=our_dble), intent(inout) :: wkrA(:)
      real(kind=our_dble), intent(inout) :: wkrB(:)
      ! declare locals
      integer(kind=our_int) :: n, p, r, i, j, posn, k
      real(kind=our_dble) :: sum, eta_max
      character(len=*), parameter :: &
           subname = "compute_pi_mat"
      ! begin
      answer = RETURN_FAIL
      !####
      n = size(x,1)
      p = size(x,2)
      r = size(pi_mat,2)
      do i = 1, n
         posn = 0
         do j = 1, r
            if( j == baseline ) then
               wkrA(j) = 0.D0
               cycle
            end if
            sum = 0.D0
            do k = 1, p
               posn = posn + 1
               sum = sum + x(i,k) * beta_vec(posn)
            end do
            wkrA(j) = sum
         end do
         eta_max = log_tiny
         do j = 1, r
            eta_max = max(eta_max, wkrA(j))
         end do
         wkrA(:) = wkrA(:) - eta_max
         sum = 0.D0
         do j = 1, r
            if( wkrA(j) < log_tiny ) then
               wkrB(j) = 0.D0
            else if( wkrA(j) > log_huge ) then
               goto 100
            else
               wkrB(j) = exp( wkrA(j) )
            end if
            sum = sum + wkrB(j)
         end do
         if( sum == 0.D0 ) goto 200
         do j = 1, r
            pi_mat(i,j) = wkrB(j) / sum
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Overflow; fitted value became too large" )
      call err_handle(err, 3, iobs=i)
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      call err_handle(err, 3, iobs=i)
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_pi_mat
    !##################################################################
    integer(kind=our_int) function run_mlogit_loglik_derivs( n, p, r, &
         x, y, baseline, beta_vec, &
         loglik, score, hess, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: y(:,:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
!      real(kind=our_dble), intent(out) :: neg_hess_inv(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      real(kind=our_dble), allocatable :: wkrA(:), wkrB(:), &
           wkprprA(:,:), nvec(:), pi_mat(:,:)
      integer(kind=our_int) :: status
      character(len=*), parameter :: &
           subname = "run_mlogit_loglik_derivs"
      ! begin
      answer = RETURN_FAIL
      !####
      if( check_mlogit_loglik_derivs_args( n, p, r, x, y, &
         baseline, beta_vec, score, hess, &
         err ) == RETURN_FAIL ) goto 800
      allocate( wkrA(r), wkrB(r), &
           wkprprA( p*(r-1), p*(r-1) ), &
           nvec(n), pi_mat(n,r), &
           stat=status )
      if( status /= 0 ) goto 100
      if( create_nvec( y, nvec, err ) == RETURN_FAIL ) goto 800
      if( compute_score_hess_mlogit(x, y, baseline, nvec, &
           beta_vec, loglik, score, wkprprA, pi_mat, &
           err, wkrA, wkrB ) == RETURN_FAIL ) goto 800
      hess(:,:) = - wkprprA(:,:)
!      neg_hess_inv(:,:) = 0.D0
!      if( cholesky_in_place( wkprprA, err ) == RETURN_FAIL ) then
!         call err_handle(err, 1, &
!              comment = "Hessian matrix not neg-def" )
!         goto 20
!      end if
!      if( invert_lower( wkprprA, err ) == RETURN_FAIL ) then
!         call err_handle(err, 1, &
!              comment = "Hessian matrix apparently singular" )
!         goto 20
!      end if
!      if( premult_lower_by_transpose( wkprprA, neg_hess_inv, &
!           err) == RETURN_FAIL ) goto 800
!20    continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(wkrA) ) deallocate(wkrA)
      if( allocated(wkrB) ) deallocate(wkrB)
      if( allocated(wkprprA) ) deallocate(wkprprA)
      if( allocated(nvec) ) deallocate(nvec)
      if( allocated(pi_mat) ) deallocate(pi_mat)
    end function run_mlogit_loglik_derivs
    !##################################################################
    integer(kind=our_int) function check_mlogit_loglik_derivs_args( &
         n, p, r, x, y, baseline, beta_vec, score, hess, &
         err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: y(:,:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      character(len=*), parameter :: &
           subname = "check_mlogit_loglik_derivs_args"
      ! begin
      answer = RETURN_FAIL
      !####
      if( ( size(x,1) /= n ) .or. ( size(x,2) /= p ) ) goto 20
      if( ( size(y,1) /= n ) .or. ( size(y,2) /= r ) ) goto 30
      if( ( baseline < 1 ) .or. ( baseline > r ) ) goto 40
      if( size(score) /= p*(r-1) ) goto 70
      if( size(hess,1) /= p*(r-1) ) goto 80
      if( size(hess,2) /= p*(r-1) ) goto 80
      if( size(beta_vec) /= p*(r-1) ) goto 100
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "Argument x has incorrect size" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument y has incorrect size" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument baseline is out of bounds" )
      goto 800
70    call err_handle(err, 1, &
            comment = "Argument score has incorrect size" )
      goto 800
80    call err_handle(err, 1, &
            comment = "Argument hess has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument beta_vec has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function check_mlogit_loglik_derivs_args
    !##################################################################
    integer(kind=our_int) function run_lcprev_loglik_derivs( n, p, r, &
         x, lik_mat, freq, baseline, beta_vec, &
         loglik, score, hess, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: lik_mat(:,:)
      real(kind=our_dble), intent(in) :: freq(:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
!      real(kind=our_dble), intent(out) :: neg_hess_inv(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      real(kind=our_dble), allocatable :: wkrA(:), wkrB(:), &
           wkprprA(:,:), pi_mat(:,:), pi_star(:,:)
      integer(kind=our_int) :: status
      character(len=*), parameter :: &
           subname = "run_lcprev_loglik_derivs"
      ! begin
      answer = RETURN_FAIL
      !####
      if( check_lcprev_loglik_derivs_args( n, p, r, x, lik_mat, freq, &
         baseline, beta_vec, score, hess, &
         err ) == RETURN_FAIL ) goto 800
      allocate( wkrA(r), wkrB(r), &
           wkprprA( p*(r-1), p*(r-1) ), &
           pi_mat(n,r), pi_star(n,r), &
           stat=status )
      if( status /= 0 ) goto 100
      if( compute_score_hess_lcprev(x, lik_mat, baseline, freq, &
           beta_vec, loglik, score, wkprprA, pi_mat, pi_star, &
           err, wkrA, wkrB ) == RETURN_FAIL ) goto 800
      hess(:,:) = - wkprprA(:,:)
!      neg_hess_inv(:,:) = 0.D0
!      if( cholesky_in_place( wkprprA, err ) == RETURN_FAIL ) then
!         call err_handle(err, 1, &
!              comment = "Hessian matrix not neg-def" )
!         goto 20
!      end if
!      if( invert_lower( wkprprA, err ) == RETURN_FAIL ) then
!         call err_handle(err, 1, &
!              comment = "Hessian matrix apparently singular" )
!         goto 20
!      end if
!      if( premult_lower_by_transpose( wkprprA, neg_hess_inv, &
!           err) == RETURN_FAIL ) goto 800
!20    continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated(wkrA) ) deallocate(wkrA)
      if( allocated(wkrB) ) deallocate(wkrB)
      if( allocated(wkprprA) ) deallocate(wkprprA)
      if( allocated(pi_mat) ) deallocate(pi_mat)
      if( allocated(pi_star) ) deallocate(pi_star)
    end function run_lcprev_loglik_derivs
    !##################################################################
    integer(kind=our_int) function check_lcprev_loglik_derivs_args( &
         n, p, r, x, lik_mat, freq, baseline, &
         beta_vec, score, hess, &
         err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n   ! number of observations
      integer(kind=our_int), intent(in) :: p   ! columns of x
      integer(kind=our_int), intent(in) :: r   ! response categories
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: lik_mat(:,:)
      real(kind=our_dble), intent(in) :: freq(:)
      integer(kind=our_int), intent(in) :: baseline   ! baseline category
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      ! declare locals
      character(len=*), parameter :: &
           subname = "check_lcprev_loglik_derivs_args"
      ! begin
      answer = RETURN_FAIL
      !####
      if( ( size(x,1) /= n ) .or. ( size(x,2) /= p ) ) goto 20
      if( ( size(lik_mat,1) /= n ) .or. ( size(lik_mat,2) /= r ) ) goto 30
      if( size(freq) /= n ) goto 35
      if( ( baseline < 1 ) .or. ( baseline > r ) ) goto 40
      if( size(score) /= p*(r-1) ) goto 70
      if( size(hess,1) /= p*(r-1) ) goto 80
      if( size(hess,2) /= p*(r-1) ) goto 80
      if( size(beta_vec) /= p*(r-1) ) goto 100
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "Argument x has incorrect size" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument lik_mat has incorrect size" )
      goto 800
35    call err_handle(err, 1, &
            comment = "Argument freq has incorrect size" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument baseline is out of bounds" )
      goto 800
70    call err_handle(err, 1, &
            comment = "Argument score has incorrect size" )
      goto 800
80    call err_handle(err, 1, &
            comment = "Argument hess has incorrect size" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Argument beta_vec has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function check_lcprev_loglik_derivs_args
    !##################################################################
    integer(kind=our_int) function compute_score_hess_lcprev( &
         x, lik_mat, baseline, freq, beta_vec, loglik, &
         score, neg_hess, pi_mat, pi_star, &
         err, wkrA, wkrB ) &
         result(answer)
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(in) :: lik_mat(:,:)
      integer(kind=our_int), intent(in) :: baseline
      real(kind=our_dble), intent(in) :: freq(:)
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: neg_hess(:,:)
      real(kind=our_dble), intent(out) :: pi_mat(:,:)
      real(kind=our_dble), intent(out) :: pi_star(:,:)
      ! declare workspaces
      type(error_type), intent(inout) :: err
      real(kind=our_dble), intent(inout) :: wkrA(:)
      real(kind=our_dble), intent(inout) :: wkrB(:)
      ! declare locals
      integer(kind=our_int) :: n, p, r, i, j, k, l, m, posnA, posnB, c
      real(kind=our_dble) :: sum, q1, q2, q3, ifunA, ifunB
      character(len=*), parameter :: &
           subname = "compute_score_hess_lcprev"
      ! begin
      answer = RETURN_FAIL
      !####
      n = size(x,1)
      p = size(x,2)
      r = size(lik_mat,2)
      loglik = 0.D0
      score(:) = 0.D0
      neg_hess(:,:) = 0.D0
      if( compute_pi_mat( x, baseline, beta_vec, pi_mat, err, wkrA, wkrB ) &
           == RETURN_FAIL ) goto 800
      do i = 1, n
         sum = 0.D0
         do j = 1, r
            if( pi_mat(i,j) < 0.D0 ) goto 100
            if( lik_mat(i,j) < 0.D0 ) goto 150
            wkrA(j) = pi_mat(i,j) * lik_mat(i,j)
            sum = sum + wkrA(j)
         end do
         if( sum == 0.D0 ) goto 200
         if( sum <= 0.D0 ) goto 250
         do j = 1, r
            pi_star(i,j) = wkrA(j) / sum
         end do
         loglik = loglik + freq(i) * log(sum)
         posnA = 0
         do j = 1, r
            if( j == baseline ) cycle
            do k = 1, p
               posnA = posnA + 1
               score(posnA) = score(posnA) + freq(i) * &
                    ( pi_star(i,j) - pi_mat(i,j) ) * x(i,k)
            end do
         end do
         posnA = 0
         do j = 1, r
            if( j == baseline ) cycle
            do k = 1, p
               posnA = posnA + 1
               posnB = 0
               do l = 1, r
                  if( l == baseline ) cycle
                  do m = 1, p
                     posnB = posnB + 1
                     if( posnB <= posnA ) then
                        ifunA = 0.D0
                        if( l == j ) ifunA = 1.D0
                        q1 = pi_mat(i,j) * ( ifunA - pi_mat(i,l) )
                        q2 = 0.D0
                        do c = 1, r
                           ifunA = 0.D0
                           if( j == c ) ifunA = 1.D0
                           ifunB = 0.D0
                           if( l == c ) ifunB = 1.D0
                           q2 = q2 + pi_star(i,c) * &
                                ( ifunA - pi_mat(i,j) ) * &
                                ( ifunB - pi_mat(i,l) )
                        end do
                        q3 = ( pi_star(i,j) - pi_mat(i,j) ) * &
                             ( pi_star(i,l) - pi_mat(i,l) )
                        neg_hess( posnA, posnB ) = neg_hess( posnA, posnB ) &
                             + freq(i) * ( q1 - q2 + q3 ) * x(i,k) * x(i,m)
                     end if
                  end do
               end do
            end do
         end do
      end do
      do posnA = 1, p*(r-1)
         do posnB = 1, posnA
            neg_hess(posnB, posnA) = neg_hess(posnA, posnB)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps

100   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      call err_handle(err, 3, iobs=i)
      goto 800
150   call err_handle(err, 1, &
            comment = "Negative likelihood encountered" )
      call err_handle(err, 3, iobs=i)
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      call err_handle(err, 3, iobs=i)
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted logarithn of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_score_hess_lcprev
    !##################################################################
 end module cvam_engine
!#####################################################################
