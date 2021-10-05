!#####################################################################
!# wrapper functions for cvam shared object
!#####################################################################
subroutine fit_cvam_model( &
     model_type_int, method_int, &
     dim_vec, input_data, input_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     flatten, ridge, prior_data, prior_data_freq_int, &
     ctrl_real, ctrl_int, omit_data_int, &
     dim_vec_est, estimate_info, estimate_var_info, &
     ctrl_mcmc_int, ctrl_mcmc_real, dim_vec_mcmc, &
     prob, beta, beta_hat, vhat_beta_rwm, &
     iter, converged, max_diff, loglik_logP, lambda, &
     freq, freq_mean, freq_int, &
     score, vhat_beta, prob_mean, beta_mean, beta_cov_mat, &
     total_freq_use_prior, total_freq_use_data_int, &
     degrees_of_freedom, &
     packed_estimates, packed_estimates_mean, packed_SEs, &
     beta_series, prob_series, logp_series, imputed_freq_int, &
     packed_estimates_series, &
     n_actual, mh_accept_rate, &
     start_logP, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: input_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: input_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: flatten
   real(kind=our_dble), intent(in) :: ridge
   integer(kind=our_int), intent(in) :: prior_data( dim_vec(7), dim_vec(2) )
   integer(kind=our_int), intent(in) :: prior_data_freq_int( dim_vec(7) )
   real(kind=our_dble), intent(in) :: ctrl_real(4)
   integer(kind=our_int), intent(in) :: ctrl_int(5)
   integer(kind=our_int), intent(in) :: omit_data_int
   integer(kind=our_int), intent(in) :: dim_vec_est(4)
   integer(kind=our_int), intent(in) :: estimate_info( dim_vec_est(1), 4) 
   integer(kind=our_int), intent(in) :: estimate_var_info( dim_vec_est(2), 3) 
   integer(kind=our_int), intent(in) :: ctrl_mcmc_int(9)
   real(kind=our_dble), intent(in) :: ctrl_mcmc_real(5)
   integer(kind=our_int), intent(in) :: dim_vec_mcmc(3)
   ! inout arguments
   real(kind=our_dble), intent(inout) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(inout) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(inout) :: beta_hat( dim_vec(6) )
   real(kind=our_dble), intent(inout) :: vhat_beta_rwm( dim_vec(6), dim_vec(6) )
   ! output arguments
   integer(kind=our_int), intent(out) :: iter
   integer(kind=our_int), intent(out) :: converged
   real(kind=our_dble), intent(out) :: max_diff
   real(kind=our_dble), intent(out) :: loglik_logP( dim_vec(8), 2 )
   real(kind=our_dble), intent(out) :: lambda( dim_vec(6) )
   real(kind=our_dble), intent(out) :: freq( dim_vec(4) )
   real(kind=our_dble), intent(out) :: freq_mean( dim_vec(4) )
   integer(kind=our_int), intent(out) :: freq_int( dim_vec(4) )
   real(kind=our_dble), intent(out) :: score( dim_vec(6) )
   real(kind=our_dble), intent(out) :: vhat_beta( dim_vec(6), dim_vec(6) )
   real(kind=our_dble), intent(out) :: prob_mean( dim_vec(4) )
   real(kind=our_dble), intent(out) :: beta_mean( dim_vec(6) )
   real(kind=our_dble), intent(out) :: beta_cov_mat( dim_vec(6), dim_vec(6) )
   real(kind=our_dble), intent(out) :: total_freq_use_prior
   integer(kind=our_int), intent(out) :: total_freq_use_data_int
   integer(kind=our_int), intent(out) :: degrees_of_freedom
   real(kind=our_dble), intent(out) :: packed_estimates( dim_vec_est(3) )
   real(kind=our_dble), intent(out) :: packed_estimates_mean( dim_vec_est(3) )
   real(kind=our_dble), intent(out) :: packed_SEs( dim_vec_est(4) )
   real(kind=our_dble), intent(out) :: &
        beta_series( dim_vec_mcmc(1), dim_vec(6) )
   real(kind=our_dble), intent(out) :: &
        prob_series( dim_vec_mcmc(2), dim_vec(4) )
   real(kind=our_dble), intent(out) :: logp_series( dim_vec_mcmc(1) )
   integer(kind=our_int), intent(out) :: &
        imputed_freq_int( dim_vec(4), dim_vec_mcmc(3) )
   real(kind=our_dble), intent(out) :: &
        packed_estimates_series( dim_vec_mcmc(1), dim_vec_est(3) )
   integer(kind=our_int), intent(out) :: n_actual(3)
   real(kind=our_dble), intent(out) :: mh_accept_rate
   real(kind=our_dble), intent(out) :: start_logP
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   logical :: converged_logical
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_model( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        flatten, ridge, prior_data, prior_data_freq_int, &
        ctrl_real(1), &  ! start_val_jitter
        ctrl_real(2), &  ! crit_EM
        ctrl_real(3), &  ! crit_NR
        ctrl_real(4), &  ! crit_boundary
        ctrl_int(1),  &  ! iter_max_EM
        ctrl_int(2),  &  ! iter_max_NR
        ctrl_int(3),  &  ! start_val_use_int
        ctrl_int(4),  &  ! start_val_default_int
        ctrl_int(5),  &  ! exclude_all_na_int
        omit_data_int, &
        estimate_info, estimate_var_info, &
        ctrl_mcmc_int(1), &  ! iter_mcmc
        ctrl_mcmc_int(2), &  ! burn_mcmc
        ctrl_mcmc_int(3), &  ! thin_mcmc
        ctrl_mcmc_int(4), &  ! impute_every
        ctrl_mcmc_int(5), &  ! save_prob_series_int
        ctrl_mcmc_int(6), &  ! type_mcmc_int
        ctrl_mcmc_int(7), &  ! stuck_limit
        ctrl_mcmc_int(8), &  ! iter_approx_bayes
        ctrl_mcmc_int(9), &  ! impute_approx_bayes_int
        ctrl_mcmc_real(1), &  ! df_da
        ctrl_mcmc_real(2), &  ! step_size_da
        ctrl_mcmc_real(3), &  ! scale_fac_da
        ctrl_mcmc_real(4), &  ! df_rwm
        ctrl_mcmc_real(5), &  ! scale_fac_rwm
        prob, beta, beta_hat, vhat_beta_rwm, &
        work, err, &
        iter, converged_logical, max_diff, &
        loglik_logP(:,1), loglik_logP(:,2), &
        lambda, freq, freq_mean, freq_int, &
        score, vhat_beta, prob_mean, beta_mean, beta_cov_mat, &
        total_freq_use_prior, total_freq_use_data_int, &
        degrees_of_freedom, &
        packed_estimates, packed_estimates_mean, packed_SEs, &
        beta_series, prob_series, logp_series, imputed_freq_int, &
        packed_estimates_series, &
        n_actual(1), & ! n_iter_actual
        n_actual(2), & ! n_sample_actual
        n_actual(3), & ! n_imp_actual
        mh_accept_rate, &
        start_logP &
        ) == RETURN_FAIL ) goto 800
   !
   if( converged_logical ) then
      converged = 1
   else
      converged = 0
   end if
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine fit_cvam_model
!#####################################################################
subroutine cvam_estimate_em( &
     model_type_int, method_int, &
     dim_vec, input_data, input_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     prob, beta, vhat_beta, &
     dim_vec_est, estimate_info, estimate_var_info, skip_SEs_int, &
     packed_estimates, packed_SEs, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: input_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: input_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(in) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(out) :: vhat_beta( dim_vec(6), dim_vec(6) )
   integer(kind=our_int), intent(in) :: dim_vec_est(4)
   integer(kind=our_int), intent(in) :: estimate_info( dim_vec_est(1), 4) 
   integer(kind=our_int), intent(in) :: estimate_var_info( dim_vec_est(2), 3) 
   integer(kind=our_int), intent(in) :: skip_SEs_int
   ! outputs
   real(kind=our_dble), intent(out) :: packed_estimates( dim_vec_est(3) )
   real(kind=our_dble), intent(out) :: packed_SEs( dim_vec_est(4) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_estimate_em( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        estimate_info, estimate_var_info, skip_SEs_int, &
        work, err, &
        packed_estimates, packed_SEs &
        ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine cvam_estimate_em
!#####################################################################
subroutine cvam_predict_em( &
     model_type_int, method_int, &
     dim_vec, pred_data, pred_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     prob, beta, vhat_beta, &
     dim_vec_pred, predict_var_info, &
     pred_mat, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: pred_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: pred_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(in) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(in) :: vhat_beta( dim_vec(6), dim_vec(6) )
   integer(kind=our_int), intent(in) :: dim_vec_pred(3)
   integer(kind=our_int), intent(in) :: predict_var_info( dim_vec_pred(3), 2) 
   ! outputs
   real(kind=our_dble), intent(out) :: pred_mat( dim_vec_pred(1), &
        dim_vec_pred(2) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_predict_em( &
        model_type_int, method_int, &
        pred_data, pred_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        predict_var_info, &
        work, err, &
        pred_mat &
        ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine cvam_predict_em
!#####################################################################
subroutine cvam_impute_freq( &
     model_type_int, method_int, &
     dim_vec, input_data, input_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     prob, beta, vhat_beta, synthetic_int, &
     result_mat, result_freq_int, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: input_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: input_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(in) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(in) :: vhat_beta( dim_vec(6), dim_vec(6) )
   integer(kind=our_int), intent(in) :: synthetic_int
   ! outputs
   integer(kind=our_int), intent(out) :: result_mat( dim_vec(4), dim_vec(2) )
   integer(kind=our_int), intent(out) :: result_freq_int( dim_vec(4) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_impute_freq( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, synthetic_int, &
        work, err, &
        result_mat, result_freq_int &
        ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine cvam_impute_freq
!#####################################################################
subroutine cvam_impute_microdata( &
     model_type_int, method_int, &
     dim_vec, input_data, input_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     prob, beta, vhat_beta, synthetic_int, &
     result_mat, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: input_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: input_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(in) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(in) :: vhat_beta( dim_vec(6), dim_vec(6) )
   integer(kind=our_int), intent(in) :: synthetic_int
   ! outputs
   integer(kind=our_int), intent(out) :: result_mat( dim_vec(1), dim_vec(2) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_impute_microdata( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, synthetic_int, &
        work, err, &
        result_mat &
        ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine cvam_impute_microdata
!#####################################################################
subroutine cvam_lik( &
     model_type_int, method_int, &
     dim_vec, input_data, input_data_freq_int, &
     n_levels_matrix, packed_map, &
     model_matrix, offset, str_zero_int, &
     prob, beta, vhat_beta, &
     dim_vec_lik, lik_var_info, &
     lik_values, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_cvam_model
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: model_type_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: dim_vec(8)
   integer(kind=our_int), intent(in) :: input_data( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: input_data_freq_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: n_levels_matrix( dim_vec(2), 4 )
   integer(kind=our_int), intent(in) :: packed_map( dim_vec(3) )
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(5), dim_vec(6) )
   real(kind=our_dble), intent(in) :: offset( dim_vec(5) )
   integer(kind=our_int), intent(in) :: str_zero_int( dim_vec(4) )
   real(kind=our_dble), intent(in) :: prob( dim_vec(4) )
   real(kind=our_dble), intent(in) :: beta( dim_vec(6) )
   real(kind=our_dble), intent(in) :: vhat_beta( dim_vec(6), dim_vec(6) )
   integer(kind=our_int), intent(in) :: dim_vec_lik(2)
   integer(kind=our_int), intent(in) :: lik_var_info( dim_vec_lik(2), 3) 
   ! outputs
   real(kind=our_dble), intent(out) :: lik_values( dim_vec(1) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_cvam) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_cvam_lik( &
        model_type_int, method_int, &
        input_data, input_data_freq_int, &
        n_levels_matrix, packed_map, &
        model_matrix, offset, str_zero_int, &
        prob, beta, vhat_beta, &
        lik_var_info, &
        work, err, &
        lik_values &
        ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_cvam( work, err )
   ijunk = put_randgen_state_R(err)
 end subroutine cvam_lik
!#####################################################################
subroutine cvam_mlogit( n, p, r, x, y, baseline, iter_max, criterion, &
     iter, converged_int, loglik, &
     score, hess, beta, beta_vec, vhat_beta_vec, pi_mat, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_mlogit
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n
   integer(kind=our_int), intent(in) :: p
   integer(kind=our_int), intent(in) :: r
   real(kind=our_dble), intent(in) :: x(n,p)
   real(kind=our_dble), intent(in) :: y(n,r)
   integer(kind=our_int), intent(in) :: baseline
   integer(kind=our_int), intent(in) :: iter_max
   real(kind=our_dble), intent(in) :: criterion
   ! outputs
   integer(kind=our_int), intent(out) :: iter 
   integer(kind=our_int), intent(out) :: converged_int
   real(kind=our_dble), intent(out) :: loglik
   real(kind=our_dble), intent(out) :: score( p*(r-1) )
   real(kind=our_dble), intent(out) :: hess( p*(r-1), p*(r-1) )
   real(kind=our_dble), intent(out) :: beta(p,r)
   real(kind=our_dble), intent(out) :: beta_vec( p*(r-1) )
   real(kind=our_dble), intent(out) :: vhat_beta_vec( p*(r-1), p*(r-1) )
   real(kind=our_dble), intent(out) :: pi_mat(n,r)
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( run_mlogit( n, p, r, x, y, baseline, iter_max, criterion, &
        iter, converged_int, loglik, score, hess, &
        beta, beta_vec, vhat_beta_vec, pi_mat, &
        err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
 end subroutine cvam_mlogit
!#####################################################################
subroutine cvam_mlogit_loglik_derivs( n, p, r, x, y, baseline, &
     beta_vec, loglik, score, hess, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_mlogit_loglik_derivs
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n
   integer(kind=our_int), intent(in) :: p
   integer(kind=our_int), intent(in) :: r
   real(kind=our_dble), intent(in) :: x(n,p)
   real(kind=our_dble), intent(in) :: y(n,r)
   integer(kind=our_int), intent(in) :: baseline
   real(kind=our_dble), intent(in) :: beta_vec( p*(r-1) )
   ! outputs
   real(kind=our_dble), intent(out) :: loglik
   real(kind=our_dble), intent(out) :: score( p*(r-1) )
   real(kind=our_dble), intent(out) :: hess( p*(r-1), p*(r-1) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( run_mlogit_loglik_derivs( n, p, r, x, y, baseline, beta_vec, &
        loglik, score, hess, err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
 end subroutine cvam_mlogit_loglik_derivs
!#####################################################################
subroutine cvam_lcprev_loglik_derivs( n, p, r, x, lik_mat, &
     freq, baseline, beta_vec, &
     loglik, score, hess, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_mlogit_loglik_derivs
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use cvam_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n
   integer(kind=our_int), intent(in) :: p
   integer(kind=our_int), intent(in) :: r
   real(kind=our_dble), intent(in) :: x(n,p)
   real(kind=our_dble), intent(in) :: lik_mat(n,r)
   real(kind=our_dble), intent(in) :: freq(n)
   integer(kind=our_int), intent(in) :: baseline
   real(kind=our_dble), intent(in) :: beta_vec( p*(r-1) )
   ! outputs
   real(kind=our_dble), intent(out) :: loglik
   real(kind=our_dble), intent(out) :: score( p*(r-1) )
   real(kind=our_dble), intent(out) :: hess( p*(r-1), p*(r-1) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( run_lcprev_loglik_derivs( n, p, r, x, lik_mat, freq, &
        baseline, beta_vec, &
        loglik, score, hess, err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
 end subroutine cvam_lcprev_loglik_derivs
!#####################################################################
