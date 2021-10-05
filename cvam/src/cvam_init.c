#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cvam_estimate_em)(int *model_type_int, int *method_int, int *dim_vec, int *input_data, int *input_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *prob, double *beta, double *vhat_beta, int *dim_vec_est, int *estimate_info, int *estimate_var_info, int *skip_SEs_int, double *packed_estimates, double *packed_SEs, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_impute_freq)(int *model_type_int, int *method_int, int *dim_vec, int *imput_data, int *input_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *prob, double *beta, double *vhat_beta, int *synthetic_int, int *result_mat, int *result_freq_int, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_impute_microdata)(int *model_type_int, int *method_int, int *dim_vec, int *input_data, int *input_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *prob, double *beta, double *vhat_beta, int *synthetic_int, int *result_mat, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_lcprev_loglik_derivs)(int *n, int *p, int *r, double *x, double *lik_mat, double *freq, int *baseline, double *beta_vec, double *loglik, double *score, double *hess, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_lik)(int *model_type_int, int *method_int, int *dim_vec, int *input_data, int *input_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *prob, double *beta, double *vhat_beta, int *dim_vec_lik, int *lik_var_info, double *lik_values, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
			      extern void F77_NAME(cvam_mlogit)(int *n, int *p, int *r, double *x, double *y, int *baseline, int *iter_max, double *criterion, int *iter, int *converged_int, double *loglik, double *score, double *hess, double *beta, double *beta_vec, double *vhat_beta_vec, double *pi_mat, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_mlogit_loglik_derivs)(int *n, int *p, int *r, double *x, double *y, int *baseline, double *beta_vec, double *loglik, double *score, double *hess, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(cvam_predict_em)(int *model_type_int, int *method_int, int *dim_vec, int *pred_data, int *pred_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *prob, double *beta, double *vhat_beta, int *dim_vec_pred, int *predict_var_info, double *pred_mat, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);
extern void F77_NAME(fit_cvam_model)(int *model_type_int, int *method_int, int *dim_vec, int *input_data, int *input_data_freq_int, int *n_levels_matrix, int *packed_map, double *model_matrix, double *offset, int *str_zero_int, double *flatten, double *ridge, int *prior_data, int *prior_data_freq_int, double *ctrl_real, int *ctrl_int, int *omit_data_int, int *dim_vec_est, int *estimate_info, int *estimate_var_info, int *ctrl_mcmc_int, double *ctrl_mcmc_real, int *dim_vec_mcmc, double *prob, double *beta, double *beta_hat, double *vhat_beta_rwm, int *iter, int *converged, double *max_diff, double *loglik_logP, double *lambda, double *freq, double *freq_mean, int *freq_int, double *score, double *vhat_beta, double *prob_mean, double *beta_mean, double *beta_cov_mat, double *total_freq_use_prior, int *total_freq_use_data, int *degrees_of_freedom, double *packed_estimates, double *packed_estimates_mean, double *packed_SEs, double *beta_series, double *prob_series, double *logp_series, int *imputed_freq_int, double *packed_estimates_series, int *n_actual, double *mh_accept_rate, double *start_logP, int *status, int *msg_len_max, int *msg_codes, int *msg_len_actual);

static const R_FortranMethodDef FortranEntries[] = {
    {"cvam_estimate_em",          (DL_FUNC) &F77_NAME(cvam_estimate_em),          23},
    {"cvam_impute_freq",          (DL_FUNC) &F77_NAME(cvam_impute_freq),          20},
    {"cvam_impute_microdata",     (DL_FUNC) &F77_NAME(cvam_impute_microdata),     19},
    {"cvam_lcprev_loglik_derivs", (DL_FUNC) &F77_NAME(cvam_lcprev_loglik_derivs), 15},
    {"cvam_lik",                  (DL_FUNC) &F77_NAME(cvam_lik),                  20},
    {"cvam_mlogit",               (DL_FUNC) &F77_NAME(cvam_mlogit),               21},
    {"cvam_mlogit_loglik_derivs", (DL_FUNC) &F77_NAME(cvam_mlogit_loglik_derivs), 14},
    {"cvam_predict_em",           (DL_FUNC) &F77_NAME(cvam_predict_em),           20},
    {"fit_cvam_model",            (DL_FUNC) &F77_NAME(fit_cvam_model),            58},
    {NULL, NULL, 0}
};

void R_init_cvam(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
