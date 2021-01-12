#' M-estimator for threshold spatial dynamic panel data model
#' 
#' @import rJava
#' @import rCMA
#' @importFrom matrixcalc is.singular.matrix
#' 
#' @param y matrix, containing regional index (first column), time index (second column, numeric) and dependent variable (third column, numeric).
#' @param x matrix, containing regional index (first column), time index (second column, numeric) and regressors (numeric).
#' @param w1 matrix, the spatial weight matrix. If w2 and w3 are supplied, the spatial weight matrix for spatial lag.
#' @param th logical, 
#' @param correction logical, whether to use adjusted score function. Default value is True.
#' @param max_try integer, maximum attempt for the solver. Default value is 5.
#' @param all_er logical, whether to output Hessian and Gamma matrix based se. Ignored if correction is set to False. Default value is False.
#' @param true_range logical, whether to used the accurate stationary check. Default value is False due to performance reasons.
#' @param residual logical, whether to output the residual. Default value is False.
#' @param w2 matrix, the spatial weight matrix for spatio-temporal lag. Default value is the same as w1.
#' @param w3 matrix, the spatial weight matrix for spatial error. Default value is the same as w1.
#' @param no_tf logical, whether to account for time effect. Default value is True.
#' @param model character, indicates the model used for estimation, can be "full", "slm", "sem", "sltl". See Details.
#' @param rcpp logical, whether to use the rcpp implementation to calculate the score function. Default value is True.
#' @param cma_pop_multi integer, multiplier for the population size used in cmaes. Default value is 1.
#' 
#' @details Estimating threshold spatial dynamic panel data model with extended Yang(2018)'s M-estimator
#' \deqn{y_{ti} = \mu_{i} +\alpha_t+ x_{ti}\beta_{q} +\rho_{q} y_{t-1,i} + \lambda_{1q}\sum_{j=1}^{n}w_{1,ij}y_{tj} \\ 
#' \qquad + \lambda_{2q}\sum_{j=1}^{n}w_{2,ij}y_{t-1,i}+ u_{ti},\\
#' u_{ti} = \lambda_{3q}\sum_{j=1}^{n}w_{3,ij}u_{tj}+ v_{ti},i=1,\ldots,n,t=1,\ldots,T, q = 1,2}
#' Sub-models can be estimated by changing argument "model"
#' \itemize{
#' \item{"full"} {Full model}
#' \item{"slm"} {\eqn{\lambda_{2q} = \lambda_{3q} = 0}}
#' \item{"sem"} {\eqn{\lambda_{1q} = \lambda_{2q} = 0}}
#' \item{"sltl"} {\eqn{\lambda_{3q} = 0}}
#' }
#' 
#' 
#' @return A list of estimation results
#' \itemize{
#' \item{"coefficient"} {list, coefficients and standard errors}
#' \item{"vc_mat"} {matrix , variance-covariance matrix}
#' \item{"hes_mat"} {matrix, Hessian matrix}
#' \item{"gamma_mat"} {matrix, Gamma matrix}
#' \item{"model"} {character, model used for estimation}
#' \item{"residual"} {numeric, residual}
#' }
#' 
#' @references Yang, Z. (2018). Unified M-estimation of fixed-effects spatial dynamic models with short panels. Journal of Econometrics, 205(2), 423-447.
#' 
#' @examples
#' data(data_th, data_w)
#' msdpdth(y = data_th$y, x = data_th$x, w1 = data_w, th = data_th$th)
#' 
#' @export
#' 
#' 


msdpdth = function(y, x, w1, th, correction = T, max_try = 5, all_er = F, true_range = F, residual = F, w3 = w1, w2 = w1, no_tf = F, model = "full", rcpp = T, cma_pop_multi = 1){
  if (!correction) all_er = T
  #fit non-threshold model for initial value
  for (i in 1:max_try){
    nth_res = msdpd(y = y, x = x, w1 = w1, correction = correction, true_range = true_range, max_try = max_try, w3 = w3, w2 = w2, no_tf = no_tf, model = model, rcpp = rcpp, cma_pop_multi = cma_pop_multi)
    if (nth_res$solve) break
    if (i == max_try) message("Failed to solve non-threshold model")
  }
  p = length(w1[1,])
  y = df_data(input_order(y))
  x = df_data(input_order(x))
  tp = dim(y)[1]
  t = tp / p
  if (t < 3) stop("Time period less than 3") 
  t1 = t-1
  tp1 = t1*p
  y_1 = as.matrix(y[-c((tp1+1):tp),-c(1,2)])
  y = as.matrix(y[-c(1:p),-c(1,2)])
  x = as.matrix(x[,-c(1,2)])
  x = x[-c(1:p),,drop = F]
  mat_c = diag(2, t1, t1)
  mat_c[(row(mat_c)==col(mat_c)+1)|(row(mat_c)==col(mat_c)-1)] = -1
  inv_c = solve(mat_c)
  x_t1 = x
  x_t2 = x
  x_t1[!th,] = 0
  x_t2[th,] = 0
  if (no_tf){
    x_mt = cbind(x_t1, x_t2)
    k = dim(x_mt)[2]
    k_o = dim(x_mt)[2]
  } else {
    #individual (time)
    tif = as.matrix(rep(1, p))
    tif = kronecker(diag(t1), tif)
    # tif = tif[, -(length(tif[1, ]))]
    x_mt = do.call(cbind, list(x_t1, x_t2, tif))
    k = dim(x_mt)[2]
    k_o = dim(x_mt)[2] - dim(tif)[2]
  }
  if (true_range){
    const_rcma_r = function(para, w, th, t, w_er, w_lam){
      rho1 = para[1]
      alp1 = para[2]
      theta1 = para[3]
      lam1 = para[4]
      rho2 = para[5]
      alp2 = para[6]
      theta2 = para[7]
      lam2 = para[8]
      p = dim(w)[1]
      iaws = diag(p) - merge_th_diag(alp1, alp2, th) %*% w_er
      irws = diag(p) - merge_th_diag(rho1, rho2, th) %*% w
      iirws = solve(irws)
      lws = merge_th_diag(lam1, lam2, th) %*% w_lam
      dthetas = merge_th_diag(theta1, theta2, th)
      if (is.singular.matrix(iaws)| norm((dthetas + lws)%*%iirws) >= 1){
        return(F)
      } else{
        return(T)
      }
    }
    const_rcma = function(x) {const_rcma_r(para = x, w = w1, t = t1, th = th, w_er = w3, w_lam = w2)}
  } else {
    alp_lbound = min(Re(eigen(w3)$values[abs(Im(eigen(w3)$values)) < 1e-6]))
    const_rcma_r = function(para, alp_lbound){
      rho1 = para[1]
      alp1 = para[2]
      theta1 = para[3]
      lam1 = para[4]
      rho2 = para[5]
      alp2 = para[6]
      theta2 = para[7]
      lam2 = para[8]
      if (alp1 >= 1| alp1 <= 1/alp_lbound| alp2 >= 1| alp2 <= 1/alp_lbound| abs(theta1) + abs(rho1) + abs(lam1) >= 1|abs(theta2) + abs(rho2) + abs(lam2) >= 1){
        return(F)
      } else{
        return(T)
      }
    }
    const_rcma = function(x) {const_rcma_r(para = x, alp_lbound = alp_lbound)}
  }
  switch (model,
          "full" = {
            init_value = rep(c(nth_res$coefficient$lambda1, nth_res$coefficient$lambda3, nth_res$coefficient$rho, nth_res$coefficient$lambda2), 2)
            if(rcpp){
              objfun_rcma = function(par) {msdpdth_aqs(para = par, y = y, x_ = x_mt, w = w1, th_e = rep(th, t1)+0, y1 = y_1, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2)}
            }else{
              objfun_rcma = function(par) {th_full_aqs(para = par, y = y, x_ = x_mt, w = w1, th = th, y1 = y_1, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2)}
            }
            for (i in 1:max_try){
              cma_obj = cmaNew()
              cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
              #with non-threshold as initial value
              cmaInit(cma_obj, dimension = 8, initialX = init_value)
              optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma, maxDimPrint = 8)
              if (optim_res$bestFitness < 1e-12) break
            }
            output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
            output$coefficient = list(lambda1_1 = optim_res$bestX[1],
                                      lambda3_1 = optim_res$bestX[2],
                                      rho_1 = optim_res$bestX[3],
                                      lambda2_1 = optim_res$bestX[4],
                                      lambda1_2 = optim_res$bestX[5],
                                      lambda3_2 = optim_res$bestX[6],
                                      rho_2 = optim_res$bestX[7],
                                      lambda2_2 = optim_res$bestX[8])
            output$coefficient = c(output$coefficient, th_full_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "beta_sigs",  w_er = w3, w_lam = w2, inv_c = inv_c))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (all_er){
              all_se = th_full_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er,  w_er = w3, w_lam = w2, inv_c = inv_c)
              if (correction){
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
                output$coefficient[["lambda3_1_se"]] = vc_er[k+5]
                output$coefficient[["rho_1_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_1_se"]] = vc_er[k+4]
                output$coefficient[["lambda1_2_se"]] = vc_er[k+7]
                output$coefficient[["lambda3_2_se"]] = vc_er[k+9]
                output$coefficient[["rho_2_se"]] = vc_er[k+6]
                output$coefficient[["lambda2_2_se"]] = vc_er[k+8]
                if (no_tf){
                  output[["vc_mat"]] = all_se$vc_mat
                } else {
                  output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
                gamma_er = sqrt(diag(all_se$gamma_mat))
                output$coefficient[["beta_se_gamma"]] = gamma_er[1:k_o]
                output$coefficient[["sigma2_se_gamma"]] = gamma_er[k+1]
                output$coefficient[["lambda1_1_se_gamma"]] = gamma_er[k+3]
                output$coefficient[["lambda3_1_se_gamma"]] = gamma_er[k+5]
                output$coefficient[["rho_1_se_gamma"]] = gamma_er[k+2]
                output$coefficient[["lambda2_1_se_gamma"]] = gamma_er[k+4]
                output$coefficient[["lambda1_2_se_gamma"]] = gamma_er[k+7]
                output$coefficient[["lambda3_2_se_gamma"]] = gamma_er[k+9]
                output$coefficient[["rho_2_se_gamma"]] = gamma_er[k+6]
                output$coefficient[["lambda2_2_se_gamma"]] = gamma_er[k+8]
                if (no_tf){
                  output[["gamma_mat"]] = all_se$gamma_mat
                } else {
                  output[["gamma_mat"]] = all_se$gamma_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
              hes_er = sqrt(diag(all_se$hes_mat))
              output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
              output$coefficient[["lambda1_1_se_hes"]] = hes_er[k+3]
              output$coefficient[["lambda3_1_se_hes"]] = hes_er[k+5]
              output$coefficient[["rho_1_se_hes"]] = hes_er[k+2]
              output$coefficient[["lambda2_1_se_hes"]] = hes_er[k+4]
              output$coefficient[["lambda1_2_se_hes"]] = hes_er[k+7]
              output$coefficient[["lambda3_2_se_hes"]] = hes_er[k+9]
              output$coefficient[["rho_2_se_hes"]] = hes_er[k+6]
              output$coefficient[["lambda2_2_se_hes"]] = hes_er[k+8]
              if (no_tf){
                output[["hes_mat"]] = all_se$hes_mat
              } else {
                output[["hes_mat"]] = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }else{
              all_se = th_full_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er,  w_er = w3, w_lam = w2, inv_c = inv_c)
              vc_er = sqrt(diag(all_se$vc_mat))
              output$coefficient[["beta_se"]] = vc_er[1:k_o]
              output$coefficient[["sigma2_se"]] = vc_er[k+1]
              output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
              output$coefficient[["lambda3_1_se"]] = vc_er[k+5]
              output$coefficient[["rho_1_se"]] = vc_er[k+2]
              output$coefficient[["lambda2_1_se"]] = vc_er[k+4]
              output$coefficient[["lambda1_2_se"]] = vc_er[k+7]
              output$coefficient[["lambda3_2_se"]] = vc_er[k+9]
              output$coefficient[["rho_2_se"]] = vc_er[k+6]
              output$coefficient[["lambda2_2_se"]] = vc_er[k+8]
              if (no_tf){
                output[["vc_mat"]] = all_se$vc_mat
              } else {
                output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
            if (residual) output$coefficient[["residuals"]] = th_full_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th,correction = correction, mode = "residual",  all_er = all_er, w_er = w3, w_lam = w2, inv_c = inv_c)
          },
          "slm" = {
            init_value = rep(c(nth_res$coefficient$lambda1, nth_res$coefficient$rho), 2)
            if(rcpp){
              objfun_rcma = function(par) {
                pars = numeric(8)
                pars[1] = par[1]
                pars[3] = par[2]
                pars[5] = par[3]
                pars[7] = par[4]
                msdpdth_aqs(para = pars, y = y, x_ = x_mt, w = w1, y1 = y_1, inv_c = inv_c, correction = correction, w_er = matrix(0,p,p), w_lam =  matrix(0,p,p), th_e = rep(th, t1)+0)
              }
            }else{
              objfun_rcma = function(par) {th_slm_aqs(para = par, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, inv_c = inv_c)}
            }
            const_rcma_slm = function(x){
              pars = numeric(8)
              pars[1] = x[1]
              pars[3] = x[2]
              pars[5] = x[3]
              pars[7] = x[4]
              return(const_rcma(pars)) 
            }
            for (i in 1:max_try){
              cma_obj = cmaNew()
              cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
              #with non-threshold as initial value
              cmaInit(cma_obj, dimension = 4, initialX = init_value)
              optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_slm, maxDimPrint = 4)
              if (optim_res$bestFitness < 1e-12) break
            }
            output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
            output$coefficient = list(lambda1_1 = optim_res$bestX[1],
                                      rho_1 = optim_res$bestX[2],
                                      lambda1_2 = optim_res$bestX[3],
                                      rho_2 = optim_res$bestX[4]
            )
            output$coefficient = c(output$coefficient, th_slm_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "beta_sigs", inv_c = inv_c))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            
            if (all_er){
              all_se = th_slm_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er, inv_c = inv_c)
              if (correction){
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
                output$coefficient[["rho_1_se"]] = vc_er[k+2]
                output$coefficient[["lambda1_2_se"]] = vc_er[k+5]
                output$coefficient[["rho_2_se"]] = vc_er[k+4]
                if (no_tf){
                  output[["vc_mat"]] = all_se$vc_mat
                } else {
                  output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
                gamma_er = sqrt(diag(all_se$gamma_mat))
                output$coefficient[["beta_se_gamma"]] = gamma_er[1:k_o]
                output$coefficient[["sigma2_se_gamma"]] = gamma_er[k+1]
                output$coefficient[["lambda1_1_se_gamma"]] = gamma_er[k+3]
                output$coefficient[["rho_1_se_gamma"]] = gamma_er[k+2]
                output$coefficient[["lambda1_2_se_gamma"]] = gamma_er[k+5]
                output$coefficient[["rho_2_se_gamma"]] = gamma_er[k+4]
                if (no_tf){
                  output[["gamma_mat"]] = all_se$gamma_mat
                } else {
                  output[["gamma_mat"]] = all_se$gamma_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
              hes_er = sqrt(diag(all_se$hes_mat))
              output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
              output$coefficient[["lambda1_1_se_hes"]] = hes_er[k+3]
              output$coefficient[["rho_1_se_hes"]] = hes_er[k+2]
              output$coefficient[["lambda1_2_se_hes"]] = hes_er[k+5]
              output$coefficient[["rho_2_se_hes"]] = hes_er[k+4]
              if (no_tf){
                output[["hes_mat"]] = all_se$hes_mat
              } else {
                output[["hes_mat"]] = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }else{
              all_se = th_slm_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er, inv_c = inv_c)
              vc_er = sqrt(diag(all_se$vc_mat))
              output$coefficient[["beta_se"]] = vc_er[1:k_o]
              output$coefficient[["sigma2_se"]] = vc_er[k+1]
              output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
              output$coefficient[["rho_1_se"]] = vc_er[k+2]
              output$coefficient[["lambda1_2_se"]] = vc_er[k+5]
              output$coefficient[["rho_2_se"]] = vc_er[k+4]
              if (no_tf){
                output[["vc_mat"]] = all_se$vc_mat
              } else {
                output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
            if (residual) output$coefficient[["residuals"]] = th_slm_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "residual",  all_er = all_er, inv_c = inv_c)
          },
          "sem" = {
            init_value = rep(c(nth_res$coefficient$lambda3, nth_res$coefficient$rho), 2)
            if(rcpp){
              objfun_rcma = function(par) {
                pars = numeric(8)
                pars[2] = par[1]
                pars[3] = par[2]
                pars[6] = par[3]
                pars[7] = par[4]
                msdpdth_aqs(para = pars, y = y, x_ = x_mt, y1 = y_1, inv_c = inv_c, correction = correction, w = matrix(0,p,p), w_er = w3, w_lam = matrix(0,p,p), th_e = rep(th, t1)+0)
              }
            }else{
              objfun_rcma = function(par) {th_sem_aqs(para = par, y = y, x_ = x_mt, y1 = y_1, th = th, inv_c = inv_c, correction = correction, w_er = w3)}
            }
            const_rcma_sem = function(x){
              pars = numeric(8)
              pars[2] = x[1]
              pars[3] = x[2]
              pars[6] = x[3]
              pars[7] = x[4]
              return(const_rcma(pars)) 
            }
            for (i in 1:max_try){
              cma_obj = cmaNew()
              cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
              #with non-threshold as initial value
              cmaInit(cma_obj, dimension = 4, initialX = init_value)
              optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_sem, maxDimPrint = 4)
              if (optim_res$bestFitness < 1e-12) break
            }
            output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
            output$coefficient = list(lambda3_1 = optim_res$bestX[1],
                                      rho_1 = optim_res$bestX[2],
                                      lambda3_2 = optim_res$bestX[3],
                                      rho_2 = optim_res$bestX[4])
            output$coefficient = c(output$coefficient, th_sem_aqs(para = optim_res$bestX, y = y, x_ = x_mt, y1 = y_1, th = th, correction = correction, mode = "beta_sigs", w_er = w3, inv_c = inv_c))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (all_er){
              all_se = th_sem_aqs(para = optim_res$bestX, y = y, x_ = x_mt, y1 = y_1, th = th,  correction = correction, mode = "opmd",  all_er = all_er, w_er = w3, inv_c = inv_c)
              if (correction){
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda3_1_se"]] = vc_er[k+3]
                output$coefficient[["rho_1_se"]] = vc_er[k+2]
                output$coefficient[["lambda3_2_se"]] = vc_er[k+5]
                output$coefficient[["rho_2_se"]] = vc_er[k+4]
                if (no_tf){
                  output[["vc_mat"]] = all_se$vc_mat
                } else {
                  output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
                gamma_er = sqrt(diag(all_se$gamma_mat))
                output$coefficient[["beta_se_gamma"]] = gamma_er[1:k_o]
                output$coefficient[["sigma2_se_gamma"]] = gamma_er[k+1]
                output$coefficient[["lambda3_1_se_gamma"]] = gamma_er[k+3]
                output$coefficient[["rho_1_se_gamma"]] = gamma_er[k+2]
                output$coefficient[["lambda3_2_se_gamma"]] = gamma_er[k+5]
                output$coefficient[["rho_2_se_gamma"]] = gamma_er[k+4]
                if (no_tf){
                  output[["gamma_mat"]] = all_se$gamma_mat
                } else {
                  output[["gamma_mat"]] = all_se$gamma_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
              hes_er = sqrt(diag(all_se$hes_mat))
              output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
              output$coefficient[["lambda3_1_se_hes"]] = hes_er[k+3]
              output$coefficient[["rho_1_se_hes"]] = hes_er[k+2]
              output$coefficient[["lambda3_2_se_hes"]] = hes_er[k+5]
              output$coefficient[["rho_2_se_hes"]] = hes_er[k+4]
              if (no_tf){
                output[["hes_mat"]] = all_se$hes_mat
              } else {
                output[["hes_mat"]] = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }else{
              all_se = th_sem_aqs(para = optim_res$bestX, y = y, x_ = x_mt, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er, w_er = w3, inv_c = inv_c)
              vc_er = sqrt(diag(all_se$vc_mat))
              output$coefficient[["beta_se"]] = vc_er[1:k_o]
              output$coefficient[["sigma2_se"]] = vc_er[k+1]
              output$coefficient[["lambda3_1_se"]] = vc_er[k+3]
              output$coefficient[["rho_1_se"]] = vc_er[k+2]
              output$coefficient[["lambda3_2_se"]] = vc_er[k+5]
              output$coefficient[["rho_2_se"]] = vc_er[k+4]
              if (no_tf){
                output[["vc_mat"]] = all_se$vc_mat
              } else {
                output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
            if (residual) output$coefficient[["residuals"]] = th_sem_aqs(para = optim_res$bestX, y = y, x_ = x_mt, y1 = y_1, th = th, correction = correction, mode = "residual",  all_er = all_er, w_er = w3, inv_c = inv_c)
          },
          "sltl" = {
            init_value = rep(c(nth_res$coefficient$lambda1, nth_res$coefficient$rho, nth_res$coefficient$lambda2), 2)
            if(rcpp){
              objfun_rcma = function(par) {
                pars = numeric(8)
                pars[1] = par[1]
                pars[3] = par[2]
                pars[4] = par[3]
                pars[5] = par[4]
                pars[7] = par[5]
                pars[8] = par[6]
                msdpdth_aqs(para = pars, y = y, x_ = x_mt, w = w1, y1 = y_1, th_e = rep(th, t1)+0, inv_c = inv_c, correction = correction, w_er = matrix(0, p, p), w_lam = w2)
              }
            }else{
              objfun_rcma = function(par) {th_sltl_aqs(para = par, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, inv_c = inv_c, correction = correction, w_lam = w2)}
            }
            const_rcma_sltl = function(x){
              pars = numeric(8)
              pars[1] = x[1]
              pars[3] = x[2]
              pars[4] = x[3]
              pars[5] = x[4]
              pars[7] = x[5]
              pars[8] = x[6]
              return(const_rcma(pars)) 
            }
            for (i in 1:max_try){
              cma_obj = cmaNew()
              cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
              #with non-threshold as initial value
              cmaInit(cma_obj, dimension = 6, initialX = init_value)
              optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_sltl, maxDimPrint = 6)
              if (optim_res$bestFitness < 1e-12) break
            }
            output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
            output$coefficient = list(lambda1_1 = optim_res$bestX[1],
                                      rho_1 = optim_res$bestX[2],
                                      lambda2_1 = optim_res$bestX[3],
                                      lambda1_2 = optim_res$bestX[4],
                                      rho_2 = optim_res$bestX[5],
                                      lambda2_2 = optim_res$bestX[6])
            output$coefficient = c(output$coefficient, th_sltl_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "beta_sigs", w_lam = w2, inv_c = inv_c))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (all_er){
              all_se = th_sltl_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er, w_lam = w2, inv_c = inv_c)
              if (correction){
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
                output$coefficient[["rho_1_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_1_se"]] = vc_er[k+4]
                output$coefficient[["lambda1_2_se"]] = vc_er[k+6]
                output$coefficient[["rho_2_se"]] = vc_er[k+5]
                output$coefficient[["lambda2_2_se"]] = vc_er[k+7]
                if (no_tf){
                  output[["vc_mat"]] = all_se$vc_mat
                } else {
                  output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
                gamma_er = sqrt(diag(all_se$gamma_er_mat))
                output$coefficient[["beta_se_gamma"]] = gamma_er[1:k_o]
                output$coefficient[["sigma2_se_gamma"]] = gamma_er[k+1]
                output$coefficient[["lambda1_1_se_gamma"]] = gamma_er[k+3]
                output$coefficient[["rho_1_se_gamma"]] = gamma_er[k+2]
                output$coefficient[["lambda2_1_se_gamma"]] = gamma_er[k+4]
                output$coefficient[["lambda1_2_se_gamma"]] = gamma_er[k+6]
                output$coefficient[["rho_2_se_gamma"]] = gamma_er[k+5]
                output$coefficient[["lambda2_2_se_gamma"]] = gamma_er[k+7]
                if (no_tf){
                  output[["gamma_mat"]] = all_se$gamma_mat
                } else {
                  output[["gamma_mat"]] = all_se$gamma_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
              hes_er = sqrt(diag(all_se$hes_mat))
              output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
              output$coefficient[["lambda1_1_se_hes"]] = hes_er[k+3]
              output$coefficient[["rho_1_se_hes"]] = hes_er[k+2]
              output$coefficient[["lambda2_1_se_hes"]] = hes_er[k+4]
              output$coefficient[["lambda1_2_se_hes"]] = hes_er[k+6]
              output$coefficient[["rho_2_se_hes"]] = hes_er[k+5]
              output$coefficient[["lambda2_2_se_hes"]] = hes_er[k+7]
              if (no_tf){
                output[["hes_mat"]] = all_se$hes_mat
              } else {
                output[["hes_mat"]] = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }else{
              all_se = th_sltl_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "opmd",  all_er = all_er, w_lam = w2, inv_c = inv_c)
              vc_er = sqrt(diag(all_se$vc_mat))
              output$coefficient[["beta_se"]] = vc_er[1:k_o]
              output$coefficient[["sigma2_se"]] = vc_er[k+1]
              output$coefficient[["lambda1_1_se"]] = vc_er[k+3]
              output$coefficient[["rho_1_se"]] = vc_er[k+2]
              output$coefficient[["lambda2_1_se"]] = vc_er[k+4]
              output$coefficient[["lambda1_2_se"]] = vc_er[k+6]
              output$coefficient[["rho_2_se"]] = vc_er[k+5]
              output$coefficient[["lambda2_2_se"]] = vc_er[k+7]
              if (no_tf){
                output[["vc_mat"]] = all_se$vc_mat
              } else {
                output[["vc_mat"]] = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
            if (residual) output$residual = th_sltl_aqs(para = optim_res$bestX, y = y, x_ = x_mt, w = w1, y1 = y_1, th = th, correction = correction, mode = "residual",  all_er = all_er, w_lam = w2, inv_c = inv_c)
          },
          stop("Undefined model")
  )
  output$model = model
  class(output) = "msdpdth"
  return(output)
}