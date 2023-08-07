
# Function to return predicted values across the x range of the observed data
# for brmsfit objects with a grouping variable.
# If by_group = FALSE returns predicted values marginalized across the groups, 
# if by_group = TRUE predictions are returned for each group.
pred_out <- function(x, x_var, group_var, 
                     by_group = TRUE,
                     probs = c(0.025, 0.5, 0.975),
                     precision = 1000){
  dat_list <- list(seq(min(x$data[x_var]), 
                       max(x$data[x_var]), 
                       length=precision),
                   unlist(unique(x$data[group_var]))) 
  names(dat_list) <- c(x_var, group_var)
  pred_dat <- expand.grid(dat_list)
  
  if(by_group){
    dat_out <- 
      cbind(pred_dat,
            t(apply(posterior_epred(x, newdata = pred_dat), MARGIN = 2, 
                    FUN = quantile, probs = probs))) |>
      data.frame()   
  } else {
    x_vec <- unique(pred_dat$log_x)
    post_mat <- posterior_epred(x, newdata = pred_dat)
    dat_out <-
      cbind(x_vec, lapply(x_vec, FUN = function(v){
        index_v <- which(v==pred_dat[, x_var])
        quantile(unlist(post_mat[, index_v]), probs = probs)
      }) |> bind_rows() |> data.frame())
    colnames(dat_out) <- c( x_var, colnames(dat_out)[-1])
  }
  
  return(dat_out)
}

# Function to extract NEC from a fitted model containing an NEC parameter,
# for brmsfit objects with a grouping variable.
# If by_group = FALSE returns posterior values marginalized across the groups, 
# if by_group = TRUE posterior values are returned for each group.
nec.brmsfit <- function(x, x_var, group_var, 
         by_group = TRUE,
         probs = c(0.025, 0.5, 0.975),
         posterior = FALSE){
  dat_list <- list(median(x$data[[x_var]]),
                   unique(x$data[[group_var]]))
  names(dat_list) <- c(x_var, group_var)
  pred_dat <- expand.grid(dat_list)
  
  if(by_group & !posterior){
    nec_out <- 
      cbind(pred_dat,
            t(apply(posterior_epred(x, newdata = pred_dat, nlpar = "nec"), MARGIN = 2, 
                    FUN = quantile, probs = probs))) |>
      data.frame()
  }

  if(by_group & posterior){
    post_mat <- posterior_epred(x, newdata = pred_dat, nlpar = "nec")
    colnames(post_mat) <- unique(x$data[[group_var]])
    nec_out <- post_mat  |>  
      data.frame() |> 
      pivot_longer(everything(), names_to = group_var, values_to = "NEC")
  }
  
  if(!by_group & !posterior){
    nec_out <- quantile(unlist(posterior_epred(x, newdata = pred_dat, nlpar = "nec")),
                        probs = probs)
  }
  
  if(!by_group & posterior){
    nec_out <- unlist(posterior_epred(x, newdata = pred_dat, nlpar = "nec"))
  }

  return(nec_out)
}











