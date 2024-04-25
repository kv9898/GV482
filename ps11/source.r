coefplot_prms = function(object, ..., sd, ci_low, ci_high, x, x.shift = 0, dict,
                         keep, drop, order, ci_level = 0.95, df.t = NULL, ref = "auto",
                         only.i = TRUE, sep, as.multiple = FALSE){
  
  # get the default for:
  # dict, ci.level, ref
  
  dots = list(...)
  is_internal = isTRUE(dots$internal__)
  if("i.select" %in% names(dots)){
    i.select = dots$i.select
    dots$i.select = NULL
  }
  
  
  varlist = list(ci = paste0(ci_level * 100, "%"))
  dots_drop = c()
  
  suggest_ref_line = FALSE
  multiple_est = FALSE
  NO_NAMES = FALSE
  is_iplot = only.i
  
  set_up(1 + is_iplot)
  
  AXIS_AS_NUM = FALSE
  if(is_internal == FALSE && ((is.list(object) && class(object)[1] == "list") || "fixest_multi" %in% class(object))){
    # This is a list of estimations
    
    #
    # Multiple estimations ####
    #
    
    multiple_est = TRUE
    
    mc = match.call()
    mc$object = as.name("my__object__")
    mc$only.params = TRUE
    mc$internal__ = TRUE
    
    nb_est = length(object)
    
    rerun = FALSE
    first = TRUE
    
    while(first || rerun){
      first = FALSE
      
      res = varlist = list()
      all_inter = c()
      all_inter_root = c()
      dots_drop = c()
      num_axis = TRUE
      suggest_ref_line = TRUE
      for(i in 1:nb_est){
        # cat("Eval =", i, "\n")
        my__object__ = object[[i]]
        prms = try(eval(mc), silent = TRUE)
        
        if("try-error" %in% class(prms)){
          stop_up("The {nth ? i} element of 'object' raises and error:\n", prms)
        }
        
        # Some meta variables
        varlist$ci = unique(prms$varlist$ci)
        varlist$depvar = unique(prms$varlist$depvar)
        varlist$var = unique(prms$varlist$var)
        varlist$fe = unique(prms$varlist$fe)
        varlist$i = unique(prms$varlist$i)
        
        dots_drop = unique(c(dots_drop, prms$dots_drop))
        suggest_ref_line = suggest_ref_line && prms$suggest_ref_line
        
        num_axis = num_axis && prms$num_axis
        
        df_prms = prms$prms
        df_prms$est_nb = i
        res[[i]] = df_prms
      }
    }
    
    AXIS_AS_NUM = num_axis
    
    all_estimates = do.call("rbind", res)
    
    # We respect the order provided by the user
    my_names_order = unique(all_estimates$estimate_names)
    my_order = 1:length(my_names_order)
    names(my_order) = my_names_order
    all_estimates$id_order = my_order[as.character(all_estimates$estimate_names)]
    all_estimates = all_estimates[base::order(all_estimates$id_order, all_estimates$est_nb), ]
    
    # we rescale -- beware of multiple est whose x is numeric!!!
    
    if(nb_est > 1){
      if(missing(sep)){
        all_sep = c(0.2, 0.2, 0.18, 0.16)
        if(length(all_sep) < nb_est - 1){
          sep = 1 / (nb_est-1) * 0.7
        } else {
          sep = all_sep[nb_est - 1]
        }
        
      } else {
        n_sep = length(sep)
        if(n_sep > 1){
          if(n_sep < nb_est - 1){
            sep = sep[n_sep]
          } else {
            sep = sep[nb_est - 1]
          }
        }
      }
      
      if(AXIS_AS_NUM){
        all_estimates$x_new = all_estimates$x + ((all_estimates$est_nb - 1) / (nb_est - 1) - 0.5) * ((nb_est - 1) * sep)
      } else {
        all_estimates$x_new = all_estimates$id_order + ((all_estimates$est_nb - 1) / (nb_est - 1) - 0.5) * ((nb_est - 1) * sep)
      }
      
    } else {
      sep = 0
      all_estimates$x_new = all_estimates$id_order
    }
    
    # The coefficients
    
    estimate = all_estimates$y
    names(estimate) = all_estimates$estimate_names
    ci_low = all_estimates$ci_low
    ci_high = all_estimates$ci_high
    x = all_estimates$x_new
    
    prms = data.frame(estimate = estimate, ci_low = ci_low, ci_high = ci_high, 
                      x = x, id = all_estimates$est_nb, 
                      estimate_names = all_estimates$estimate_names, 
                      estimate_names_raw = all_estimates$estimate_names_raw)
    if(!is.null(all_estimates$is_ref)){
      prms$is_ref = all_estimates$is_ref
    }
    
  } else {
    
    #
    # Single estimation ####
    #
    
    
    if(is.list(object)){
      
      sum_exists = FALSE
      for(c_name in class(object)){
        if(exists(paste0("summary.", c_name), mode = "function")){
          sum_exists = TRUE
          break
        }
      }
      
      if(!sum_exists){
        # stop("There is no summary method for objects of class ", c_name, ". 'coefplot' applies summary to the object to extract the coeftable. Maybe add directly the coeftable in object instead?")
        mat = coeftable(object)
      } else {
        fun_name = paste0("summary.", c_name)
        args_name_sum = names(formals(fun_name))
        args_sum = intersect(names(dots), args_name_sum)
        
        dots_drop = args_sum
        
        # we kick out the summary arguments from dots
        dots[args_sum] = NULL
        
        # We reconstruct a call to coeftable
        mc_coeftable = match.call(expand.dots = TRUE)
        mc_coeftable[[1]] = as.name("coeftable")
        mc_coeftable[setdiff(names(mc_coeftable), c(args_sum, "object", ""))] = NULL
        
        mat = eval(mc_coeftable, parent.frame())
      }
      
      sd = mat[, 2]
      estimate = mat[, 1]
      
      names(estimate) = rownames(mat)
      
      if("fml" %in% names(object)){
        depvar = gsub(" ", "", as.character(object$fml)[[2]])
        if(depvar %in% names(dict)) depvar = dict[depvar]
        varlist$depvar = depvar
      }
      
    } else if(is.matrix(object)){
      # object is a matrix containing the coefs and SEs
      
      m_names = tolower(colnames(object))
      if(ncol(object) == 4 || (grepl("estimate", m_names[1]) && grepl("std\\.? error", m_names[1]))){
        sd = object[, 2]
        estimate = object[, 1]
        
        names(estimate) = rownames(object)
        
      } else {
        stop_up("Argument 'object' is a matrix but it should contain 4 columns (the two first ones should be reporting the estimate and the standard-error). Either provide an appropriate matrix or give directly the vector of estimated coefficients in arg. estimate.")
      }
      
    } else if(length(object[1]) > 1 || !is.null(dim(object)) || !is.numeric(object)){
      stop_up("Argument 'object' must be either: i) an estimation object, ii) a matrix of coefficients table, or iii) a numeric vector of the point estimates. Currently it is neither of the three.")
    } else {
      # it's a numeric vector
      estimate = object
    }
    
    n = length(estimate)
    
    if(missing(sd)){
      if(missing(ci_low) || missing(ci_high)) stop_up("If 'sd' is not provided, you must provide the arguments 'ci_low' and 'ci_high'.")
      
      varlist$ci = NULL
    } else {
      if(!missing(ci_low) || !missing(ci_high)) warning("Since 'sd' is provided, arguments 'ci_low' or 'ci_high' are ignored.")
      
      # We compute the CI
      check_arg(df.t, "NULL | numeric scalar GE{1}")
      if(inherits(object, "fixest")){
        if(is.null(df.t)){
          df.t = degrees_freedom(object, "t")
          if(is.null(df.t)){
            # may for user defined functions using fixest
            df.t = degrees_freedom(object, "resid")
            if(is.null(df.t)){
              # let's be safe
              df.t = Inf
            }
          }
        }
        
        nb = fixest_CI_factor(object, ci_level, df.t = df.t)[2]
      } else {
        # non fixest objects
        if(is.null(df.t)){
          df.t = tryCatch(nobs(object) - length(coef(object)), 
                          error = function(e) NULL)
          if(is.null(df.t)){
            if(not_too_many_messages("coefplot_df")){
              message("The degrees of freedom for the t distribution could",
                      " not be deduced. Using a Normal distribution instead.\n",
                      "Note that you can provide the argument `df.t` directly.")
            }
            df.t = Inf
          }
        }
        
        nb = abs( qt( (1-ci_level) / 2, df.t) )
      }
      
      ci_high = estimate + nb*sd
      ci_low = estimate - nb*sd
    }
    
    #
    # iplot ####
    #
    
    ref_id = NA
    xlab_suggest = NULL
    if(is_iplot){
      if(is.null(names(estimate))){
        stop_up("'iplot' must be used only with fixest objects containing variables created with i(). Currently it does not seem to be the case.")
      }
      
      all_vars = names(estimate)
      
      if(!any(grepl("::", all_vars))){
        stop_up("'iplot' must be used only with fixest objects containing variables created with i(). Currently it does not seem to be the case.")
      }
      
      # Four cases:
      # - factor_var::value
      # - factor_var::value:xnum
      # - xnum:factor_var::value
      # - factor_var::value:xfact::value
      
      # Restriction:
      # it only accepts "pure" i() variables
      
      # We can handle only case 1, 2, and 3
      # case 4 is too messy
      # case 4: multiple lines + legend ?
      
      # We take the first i() in the list
      # after having applied keep_apply
      
      # avoids bug with IVs => problem if user names the variables that way
      is_IV = FALSE
      if(isTRUE(object$iv) && identical(object$iv_stage, 2)){
        all_vars = gsub("^fit_", "", all_vars)
        names(estimate) = all_vars
      }
      
      all_vars = keep_apply(all_vars, keep)
      
      mm_info = object$model_matrix_info
      ANY_TWO_FACTORS = FALSE
      
      # Finding out which to display
      i_running = 0
      for(i in seq_along(mm_info)){
        info = mm_info[[i]]
        if(isFALSE(info$is_inter_fact)){
          
          # First case: regular variable. species::setosa
          if(any(info$coef_names %in% all_vars)){
            i_running = i_running + 1
            
            if(i.select == i_running){
              # That's the one
              break
            }
            
          }
        }
      }
      
      if(i_running < i.select){
        # not found, maybe an interaction? second round
        
        for(i in seq_along(mm_info)){
          info = mm_info[[i]]
          ANY_TWO_FACTORS = ANY_TWO_FACTORS || info$is_inter_fact
          if(isFALSE(info$is_inter_fact)){
            
            # Second case: interacted variable. species::setosa:x1 ou x1:species::setosa
            if(!isTRUE(info$is_inter_num)){
              pattern = paste0("(:", info$f_name, "::.+)|(", info$f_name, "::[^:]+:)")
              vars_inter = grep(pattern, all_vars, value = TRUE)
              if(length(vars_inter) > 0){
                vars_inter_unik = unique(gsub(pattern, "", vars_inter))
                ok = FALSE
                for(v in vars_inter_unik){
                  
                  inter_before = paste0(v, ":", info$coef_names_full)
                  inter_after  = paste0(info$coef_names_full, ":", v)
                  
                  if(any(c(inter_before, inter_after) %in% all_vars)){
                    i_running = i_running + 1
                    if(i.select == i_running){
                      # That's the one
                      ok = TRUE
                      break
                    }
                  }
                }
                
                if(ok){
                  if(any(inter_before %in% all_vars)){
                    info$coef_names_full = inter_before
                  } else {
                    info$coef_names_full = inter_after
                  }
                  
                  break
                }
                
              }
            }
          }
        }
      }
      
      if(i_running < i.select){
        
        if(i_running == 0){
          if(length(mm_info) == 0){
            stop_up("In iplot(), no variable created with i() found (it works only with that kind of variables)")
          }
          
          if(ANY_TWO_FACTORS){
            stop_up("In iplot(), no valid variable created with i() found. Note that it does not work when i() interacts two factors.")
          }
          
          stop_up("In iplot(), no valid variable created with i() found. Note that if there are interactions, they should be created with something of the form i(factor_var, num_var), otherwise the support is limited.")
          
        }
        
        stop_up("In iplot(), the value of 'i.select' (=", i.select, ") exceeds the number of valid i() variables found (=", i_running, "). Please give a lower number.")
      }
      
      # "keep" here works differently => new arg. i.select?
      
      ANY_AUTO_REF = length(info$ref_id) > 0
      SHOW_REF = (identical(ref, "auto") || isTRUE(ref) || identical(ref, "all")) && ANY_AUTO_REF
      SHOW_REF_FIRST = SHOW_REF && !identical(ref, "all")
      
      # Global variables for fill_coef function
      names_coef = names(estimate)
      names_all = info$coef_names_full
      new_names = info$items
      
      is_rm = FALSE
      ID_rm = NULL
      if(ANY_AUTO_REF){
        if(SHOW_REF_FIRST){
          ID_rm = info$ref_id[-1]
        } else if(SHOW_REF == FALSE){
          ID_rm = info$ref_id
        }
        is_rm = length(ID_rm) > 0
      }
      
      if(is_rm){
        names_all = names_all[-ID_rm]
        new_names = new_names[-ID_rm]
      }
      
      fill_coef = function(coef){
        # we get the vector right
        res = rep(0, length(names_all))
        names(res) = names_all
        
        # we need it for CI
        names(coef) = names_coef
        inter_names = intersect(names_all, names_coef)
        res[inter_names] = coef[inter_names]
        
        names(res) = new_names
        res
      }
      
      estimate = fill_coef(estimate)
      ci_high = fill_coef(ci_high)
      ci_low = fill_coef(ci_low)
      estimate_names = new_names
      estimate_names_raw = names_all
      
      if(isTRUE(info$is_num) && missing(x)){
        AXIS_AS_NUM = TRUE
        names(estimate) = NULL
        x = info$items
        if(is_rm) x = x[-ID_rm]
      }
      
      # ref
      if(SHOW_REF){
        suggest_ref_line = length(info$ref_id) == 1 && info$is_inter_num
        is_ref = seq_along(estimate) == info$ref_id[1]
        
      } else {
        is_ref = rep(FALSE, length(estimate))
      }
      
      n = length(estimate)
      
      varlist$i = dict_apply(gsub("::.*", "", names_all[1]), dict)
    }
    
    # The DF of all the parameters
    prms = data.frame(estimate = estimate, ci_low = ci_low, ci_high = ci_high)
    if(is_iplot){
      prms$estimate_names = estimate_names
      prms$estimate_names_raw = estimate_names_raw
      prms$is_ref = is_ref
    } else if(!is.null(names(estimate))){
      prms$estimate_names = names(estimate)
      prms$estimate_names_raw = names(estimate)
    } else {
      NO_NAMES = TRUE
      prms$estimate_names = paste0("c", 1:nrow(prms))
      prms$estimate_names_raw = prms$estimate_names
    }
    
    # setting the names of the estimate
    if(!missing(x)){
      if(length(x) != n){
        stop_up("Argument 'x' must have the same length as the number of coefficients (currently {len ? x} vs {n ? n}).")
      }
      
      if(!is.numeric(x)){
        names(estimate) = x
        prms$estimate_names = names(estimate)
      } else if(NO_NAMES) {
        AXIS_AS_NUM = TRUE
      }
      
      prms$x = x
    }
    
    # We add the reference
    if(!(identical(ref, "auto") || identical(ref, "all")) && length(ref) > 0 && !isFALSE(ref)){
      
      if(AXIS_AS_NUM){
        if(!is.numeric(ref) || length(ref) > 1){
          check_arg(ref, "numeric scalar")
          if(is.logical(ref)){
            stop_up("In this context, argument 'ref' must be a numeric scalar.")
          }
        } else {
          names(ref) = "reference"
        }
        
      } else {
        if(is.null(names(ref))){
          if(!is.character(ref) || length(ref) > 1){
            check_arg(ref, "character scalar", .message = "Argument 'ref' must be either: a single character, either a list or a named integer vector of length 1 (The integer gives the position of the reference among the coefficients).")
          } else {
            refname = ref
            ref = list()
            ref[[refname]] = 1
          }
        }
        
        ref = unlist(ref)
        
        if(!isScalar(ref, int = TRUE)){
          reason = ifelse(length(ref) == 1, " an integer", " of length 1")
          stop_up("Argument 'ref' must be either: a single character, either a list or a named integer vector of length 1. The integer gives the position of the reference among the coefficients. Currently this is not ", reason, ".")
        }
      }
      
      # we recreate the parameters
      n = nrow(prms)
      prms$is_ref = FALSE
      ref_row = data.frame(estimate = 0, ci_low = 0, ci_high = 0, 
                           estimate_names = names(ref), 
                           estimate_names_raw = names(ref), is_ref = TRUE)
      if(AXIS_AS_NUM){
        ref_row$x = unname(ref)
        prms = rbind(prms, ref_row)
        prms = prms[base::order(prms$x), ]
        x = prms$x
        
      } else {
        prms = rbind(prms, ref_row)
        if(ref > n) ref = n + 1
        ids = 1:n
        ids[ids >= ref] = ids[ids >= ref] + 1
        prms = prms[base::order(c(ids, ref)), ]
      }
      
    }
  }
  
  n = nrow(prms)
  
  #
  # order/drop/dict ####
  #
  
  if(!is.null(prms$estimate_names)){
    
    if(!AXIS_AS_NUM){
      # dict
      if(missnull(dict)){
        dict = c()
      } else {
        prms$estimate_names = dict_apply(prms$estimate_names, dict)
      }
    }
    
    # dropping some coefs
    all_vars = unique(prms$estimate_names)
    
    if(!missing(keep) && length(keep) > 0){
      all_vars = keep_apply(all_vars, keep)
      
      if(length(all_vars) == 0){
        msg = ""
        if(is_iplot){
          msg = sma(" Didn't you mean to use 'i.select'? ",
                    "\nAlso note that in `iplot`, the variables names are the **values** that should be in the x-axis.",
                    "\nFYI, valid x-axis values are: {bq, '5|etc'K, ', 'c ? valid_vars}.",
                    valid_vars = unique(prms$estimate_names))
        }
        
        stop_up("Argument 'keep' has removed all variables!", msg)
      }
      
      prms = prms[prms$estimate_names %in% all_vars,]
    }
    
    if(!missing(drop) && length(drop) > 0){
      all_vars = drop_apply(all_vars, drop)
      
      if(length(all_vars) == 0){
        stop_up("Argument 'drop' has removed all variables!")
      }
      
      prms = prms[prms$estimate_names %in% all_vars,]
    }
    
    if(AXIS_AS_NUM && !missing(order) && length(order) > 0){
      warning("Argument 'order' is ignored since the x-axis is numeric.")
    }
    
    if(!AXIS_AS_NUM){
      # ordering the coefs
      if(!missing(order) && length(order) > 0){
        all_vars = order_apply(all_vars, order)
        
        my_order = 1:length(all_vars)
        names(my_order) = all_vars
        prms$id_order = my_order[prms$estimate_names]
        
        prms = prms[base::order(prms$id_order), ]
      }
    }
    
    estimate = prms$estimate
    names(estimate) = prms$estimate_names
    ci_high = prms$ci_high
    ci_low = prms$ci_low
    if(!is.null(prms$x)) x = prms$x
    
    n = nrow(prms)
  }
  
  #
  # we create x_labels, x_value & x_at
  #
  
  # id: used for colors/lty etc
  if(!multiple_est){
    if(as.multiple){
      prms$id = 1:nrow(prms)
    } else {
      prms$id = 1
    }
  }
  
  if(multiple_est){
    # We don't allow x.shift
    
    my_xlim = range(x) + c(-1, 1) * ((nb_est - 1) * sep)
    x_value = x
    
    # We allow identical aliases to display properly (they can be identified with 'group')
    quoi = unique(prms[, c("estimate_names", "estimate_names_raw")])
    
    x_labels = quoi$estimate_names
    x_labels_raw = quoi$estimate_names_raw
    x_at = 1:length(x_labels)
    
  } else if(!missing(x) && is.numeric(x)){
    my_xlim = range(c(x + x.shift, x - x.shift))
    
    x_value = x + x.shift
    
    if(NO_NAMES){
      x_at = NULL
      x_labels = x_labels_raw = NULL
    } else {
      x_at = x
      x_labels = prms$estimate_names
      x_labels_raw = prms$estimate_names_raw
    }
    
  } else {
    x_at = 1:n
    x_value = 1:n + x.shift
    
    if(NO_NAMES){
      x_labels = x_labels_raw = 1:n
    } else {
      x_labels = prms$estimate_names
      x_labels_raw = prms$estimate_names_raw
    }
    
    my_xlim = range(c(1:n + x.shift, 1:n - x.shift)) + c(-0.5, +0.5)
  }
  
  prms$x = x_value
  prms$y = prms$estimate
  
  return(list(prms = prms, num_axis = AXIS_AS_NUM, at = x_at, labels = x_labels, 
              x_labels_raw = x_labels_raw, varlist = varlist, 
              dots_drop = dots_drop, xlim = my_xlim, 
              suggest_ref_line = suggest_ref_line, 
              multiple_est = multiple_est))
}