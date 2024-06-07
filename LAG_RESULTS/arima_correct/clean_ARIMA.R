#' ES_ARIMA function
#' This function evaluates the performance of ensembles and AUTO ARIMA models.
#' It calculates the absolute error and weighted interval scores for each U.S. state.
#' The function fits a rolling window of 2 years of data (104 weeks) to generate forecasts.
#' The user defines if it will use an AUTO ARIMA, or ensembles of 27 or 64 models and
#' choose the number of weeks ahead for each forecast.
#' 
#' @param data_ A list containing the data for each U.S. state.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#' @param my_n_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param ES27 A logical value indicating whether to use ensembles of 27 models. Default is \code{TRUE}.
#' @param ES64 A logical value indicating whether to use ensembles of 64 models. Default is \code{FALSE}.
#'
#' @return A list containing the forecast results and performance metrics.
#'
#'
#'
ES_ARIMA<-function(data_, auto=FALSE, my_n_ahead=1, ES27=TRUE, ES64=FALSE){
  # The model run the ES27 or the ES64.
  ES27=!ES64 
  
# change name
  # Empty list of final results. It will contain forecasts, prediction intervals and number of models.
  results<-listenv()   
  if(ES27){
    pdq=c(0,1,2) # 27 possible ARIMA pdq's.
    my_order_params<-permutations(3,3,pdq, repeats.allowed = TRUE) # Create 27 permutations of 0,1,2
  }
  if(ES64){
    pdq=c(0,1,2,3) # 64 possible ARIMA pdq's.
    my_order_params<-permutations(4,3,pdq, repeats.allowed = TRUE) # Create 64 permutations of 0,1,2,3
  }
# use new name
  # Apply the ROLLING_ARIMA function to get results. Window set to 104 weeks (2 years).
  
  results[[1]]%<-% ROLLING_ARIMA(data_, my_n_ahead=my_n_ahead, look_back_amount = 104, order_params=my_order_params, auto=auto) %packages% "forecast" 
  # Resolve from system environment using the future package 
# use new name  
  suppressWarnings(invisible(resolved(results[[1]]))) # 
  
  # Put forecasts, prediction intervals and number of models into separate lists. 
  list_all_preds<-list() # Empty list for forecasts
  list_all_quantiles<- list() # Empty list for prediction intervals
  list_all_models<- list() # Empty list for number of models
# use new name
  # Get forecasts and dates
  list_all_preds[[1]]<- results[[1]][[1]][[1]]
# use new name
  # Get prediction intervals
  list_all_quantiles[[1]]<-results[[1]][[2]][[1]]
# use new name
  # Get the number of models in each ensemble by date
  list_all_models[[1]]<-results[[1]][[3]][[1]]
  ############################
  # Format forecasts for WIS #
  prediction_intervals<- FormatForScoring_correct(pred_intervals=list_all_quantiles, data_, model_name = "TestModel", my_n_week_ahead = my_n_ahead) # Format for weighted interval score calculation
  ########################
  # Format truth for WIS #
  truth<-NULL
  truth<-as.data.frame(data_) # Format flu cases from current state for WIS calculation
  truth["target_variable"]<-"cases" # rename cases as target_variable
  truth["model"]<-prediction_intervals[1,"model"]# insert 1 as the index in Model, it can be just the model name
  truth<- truth %>% rename_at("cases", ~'value') # rename the column named cases to values
  ########################################################################################
  # Calculate the WIS using our prediction intervals and truth values by target_end_date #
  my_forecast_scores<-score_forecasts(prediction_intervals, truth) 
  #############################################
  # Get the number of models in the ensembles #
  all_models<-data.frame(as.Date(unique(prediction_intervals$target_end_date)),list_all_models[[1]]) #Get the number of models for each forecasted date
  colnames(all_models)<-c("Dates","Number_of_models")
  ##############################
  # Format flu cases and dates #
  all_cases<-data.frame(data_$cases, as.Date(data_$target_end_date)) # Get flu cases from current state
  colnames(all_cases)<-c("cases","Dates")
  ##############################
  # Format forecasts and dates #
  all_forecasts<-data.frame(expm1(list_all_preds[[1]]$Prediction), as.Date(list_all_preds[[1]]$Pred_For_Date)) # Get ensembled forecast from current state
  colnames(all_forecasts)<-c("forecasts","Dates")
  ################################
  # Join flu cases and forecasts #
  forecasts_and_cases<-inner_join(all_forecasts,all_cases, by="Dates") # Join cases and forecasts
  colnames(forecasts_and_cases)<-c("forecasts","Dates","cases")
  ##############################
  # Get WIS and absolute error #
  WIS_errors<-data.frame(as.Date(c(my_forecast_scores[,"target_end_date"]$target_end_date)),c(my_forecast_scores[,"abs_error"]$abs_error),c(my_forecast_scores[,"wis"]$wis))
  colnames(WIS_errors)<-c("Dates","abs_error","WIS")  
  ##########################################################
  # Join WIS, absolute error and number of models by dates # 
  WIS_error_Nmodels<-inner_join(WIS_errors,all_models, by="Dates")
  colnames(WIS_error_Nmodels)<-c("Dates","abs_error","WIS","Number_of_models")
  ################################################################################
  # Join WIS, absolute error, number of models, forecasts and ILI cases by dates # 
  final_results<-inner_join(WIS_error_Nmodels,forecasts_and_cases, by="Dates")
  colnames(final_results)<-c("target_end_date","abs_error","WIS","Number_of_models","forecasts","cases")
  return(final_results)
}


#' ROLLING_ARIMA
#'
#' This function fits ARIMA models on a rolling window of data and generates forecasts.
#'
#' @param data_ A list containing the data for each U.S. state.
#' @param my_n_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param look_back_amount An integer specifying the number of weeks to look back for the rolling window. Default is \code{104}.
#' @param order_params A list specifying the order parameters for the ARIMA model. Default is \code{NULL}.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#'
#' @return A list containing the forecast results and performance metrics.
#'
ROLLING_ARIMA <- function(data_, my_n_ahead=1, look_back_amount = 104, order_params=NULL, auto=FALSE) {
  
  # List with the data from a single state
  single_state=list(data_) # !!!!!!! CHANGE NAME
  # All models in the ensemble
  all_my_models<-c() # !!!!!!! CHANGE NAME
  # Final predictions and dates list
  prediction<-list() # empty list for forecasts !!!!!!! CHANGE NAME
  # Predictions and dates data frame 
  prediction_df<- data.frame("Pred_For_Date"= NULL, "Prediction" = NULL) # empty data frame for forecasts and dates
  # Final predictive quantiles list
  prediction_quantile<-list() # empity list for prediction intervals !!!!!!! CHANGE NAME
  # List in which each predictive quantiles data frame will be added with forecast target_end_date as its name
  prediction_quantile_ls<- list() 

  for(iter in  1:(NROW(single_state[[1]])-(look_back_amount))){ # rolling window of 104 weeks up to the end of the ILI cases dataset
    
    # Empty list for fitted ARIMAS
    fitted_models<-list()
    # Window of 104 weeks.
    sample_data<- iter:(look_back_amount+iter-1) 
    # AIC scores for weighting the ensembles
    model_aic_scores<-c() 
    # Each model id, use later on a for loop
    model_id<-1
    
##########
# Fitting 
##########   
    
    # Start if ILI data has 104 elements
    if(length(single_state[[1]]$cases[sample_data])==look_back_amount){ 
#    if(n_unique(log(single_state[[1]]$cases[sample_data]))==look_back_amount){ 
      
      # run 1 time for the auto.arima or run 27 or 64 times for the ensembles
      for(j in 1:nrow(order_params)){
        fit<- NULL # start with fit as NULL
        
        # try to fit an ARIMA model
        tryCatch(
          expr = {
            # if auto = FALSE, run the ensemble of 27 or 64
            if(!auto){ 
              # fit ensembles of ARIMAs on log1p of the data
              fit<-arima(ts(log1p(single_state[[1]]$cases[sample_data]), deltat = 1/52), order = order_params[j,], method = "CSS-ML") 
              }
            # if auto = TRUE, run auto.arima
            else{
              # fit an auto.arima on the log1p of the data
              fit<-invisible(auto.arima(ts(log1p(single_state[[1]]$cases[sample_data]), deltat = 1/52) ,stepwise=TRUE,approximation=FALSE, 
                                        allowdrift=FALSE,
                                        parallel = TRUE,  # speeds up computation
                                        trace=TRUE)) # trace my not be avaliable
              

              }
          # save each fitted ARIMA in fitted_models[[j]]
          fitted_models[[j]]<-fit
          # save the AIC of each fitted model
          model_aic_scores[model_id]<- fit$aic
          }
          # fit will be NULL if there is an error in the fitting process
          ,error = function(e){
          }
        )#end tryCatch
        # If fit == NULL, save the fitted model and the AIC as NANs  
        if(is.null(fit) || is.null(fitted_models[[j]])){ 
          fitted_models[[j]]<-NA
          model_aic_scores[model_id]<- NA
          #checker<-TRUE
        }
        # if auto==TRUE, use the auto.arima, which stop the loop and get only one value
        if(auto)
          break
        # sum 1 to model_id if using the ensembles
        model_id<-model_id+1 
        }
      
##############
# Forecasting 
##############
      
      # general initial variables
      predicted_value<- 0 # predicted values 
      pi<-numeric(my_n_ahead) # prediction intervals
      pi_lower<-numeric(23*my_n_ahead) # lower prediction interval
      pi_upper<-numeric(23*my_n_ahead) # upper prediction interval
      m<- numeric(my_n_ahead) # mean forecast value
      s<- numeric(my_n_ahead) # standard deviation
      sims<-c() # simulations for the mixture of gaussians  
      
      # Ensemble weights initial variables
      model_weights<- c()
      min_aic<- min(model_aic_scores, na.rm = TRUE) # min models' aic
      total_aic<-sum(exp(-.5*(model_aic_scores-min_aic)), na.rm =TRUE ) # sum of aics without nan values
      
      # Counts the number of models utilized in each forecast  
      my_n_models<-0 
      for(my_model in fitted_models){ # for each valid model in fitted_models sum 1 to my_n_models
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){
          my_n_models<-my_n_models+1 
        }}
      
      # Save the number of models utilized in each iteration  
      all_my_models<-append(all_my_models,my_n_models)
      
      # Generates a sequence of dates based on the last date by n week ahead
      flu_dates<- single_state[[1]]$target_end_date[sample_data] # current data inside the 104 weeks window
      last_date <- max(flu_dates) # last date of this window
      prediction_date <- seq.Date(from = last_date + 7 , by = "week", length.out = my_n_ahead) 
      
      #########################################
      # run the models saved on fitted models #
      for(my_model in fitted_models){
        ######################
        # if models are valid
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){ 
          # calculate the weights
          model_weights_<- exp(-.5*(my_model$aic - min_aic))/total_aic
          # weight the forecast
          predicted_value<- model_weights_*predict(my_model, n.ahead = my_n_ahead)$pred[my_n_ahead] + predicted_value 
          # generate prediction intervals
          pi<-forecast(my_model, h = my_n_ahead, level =  c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)) # get the levels of prediction intervals
          # weight the lower and upper bounds
          pi_lower<-model_weights_*pi[["lower"]]+pi_lower 
          pi_upper<-model_weights_*pi[["upper"]]+pi_upper
          # if predicted value is NA 
          if(is.na(predicted_value)){
            print("predicted_value is na")
          }
          #######################################################################################
          # simulate values for each model based on its weights to build the mixture of Gaussians
          new.sims<-c()
          fc <- forecast(my_model, h=my_n_ahead, level=99) ### forecast for 99% confidence
          m <- fc$mean[my_n_ahead]  ## mean forecast value 
          s <- ((fc$upper[my_n_ahead]-fc$lower[my_n_ahead])/2.58/2)  # standard deviation for for 99% confidence
          n<-ceiling(model_weights_*1e6) # number of simulations based on the weights
          new.sims <- rnorm(n, m=m, sd=s) # simulate values for each weighted model in as a gaussian
          sims <- c(sims, new.sims) # combine simulated values for each model
        }
      }
        # get the prediction dataset
        single_prediction<- data.frame("Pred_For_Date"= prediction_date[my_n_ahead], "Prediction" = predicted_value)# get the forecast for that date
        #rbind all forecasts and dates
        prediction_df<-rbind(prediction_df, single_prediction) 
        # Define the 23 quantiles
        probabilities <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99) 
        # get 23 predictive quantiles based on the mixture of gaussians 
        my_quantiles <- quantile(sims, probs=probabilities)
        
        #######################
        # Predictive quantiles
        
        # reset the initial prediction_df_quantile for each iteration
        prediction_df_quantile<- data.frame("pi_level"= NULL, "lower" = NULL, "uppper" = NULL, "quantile"= NULL, "point_forecast" = NULL)# empity dataframe of prediction intervals
        # save the 23 predictive quantiles, the upper and lower bounds, and the forecasts into a data frame
        for(j in 1:23){ 
          single_df_quantile<- data.frame("pi_level"= pi[["level"]][j]*.01, "lower" = pi_lower[[j+(23*(my_n_ahead-1))]], "uppper" = pi_upper[[j+(23*(my_n_ahead-1))]],"quantile"= my_quantiles[j], "point_forecast" = predicted_value)# fill prediction intervals dataframe with correct values for each week ahead
          prediction_df_quantile<-rbind(prediction_df_quantile, single_df_quantile) 
        }
        # put it into a list named with its forecast target_end_date
        prediction_quantile_ls[[toString(prediction_date[my_n_ahead])]]<-prediction_df_quantile 
    }
    
    # If don't have 104 weeks of data
    else
      print(paste0("Not enough values"))
  }
  
  # after everything is done
  print("Complete.")
  
  prediction[[1]]<-prediction_df # save the list of forecasts
  prediction_quantile[[1]]<- prediction_quantile_ls # save the list of predictive quantiles
  df_all_my_models<-data.frame(all_my_models) # save the number of models utilized in each forecast
  
  return(list("Point_ForeCast "=prediction, "Quantiles"=prediction_quantile, "Number_of_models"=df_all_my_models))
}

###############################################################
# This function combines the states data and the states codes # 
###############################################################

combining_states_data<-function(my_data=NULL, state_codes=NULL){
  
  my_data = subset(my_data, select = c(STATE,YEAR,EPI_WEEK,ILITOTAL))
  state_codes = subset(state_codes, select = c(location,location_name))
  names(state_codes)<- c('STATE_NUMBER','STATE')
  my_data<-cbind(my_data, MMWRweek2Date(MMWRyear=my_data$YEAR, MMWRweek=my_data$EPI_WEEK))
  suppressWarnings(invisible(my_data[apply(my_data, 1, purrr::compose(is.finite, all)),]))
  
  # Joining datasets
  final_data <- my_data %>%
    left_join(state_codes, by = "STATE")
  
  names(final_data)<- c('state_name','MMWRyear','MMWRweek','cases','target_end_date','location')
  final_data$location<-as.numeric(final_data$location)
  final_data$cases<-suppressWarnings(as.numeric(final_data$cases))
  
  final_data$target_end_date = as.Date(final_data$target_end_date,format = "%Y/%m/%d")
  final_data<-drop_na(final_data)
  grouped_data <-final_data %>% group_split(location)
  return(grouped_data)
}

#############################################################
# This function formats the dataset for calculating the WIS # correct
#############################################################

FormatForScoring_correct <- function(pred_intervals, data_, model_name, my_n_week_ahead=1, my_temporal_resolution="wk", my_target_variable="cases") {
  my_tibble<- NULL
  my_tibble<-tibble(model=c(""),forecast_date=c(as.Date(c()) ), location=c(double() ), horizon=c(double() ),
                    temporal_resolution=c(""), target_variable=c(""), target_end_date=c(as.Date(c()) ), type= c(""), quantile=c(double() ),
                    value =c(double()))
  
  for(i in 1:NROW(pred_intervals) ){
    dates_to_get<- names(pred_intervals[[i]])
    my_location<-data_$location[1]
    
    for(dates_ in dates_to_get){
      
      my_target_end_date<-as.Date(dates_)#-(7)
      my_forecast_date<-as.Date(dates_)-(7*my_n_week_ahead)
      
      my_tibble<- my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=my_n_week_ahead,
                                      temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "point", quantile=NA,
                                      value = expm1(pred_intervals[[1]][[dates_]]$point_forecast[1]) )
      for(quantile_level in pred_intervals[[i]][dates_]){
        
        my_quantile_value<-expm1(quantile_level$quantile)
        my_tibble<-my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=my_n_week_ahead,
                                       temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "quantile",
                                       quantile=quantile_level$pi_level, value = my_quantile_value)
      }
    }
  }
  return(my_tibble)
}

#############################################################
