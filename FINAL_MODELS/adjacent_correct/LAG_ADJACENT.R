#' ES_ADJACENT function
#' This function evaluates the performance of ensembles and AUTO ARIMAX models.
#' Our model calculates the absolute error and weighted interval scores for each U.S. state.
#' The ARIMAX models utilize mean cases by epidemiological weeks as exogenous variable.
#' This exogenous variable is calculated based on the 2 years before the prediction date. 
#' The function fits a rolling window of 2 years of data (104 weeks) to generate forecasts.
#' The user defines if it will use an AUTO ARIMA, or ensembles of 27 or 64 models and
#' choose the number of weeks ahead for each forecast.
#' 
#' @param current_state_tb A list containing the data for each U.S. state.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param ES27 A logical value indicating whether to use ensembles of 27 models. Default is \code{TRUE}.
#' @param ES64 A logical value indicating whether to use ensembles of 64 models. Default is \code{FALSE}.
#'
#' @return A list containing the forecast results and performance metrics.
#'
#'
#'
ES_ADJACENT<-function(current_state_tb, auto=FALSE, n_weeks_ahead=1, week_lag=0, ES27=TRUE, ES64=FALSE, window=104, list_of_states=list_of_states){
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
  
  results[[1]]%<-% ROLLING_ADJACENT(current_state_tb, n_weeks_ahead=n_weeks_ahead, week_lag=week_lag, window = window, order_params=my_order_params, auto=auto, list_of_states=list_of_states) %packages% "forecast" 
  # Resolve from system environment using the future package 
  # use new name  
  suppressWarnings(invisible(resolved(results[[1]]))) # 
  
  # Put forecasts, prediction intervals and number of models into separate lists. 
  list_all_pred<-list() # Empty list for forecasts
  list_all_pred_quantiles<- list() # Empty list for predictive quantiles
  list_number_of_models<- list() # Empty list for number of models
  # use new name
  # Get forecasts and dates
  list_all_pred[[1]]<- results[[1]][[1]][[1]]
  # use new name
  # Get prediction intervals
  list_all_pred_quantiles[[1]]<-results[[1]][[2]][[1]]
  # use new name
  # Get the number of models in each ensemble by date
  list_number_of_models[[1]]<-results[[1]][[3]][[1]]
  #######################################################################
  # Format a tibble of predictive quantiles and point forecasts for WIS #
  predictive_quantiles_tb<- FormatForWIS(list_all_pred_quantiles=list_all_pred_quantiles, current_state_tb, model_name = "TestModel", n_weeks_ahead = n_weeks_ahead) # Format for weighted interval score calculation
  ########################
  # Format truth for WIS #
  truth<-current_state_tb
  truth["target_variable"]<-"cases" # rename cases as target_variable
  truth["model"]<-predictive_quantiles_tb[1,"model"]# the model name from predictive quantiles in the truth
  truth<- truth %>% rename_at("cases", ~'value') # rename the column named cases to values
  ########################################################################################
  # Calculate the WIS using our prediction intervals and truth values by target_end_date #
  my_forecast_scores<-score_forecasts(predictive_quantiles_tb, truth) 
  #############################################
  # Get the number of models in the ensembles #
  all_models<-data.frame(as.Date(unique(predictive_quantiles_tb$target_end_date)),list_number_of_models[[1]]) #Get the number of models for each forecasted date
  colnames(all_models)<-c("Dates","Number_of_models")
  ##############################
  # Format flu cases and dates #
  all_cases<-data.frame(current_state_tb$cases, as.Date(current_state_tb$target_end_date)) # Get flu cases from current state
  colnames(all_cases)<-c("cases","Dates")
  ##############################
  # Format forecasts and dates #
  all_forecasts<-data.frame(expm1(list_all_pred[[1]]$Prediction), as.Date(list_all_pred[[1]]$predicted_date)) # Get ensembled forecast from current state
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
#' @param current_state_tb A list containing the data for each U.S. state.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param window An integer specifying the number of weeks to look back for the rolling window. Default is \code{104}.
#' @param order_params A list specifying the order parameters for the ARIMA model. Default is \code{NULL}.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#'
#' @return A list containing the forecast results and performance metrics.
#'
ROLLING_ADJACENT <- function(current_state_tb, week_lag=week_lag, n_weeks_ahead=1, window = 104, order_params=NULL, auto=FALSE, list_of_states=list_of_states) {
  
  #Creat a list with adjacent states indexes based on our list
  ##############################################################################  
  # MEAN CASES BY ADJACENT STATE
  us_adjacent_states <- list(
    Alabama = c( 8, 21, 39), # Florida, Georgia, Mississippi, Tennessee (48, 8, 21, 39)
    #correct
    Arizona = c(4, 5, 25, 28, 41), # California, Colorado, Nevada, New Mexico, Utah
    #correct
    Arkansas = c(15, 21, 22, 33, 39, 40), #  Louisiana, Mississippi, Missouri, Oklahoma, Tennessee, Texas
    #correct
    California = c(2, 25, 34), # Arizona, Nevada, Oregon
    #correct
    Colorado = c(2, 13, 24, 28, 33, 41, 47), #Arizona, Kansas, Nebraska, New Mexico, Oklahoma, Utah, Wyoming
    #correct
    Connecticut = c(29, 18, 36), # New York, Massachusetts, Rhode Island
    ## correct
    Delaware = c(17, 27, 35), # Maryland, New Jersey, Pennsylvania
    ## correct
    Georgia = c(1, 30, 37, 39), # Alabama, Florida, North Carolina, South Carolina, Tennessee(1, 48, 30, 37, 39)
    #correct
    Idaho = c(23, 25, 34, 41, 44, 47), # Montana, Nevada, Oregon, Utah, Washington, Wyoming
    # correct
    Illinois = c(11, 12, 14, 22, 46), # Indiana, Iowa, Kentucky, Missouri, Wisconsin
    # correct
    Indiana = c(10, 14, 19, 32), # Illinois, Kentucky, Michigan, and Ohio.
    # correct
    Iowa = c(10, 20, 22, 24, 38, 46), # Illinois, Minnesota, Missouri, Nebraska, South Dakota, Wisconsin
    # correct
    Kansas = c(5, 22, 24, 33), # Colorado, Missouri, Nebraska, Oklahoma.
    # correct
    Kentucky = c(10, 11, 22, 32, 39, 43, 45), # Illinois, Indiana, Missouri, Ohio, Tennessee, Virginia, West Virginia
    # correct
    Louisiana = c(3, 21, 40), # Arkansas, Mississippi, Texas
    # correct
    Maine = c(26), # New Hampshire
    # correct
    Maryland = c(7, 35, 43, 45), # Delaware, Pennsylvania, Virginia, West Virginia
    # correct
    Massachusetts = c(6, 26, 29, 36, 42), # Connecticut, New Hampshire, New York, Rhode Island, Vermont
    # correct
    Michigan = c(10, 11, 32, 46), # Illinois, Indiana, Ohio, and Wisconsin
    # correct
    Minnesota = c( 12, 31, 38, 46), # Iowa, North Dakota, South Dakota, and Wisconsin
    # correct
    Mississippi = c(1, 3, 15, 39), # Alabama, Arkansas, Louisiana, Tennessee
    # correct
    Missouri = c(12, 10, 14, 39, 3, 33, 13, 24), # Iowa, Illinois, Kentucky, Tennessee, Arkansas, Oklahoma, Kansas, Nebraska
    # correct
    Montana = c(9, 31, 38, 47), # Idaho, North Dakota, South Dakota, Wyoming
    # correct
    Nebraska = c(5, 12, 13, 22, 38, 47), # Colorado, Iowa, Kansas, Missouri, South Dakota, Wyoming
    # correct
    Nevada = c(2, 4, 9, 25, 41), # Arizona, California, Idaho, Oregon, Utah
    # correct
    `New Hampshire` = c(16, 18, 42), # Maine, Massachusetts, Vermont
    # correct
    `New Jersey` = c(7, 29, 35), # Delaware, New York, Pennsylvania
    # correct
    `New Mexico` = c(2, 5, 33, 40, 41), # Arizona, Colorado, Oklahoma, Texas, Utah
    # correct
    `New York` = c(6, 18, 27, 35, 36, 42), # Connecticut, Massachusetts, New Jersey, Pennsylvania, Rhode Island, Vermont
    # correct
    `North Carolina` = c(8, 30, 39, 43), # Georgia, South Carolina, Tennessee, Virginia
    # correct
    `North Dakota` = c(20, 23, 38),# Minnesota, Montana, South Dakota
    # correct
    Ohio = c(11, 14, 19, 35, 45),# Indiana, Kentucky, Michigan, Pennsylvania, West Virginia
    # correct
    Oklahoma = c(3, 5, 13, 22, 28, 40), # Arkansas, Colorado, Kansas, Missouri, New Mexico, Texas
    # correct
    Oregon = c(4, 9, 25, 44), # California, Idaho, Nevada, Washington
    # correct
    Pennsylvania = c(7, 17, 27, 29, 32, 45), # Delaware, Maryland, New Jersey, New York, Ohio, West Virginia
    # correct
    `Rhode Island` = c(6, 18), # Connecticut, Massachusetts
    # correct
    `South Carolina` = c(8, 30), # Georgia, North Carolina
    # correct
    `South Dakota` = c(12, 20, 23, 24, 31, 47), # Iowa, Minnesota, Montana, Nebraska, North Dakota, Wyoming
    # correct
    Tennessee = c(1, 3, 8, 14, 21, 22, 30, 43), # Alabama, Arkansas, Georgia, Kentucky, Mississippi, Missouri, North Carolina, Virginia
    # correct
    Texas = c(3, 15, 28, 33), # Arkansas, Louisiana, New Mexico, Oklahoma
    # correct
    Utah = c(2, 5, 9, 25, 28, 47), # Arizona, Colorado, Idaho, Nevada, New Mexico, Wyoming
    # correct
    Vermont = c(18, 26, 29), # Massachusetts, New Hampshire, New York
    # correct
    Virginia = c(14, 17, 30, 39, 45), # Kentucky, Maryland, North Carolina, Tennessee, West Virginia
    # correct
    Washington = c(9, 34), # Idaho, Oregon
    # correct
    `West Virginia` = c( 14, 17, 32, 35, 43), # Kentucky, Maryland, Ohio, Pennsylvania, Virginia
    # correct
    Wisconsin = c(10, 12, 19, 20), # Illinois, Iowa, Michigan, Minnesota
    # correct
    Wyoming = c(5, 9, 23, 24, 38, 41), # Colorado, Idaho, Montana, Nebraska, South Dakota, Utah
    # correct
    Florida = c(1, 8) # Alabama, Georgia
  )
  

  # Function for calculating the mean cases by adjacent states for the time series.
  # I works for a given selection of adjacent states based on the current state.
  ###################################################################################
  sum_cases<-numeric(NROW(current_state_tb))
  
  for (adjacent_state_name in us_adjacent_states[[unique(current_state_tb$state_name)]]){ # get the names of the adjacent states based on the current state name and our list of adjacent states.
    cases<- list_of_states[[adjacent_state_name]]$cases
    sum_cases<-cases + sum_cases
  }
  
  mean_cases_adjacent_states<-sum_cases/length(us_adjacent_states[[unique(current_state_tb$state_name)]]) # get the name of the state
  
  ###################################################################################
  # List with the data from a single state
  #  current_state_ls=list(current_state_tb) # !!!!!!! CHANGE NAME
  # All models in the ensemble
  N_of_models<-c() # !!!!!!! CHANGE NAME
  # Final predictions and dates list
  prediction<-list() # empty list for forecasts !!!!!!! CHANGE NAME
  # Predictions and dates data frame 
  prediction_df<- data.frame("predicted_date"= NULL, "Prediction" = NULL) # empty data frame for forecasts and dates
  # Final predictive quantiles list
  prediction_quantile<-list() # empity list for prediction intervals !!!!!!! CHANGE NAME
  # List in which each predictive quantiles data frame will be added with forecast target_end_date as its name
  prediction_quantile_ls<- list() 
  
#  for(iter in  1:(NROW(current_state_tb$cases)-(window))){ # rolling window of 104 weeks up to the end of the ILI cases dataset
    # Empty list for fitted ARIMAS
                # 2               :    105
  for(iter in  (1+week_lag):(NROW(current_state_tb$cases)-(window))){ # rolling window of 104 weeks up to the end of the ILI cases dataset
      
    current_state_rolling_window<- (iter):(window+iter-1)  
    #exog_rolling_window    2 - 1    :     104 +2 - 1 - 1 
    
    exog_rolling_window<- (iter-week_lag):(window+iter-1-week_lag)    
    
    
    #print(iter)
    fitted_models<-list()
    # Window of 104 weeks.
    # rolling_window<- iter:(window+iter-1) 
    # AIC scores for weighting the ensembles
    model_aic_scores<-c() 
    # Each model id, use later on a for loop
    model_id<-1
    
    # rolling window of exogenous time series based on the adjacent states cases
    ############################### build the ADJACENT dataframe #####################################
    adjacent_states_dataset= data.frame(mean_cases_adjacent_states[exog_rolling_window])
    
    # This is going to use the data from adjacent state for the day the forecast is being made.
    exog_var<-c(mean_cases_adjacent_states[104+(iter-1)], mean_cases_adjacent_states[104+(iter-1)], mean_cases_adjacent_states[104+(iter-1)], mean_cases_adjacent_states[104+(iter-1)])
#    exog_var<-c(adjacent_states_dataset[104,1], adjacent_states_dataset[104,1], adjacent_states_dataset[104,1], adjacent_states_dataset[104,1])
    ############################### build the ADJACENT dataframe #####################################
    
    ##########
    # Fitting 
    ##########
    
    # Start if ILI data has 104 elements
    if(length(current_state_tb$cases[current_state_rolling_window])==window){ 
      #    if(n_unique(log(current_state_ls[[1]]$cases[rolling_window]))==window){ 
      
      # run 1 time for the auto.arima or run 27 or 64 times for the ensembles
      for(j in 1:nrow(order_params)){
        fit<- NULL # start with fit as NULL
        
        # try to fit an ARIMA model
        tryCatch(
          expr = {
            if(!auto){
              # if auto = FALSE, run the ensemble of 27 or 64
              
              # fit ensembles of ARIMAs on log1p of the data
              fit<-Arima(log1p(current_state_tb$cases[current_state_rolling_window]), xreg=log1p(adjacent_states_dataset[,1]), order = order_params[j,], method = "CSS-ML") #
            }
            
            # if auto = TRUE, run auto.arima
            else{
              # fit an auto.arima on the log1p of the data
              fit<-invisible(auto.arima(log1p(current_state_tb$cases[current_state_rolling_window]), xreg=log1p(adjacent_states_dataset[,1]) ,stepwise=TRUE,approximation=FALSE,
                                        seasonal=FALSE, # This will extent to SARIMA
                                        allowdrift=FALSE, # speeds up computation
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
      pi<-numeric(n_weeks_ahead) # prediction intervals
      pi_lower<-numeric(23*n_weeks_ahead) # lower prediction intervals, it will have 23*N_week_ahead values
      pi_upper<-numeric(23*n_weeks_ahead) # upper prediction intervals, it will have 23*N_week_ahead values
      m<- numeric(n_weeks_ahead) # mean forecast value
      s<- numeric(n_weeks_ahead) # standard deviation
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
      
      ######################################################
      # Save the number of models utilized in each iteration  
      N_of_models<-append(N_of_models,my_n_models)
      
      ######################################################################
      # Generates a sequence of dates based on the last date by n week ahead
      weekly_dates<- current_state_tb$target_end_date[current_state_rolling_window] # current data inside the 104 weeks window
      last_date <- max(weekly_dates) # last date of this window
      my_predicted_dates <- seq.Date(from = last_date + 7 , by = "week", length.out = n_weeks_ahead) 
      predicted_date<-my_predicted_dates[n_weeks_ahead]
      
      ##########################################
      # run the models stored at fitted models #
      for(my_model in fitted_models){
        ######################
        # if models are valid
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){ 
          # calculate the weights based on the total AIC previously calculated
          model_weights_<- exp(-.5*(my_model$aic - min_aic))/total_aic
          
          #######################################################################################
          # simulate values for each model based on its weights to build the mixture of Gaussians
          new.sims<-c()
          # ADJACENT                   
          fc <- forecast(my_model, h=n_weeks_ahead, xreg=log1p(exog_var[1:n_weeks_ahead]), level=99)
          #fc <- forecast(my_model, h=n_weeks_ahead, level=99) ### forecast for 99% confidence
          m <- fc$mean[n_weeks_ahead]  ## mean forecast value 
          s <- ((fc$upper[n_weeks_ahead]-fc$lower[n_weeks_ahead])/2.58/2)  # standard deviation for for 99% confidence
          n<-ceiling(model_weights_*1e6) # number of simulations based on the weights
          new.sims <- rnorm(n, m=m, sd=s) # simulate values for each weighted model in as a gaussian
          sims <- c(sims, new.sims) # combine simulated values for each model
          
          #####################################
          # calculate the ensemble prediction #
          predicted_value <- model_weights_*m + predicted_value ### m = mean forecast for 99% confidence
          
          ###########################
          # if predicted value is NA 
          if(is.na(predicted_value)){
            print("predicted_value is na")
          }
          
          ##########################################################################################
          # weight lower and upper bounds and generate prediction interval levels for N weeks ahead
          #pi<-forecast(my_model, h = n_weeks_ahead, level =  c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)) # get the levels of prediction intervals
          
          # ADJACENT
          pi<-forecast(my_model, h = n_weeks_ahead, xreg=log1p(exog_var[1:n_weeks_ahead]), level =  c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99))
          
          # weight the lower and upper bounds
          pi_lower<-model_weights_*pi[["lower"]]+pi_lower 
          pi_upper<-model_weights_*pi[["upper"]]+pi_upper
        }
      }
      
      # get the prediction dataset
      single_prediction<- data.frame("predicted_date"= predicted_date, "Prediction" = predicted_value)# get the forecast for that date
      #rbind all forecasts and dates
      prediction_df<-rbind(prediction_df, single_prediction) 
      # Define the 23 quantiles
      probabilities <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99) 
      # get 23 predictive quantiles based on the mixture of gaussians 
      my_quantiles <- quantile(sims, probs=probabilities)
      
      ########################
      # Predictive quantiles #
      ########################
      
      # reset the initial prediction_df_quantile for each iteration
      prediction_df_quantile<- data.frame("pi_level"= NULL, "lower" = NULL, "uppper" = NULL, "quantile"= NULL, "point_forecast" = NULL)# empity dataframe of prediction intervals
      
      # save the 23 predictive quantiles, the upper and lower bounds for the given week ahead, and the ensemble forecasts into a data frame
      for(j in 1:23){ 
        single_df_quantile<- data.frame("pi_level"= pi[["level"]][j]*.01, "lower" = pi_lower[[j+(23*(n_weeks_ahead-1))]], "uppper" = pi_upper[[j+(23*(n_weeks_ahead-1))]],"quantile"= my_quantiles[j], "point_forecast" = predicted_value)# fill prediction intervals dataframe with correct values for each week ahead
        prediction_df_quantile<-rbind(prediction_df_quantile, single_df_quantile) 
      }
      # put it into a list named with its forecast target_end_date
      prediction_quantile_ls[[toString(predicted_date)]]<-prediction_df_quantile 
    }
    
    # If don't have 104 weeks of data
    else
      print(paste0("Not enough values"))
  }
  
  # after everything is done
  print("Complete.")
  
  prediction[[1]]<-prediction_df # save the list of forecasts
  prediction_quantile[[1]]<- prediction_quantile_ls # save the list of predictive quantiles
  df_N_of_models<-data.frame(N_of_models) # save the number of models utilized in each forecast
  
  return(list("Point_ForeCast "=prediction, "Quantiles"=prediction_quantile, "Number_of_models"=df_N_of_models))
}

##########################################################
# This function splits the ILI data into separate states # 
##########################################################

combining_states_data<-function(ILI_data=NULL, state_codes=NULL){
  
  # select only the STATE, YEAR, EPI_WEEK, ILITOTAL columns from ILI data
  ILI_data = subset(ILI_data, select = c(STATE,YEAR,EPI_WEEK,ILITOTAL))
  # add a column with weekly dates
  ILI_data<-cbind(ILI_data, MMWRweek2Date(MMWRyear=ILI_data$YEAR, MMWRweek=ILI_data$EPI_WEEK))
  # keep only finite values
  #ILI_data[apply(ILI_data, 1, purrr::compose(is.finite, all)),]
  
  # select only location and location_name from state_codes
  state_codes = subset(state_codes, select = c(location,location_name))
  names(state_codes)<- c('STATE_NUMBER','STATE')
  
  # Joining datasets
  combined_data <- ILI_data %>%
    left_join(state_codes, by = "STATE")
  
  # Renaming, organizing variables types and removing NANs
  names(combined_data)<- c('state_name','MMWRyear','MMWRweek','cases','target_end_date','location')
  combined_data$location<-as.numeric(combined_data$location)
  combined_data$cases<-suppressWarnings(as.numeric(combined_data$cases))
  combined_data$target_end_date = as.Date(combined_data$target_end_date,format = "%Y/%m/%d")
  combined_data<-drop_na(combined_data)
  
  # split into different states
  states_data <-combined_data %>% group_split(location)
  return(states_data)
}

#############################################################
# This function formats the dataset for calculating the WIS # GOOD FUNCTION
#############################################################

FormatForWIS <- function(list_all_pred_quantiles, current_state_tb, model_name, n_weeks_ahead=1, my_temporal_resolution="wk", my_target_variable="cases") {
  my_tibble<- NULL
  
  # Create an empty tibble in the correct format
  my_tibble<-tibble(model=c(""),forecast_date=c(as.Date(c())), location=c(double()), horizon=c(double() ),
                    temporal_resolution=c(""), target_variable=c(""), target_end_date=c(as.Date(c())), type= c(""), quantile=c(double() ),
                    value =c(double()))
  
  # Loop over all the quantiles
  for(single_quantile in 1:NROW(list_all_pred_quantiles) ){
    # Get the dates of each single quantiles 
    predicted_dates_ls<- names(list_all_pred_quantiles[[single_quantile]])
    # Get the location number 
    my_location<-current_state_tb$location[1]
    
    # Get predicted dates from the list of predicted dates
    for(predicted_date_ in predicted_dates_ls){
      # Save predicted date as target_end_date
      my_target_end_date<-as.Date(predicted_date_)
      # Save the date in which the prediction was made as (predicted_date - (n_weeks_ahead*7))
      my_forecast_date<-as.Date(predicted_date_)-(7*n_weeks_ahead)
      
      # Add rows with predictions to the tibble as point_forecast
      my_tibble<- my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=n_weeks_ahead,
                                      temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "point", quantile=NA,
                                      value = expm1(list_all_pred_quantiles[[1]][[predicted_date_]]$point_forecast[1])) # exponentiating the predictions back
      
      # Add rows with predictive_quantiles to the tibble as quantiles      
      for(quantile_level in list_all_pred_quantiles[[single_quantile]][predicted_date_]){
        my_quantile_value<-expm1(quantile_level$quantile) # exponentiating the predictions back
        my_tibble<-my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=n_weeks_ahead,
                                       temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "quantile",
                                       quantile=quantile_level$pi_level, value = my_quantile_value)
      }
    }
  }
  return(my_tibble)
}
