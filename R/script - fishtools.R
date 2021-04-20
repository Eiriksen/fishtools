

#' Function for getting the mean temperature of a time period
#' Dates must be supplied as a date object
tempData_mean = function(tb_temp, date_start, date_end, column_temp) {
  if(is.na(date_start) | is.na(date_end) | is.na(column_temp)) return(NA)
  return(
   mean(unlist(tb_temp[column_temp])[tb_temp$date < date_end & tb_temp$date > date_start],na.rm=T)
  )
}



#' calculates thermal growth coefficient (from Jobling: The thermal growth coefficient (TGC) model of fish growth: a cautionary note)
calc_TGC = function(W0,W1,temperature,time) 
{
  time <- as.numeric(time)
  TGC <- (((W1^(1/3))-(W0^(1/3)))/(temperature*time)*1000)
  return(TGC)
}

#' predicts growth based on TGC, same model as above
calc_growth_TGC <- function(W0, TGC, temperature, time)
{
  time <- as.numeric(time)
  growth <- as.numeric((W0^(1/3)+((TGC/1000)*(temperature*time)) )^3)
  return(growth)
}


#' For a whole dataframe, calculates Thermal Growth Coefficient TGC for a given period
calc_TGC_df <- function(tb_fish,tb_temp,period)
{
  time <- tb_fish[[glue("date.{period}")]] - tb_fish[[glue("date.{period-1}")]]
  W0   <- tb_fish[[glue("weight.{period-1}")]]
  W1   <- tb_fish[[glue("weight.{period}")]]
  
  temperature <- mapply(
    FUN       = tempData_mean,
    date_start = ymd(tb_fish[[glue("date.{period-1}")]]),
    date_end   = ymd(tb_fish[[glue("date.{period}")]]),
    column    = "temp",
    MoreArgs  = list(tb_temp=tb_temp)
  )

  TGC <- calc_TGC(W0, W1, temperature, time)
  
  return(TGC)
}


#' @title predict_weights
#' @description applies to a whole table of weight data for fish, combines growth data with temperature data to predict future weights. Uses Jobling's TGC formula. 
#' @param tb_fish The table of fish growth data to predict weight for. Table must be in "wide" format; that is, one row pr individual, and columns weight.1, weight.2, weight.3 ... date.1, date.2, date.3... etc for measured weights and their corresponding dates for different periods. Weights must be numerical. Dates must be in date format. The table can contain other columns, these will be ignored. Periods end when the fish has been weighed, so period 2, for example, are the dates between date.1 and date.2. Weight.2 is thus the weight at the end of period 2.
#' @param tb_temp A table that contains temperature data for the fish in tb_fish. Needs two columns: date and temp. Column date must be in date format. Is used to calculate mean temperatures between periods.
#' @param dates A vector of dates to predict new weights for. Must be contained in a vector and formatted as dates (see example)
#' @param TGC_baseline If fish have only been measured once (there's just one period in the table), a thermal growth coefficient can't be calculated, so a TGC needs to be assumed (check data for your species). 
#' @param period_start optional, if for example set to 3, then periods 1 and 2 will be ignored when estimating thermal growth coefficients. 
#' @param tgc_p_modifier optional, allows the thermal growth coefficients to be modified by. If tgc_p_modifier is for example set to 1.2, then all thermal growth coefficients are multiplied by 1.2
#' @param temp_abs_modifier optional, allows future temperatures to be modified. Setting this to 2 will increase predicted temperature by 2.
#' @examples tb_fish_predicted <- predict_weight(tb_fish,tb_temp,c(ymd(20210601),ymd(20211001),ymd(20211201)))
#' @examples tb_fish_predicted <- predict_weight(tb_fish,tb_temp,today()+weeks(1:5))
#' @export
predict_weights <- function(tb_fish, tb_temp, dates, TGC_baseline=NULL, period_start=NULL, tgc_p_modifier=1, temp_abs_modifier=0)
{
  require(glue)
  require(ggplot2)
  
  # Obtain all the periods so far as a vector (of numbers) from tb_fish
  # remove periods that are after the specified start period (in case user want's to exclude some early data)
  periods <- extract_periods(tb_fish)
  if(!missing(period_start)) periods <- periods[periods <= period_start]
  
  if (length(periods)==1)
  {
    if (missing(TGC_baseline)) stop("Only one period in dataset. Baseline TGC needs to be supplied.")
    else tgc_weighed = TGC_baseline * tgc_p_modifier
  }
  else
  {
    message("Calculating Thermal Growth Coefficients (TGC) for all fish for all periods...")
    # calculate TGC for all individuals for all periods (but not the first one)
    tgc <- as.data.frame(sapply(
      X = periods[2:length(periods)], 
      FUN = function(i){calc_TGC_df(tb_fish,tb_temp,i)}
    ))
    message("TGC calculated for periods: ",paste(periods[2:length(periods)],collapse=", "))
    message("See histogram for TGC values:")

    # Plot all TGC for the user
    plot_tgc <- tgc %>% 
      pivot_longer(everything(),names_to="period",values_to="tgc") %>%
      ggplot() +
      aes(x=tgc) +
      xlim(0,5) +
      geom_histogram(bins=50) +
      facet_wrap(~period)
    print(plot_tgc)
    
    # for each fish create one new tgc that is the mean of these TGC values, 
    # but weighing the last one by 2
    tgc_weighted <- apply(tgc,1,FUN=function(x){
      x=na.omit(x)
      if(length(x)==0) return(NA)
      weighted.mean(x,c(rep(1,length(x)-1),2))
    })
  
    # apply tgc modifier 
    tgc_weighted = tgc_weighted*tgc_p_modifier
  }
  
  # the first period to predict
  period_cur = max(periods+1)

  # for each date, record the predicted weight
  for (date in as.list(dates)){
    
    # get the date from the previous measurement 
    date_prev   <- glue("date.{period_cur-1}")
    # calculate the mean temp from the previous period until the new date
    # and add the temp modifier
    temperature = mapply(
      tempData_mean,
      date_start = tb_fish[[date_prev]],
      date_end   = date,
      column_temp = "temp",
      MoreArgs  = list(tb_temp=tb_temp)
      ) + temp_abs_modifier

   
    # find the last recorded weight of the fish
    W0 <- tb_fish[[glue("weight.{period_cur-1}")]]
    
    # create new column for the weight and date of the fish
    weight_new <- glue("weight.{period_cur}")
    date_new   <- glue("date.{period_cur}")
    
    # predict the growth for the newdate
    tb_fish[[date_new]]   <- date
    tb_fish[[weight_new]] <- calc_growth_TGC(
      W0 = W0,
      TGC  = tgc_weighted,
      time = date-tb_fish[[date_prev]],
      temperature = temperature
    )

    message(glue("{weight_new} predicted at {date}"))

    period_cur = period_cur+1
  }

  return(tb_fish)
}

#' Takes a wide-format dataset of individuals, with period-associated measurements denoted by dot-number format, e.g:
#' Weight.1, weight.2, weight.3 or length.1 length.2 or date.1, date2. etc
#' Returns a numeric vector of all the periods in that dataset eg: c(1,2,3,4,5) for a dataset with 5 periods
#' @export
extract_periods <- function(tb_fish){
  require(tidyverse)

  periods <- colnames(tb_fish) %>%
    str_extract("weight[.][0-9]*") %>%
    sub("weight.","",.) %>%
    str_extract("\\-*\\d+\\.*\\d*") %>%
    as.numeric() %>%
    na.omit() %>%
    unique()

  return(periods)
}


## Feed functions -------------------------------------------------------------
# Based on the input of a list of weights (the weights of all fishes in a group)
# Returns table with the different feed sizes and their respective distributions in the feed
calc_feedSizeDist = function(weights)
{
  tabl = case_when(weights <= 15 ~ 1.2,
                   weights > 15  & weights <= 30  ~ 1.7,
                   weights > 30  & weights <= 70  ~ 2.5,
                   weights > 70  & weights <= 150 ~ 3.5,
                   weights > 150 & weights <= 425 ~ 5.0,
                   weights > 425 & weights <= 2300~ 7.0,
                   weights > 2300 & weights <= 9000~9.0,
                   weights > 9000 ~ 100,
                   TRUE ~ NaN)
  return(tabl)
}


# uses the same function as above, but returns results as percentage for each feed size
calc_feedSizePercentage = function(weights) {
  weightClasses = calc_feedSizeDist(weights)
  dat = data.frame(weightClasses)
  totalLength = length(dat$weightClasses[!is.na(dat$weightClasses)])
  
  return(
    dat %>% group_by(weightClasses) %>% summarise( perc = length(weightClasses)/totalLength*100 )
  )
}


calc_feedSizePercentage_special = function(weights,prefix) {
  weightClasses = paste(prefix,calc_feedSizeDist(weights),sep="")
  dat = data.frame(weightClasses)
  totalLength = length(dat$weightClasses[!is.na(dat$weightClasses)])
  
  return(
    dat %>% group_by(weightClasses) %>% summarise( perc = length(weightClasses)/totalLength*100 )
  )
}


# calculates the total amount of feed given between two dates
# internal function
feedData_amount_fed = function(dates, feed, date1, date2) {
  return(
    sum( feed[ dates > date1 & dates < date2], na.rm=T )
  )
}


# calculates the amount of feed given a specific tank between two dates
feedData_feedGiven = function(df_feed,tankID, date1, date2) {
  data_feed_tank = df_feed[which(df_feed$tank == tankID),]
  amount = feedData_amount_fed(data_feed_tank$date, data_feed_tank$feed, date1, date2)
  return(amount)
}




# function for converting the feed datasat to a "pr-day" type dataset
feedData_convert_feedPrDate = function(df_feed,tank,date) {
  #BE AWARE: this function does not work well for the very last period of feeding
  
  #find the last date where this tank was "fed"
  prevDates = df_feed$date[df_feed$date <= date & df_feed[tank] != 0]
  prevDate  = max(prevDates,na.rm=T)
  #get the next date where this tank is being "fed"
  comingDates = df_feed$date[df_feed$date > date]
  if (length(na.omit(comingDates))==0) return(NA)
  nextDate = min(comingDates,na.rm=T)
  #get the amount of feed at the intial date
  amount = as.numeric(unlist(df_feed[tank][df_feed$date == prevDate,1]))
  
  #get the number of days between last and next date
  numDays = as.numeric(nextDate-prevDate)
  
  #divide amount fed by number of dates
  dailyAmount = amount/numDays
  return(dailyAmount[1])
}




calc_feedNeeds = function(tb_fish, period_start, period_end, round_to=0){
  # Create binnies to dump food into
  binnies = data.frame("period"=0, "size_1p2"=0, "size_1p7"=0, "size_2p5"=0, "size_3p5"=0, "size_5"=0, "size_7"=0, "size_9"=0)
  # Fill bins function
  fillBins = function(bins, bin, amount, identifier){
    # (Neccessary function for compiling feed data to binnies
    #)
    # Function that uses dataframe with columns (bins) corresponding to each feed size pr week
    # - this function then reads data about feed consumption for a list of fishsh (feed)
    # - as well as data about which size each fish uses (bin)
    # - and creates a row for that week with the appropriate amonut of food in each bin, t
    # - and then adds that row to the dataframe that was supplied to it
    size_1p2 = sum(amount[bin==1.2], na.rm=T)
    size_1p7 = sum(amount[bin==1.7], na.rm=T)
    size_2p5 = sum(amount[bin==2.5], na.rm=T)
    size_3p5 = sum(amount[bin==3.5], na.rm=T)
    size_5   = sum(amount[bin==5], na.rm=T)
    size_7   = sum(amount[bin==7], na.rm=T)
    size_9   = sum(amount[bin==9], na.rm=T)
    
    bins = rbind(bins, c(identifier, size_1p2, size_1p7, size_2p5, size_3p5, size_5, size_7, size_9))
    return(bins)
  }
  
  for (i in period_start:period_end){
    #for each period:
    #figure out how much the fishes grew since last period
    growth = tb_fish[paste("weight.",i,sep="")]-tb_fish[paste("weight.",i-1,sep="")]
    #Using their FCR, figure out how much feed they needed for that period
    feed_amount = growth/tb_fish$FCR
    #Using their size, find the size of feed they needed for that period
    feed_size <- calc_feedSizeDist(tb_fish[paste("weight.",i-1,sep="")] )
    #add the feed for all the fish for that week to binnies
    binnies <- fillBins(binnies, unlist(feed_size), unlist(feed_amount), i)
  }
  
  # Return the resulting dataset for fishes, feed, and the summary of feed needs.
  
  return(list(
    "feed"      = binnies,
    "summary"   = list(
      "Size_1p2"=(sum(binnies$size_1p2)/1000) %>% round(round_to),
      "Size_1p7"=(sum(binnies$size_1p7)/1000) %>% round(round_to),
      "Size_2p5"=(sum(binnies$size_2p5)/1000) %>% round(round_to),
      "Size_3p5"=(sum(binnies$size_3p5)/1000) %>% round(round_to),
      "Size_5"  =(sum(binnies$size_5)/1000) %>% round(round_to),
      "Size_7"  =(sum(binnies$size_7)/1000) %>% round(round_to),
      "Size_9"  =(sum(binnies$size_9)/1000) %>% round(round_to))
  ))
}




# write_size and feed table
write_size_and_feed_table = function(tb_fish,tb_temp,period_start,filename="table - predicted biomass and feed size.csv")
{
  df_results <- data.frame()
  periods <- extract_periods(tb_fish)
  
  for (period in periods[period_start:length(periods)])
  {
    
    tb_fish_c = tb_fish %>% filter(temp=="cold")
    tb_fish_w = tb_fish %>% filter(temp=="hot")
    proportions_feed_w <- tb_fish_w[[glue("weight.{period}")]] %>% na.omit() %>% calc_feedSizePercentage_special("W.") %>% funkyTranspose()
    proportions_feed_c <- tb_fish_c[[glue("weight.{period}")]] %>% na.omit() %>% calc_feedSizePercentage_special("C.") %>% funkyTranspose()
    
    tank_biomass     <- data.frame(tank=sort(unique(tb_fish$tank.2))) %>% na.omit()
    tank_biomass$biomass <- apply(tank_biomass,1,FUN=function(x){
      tank_name <- x[["tank"]]
      tank_this <- tb_fish[tb_fish$tank.2==tank_name,]
      mean_mass <- tank_this[[glue("weight.{period}")]] %>% na.omit() %>% mean(na.rm=T)
      return(mean_mass)
    })
    tank_biomass <- tank_biomass %>% funkyTranspose()
    
    week             <- data.frame(week=week(min(tb_fish[[glue("date.{period}")]])))
    
    row              <- cbind(week,proportions_feed_c,proportions_feed_w,tank_biomass)
    df_results <- df_results %>% bind_rows(row) %>% round_df(1)
    
  }
  filename=glue(filename)
  save.data(df_results, filename)
  message("Table saved to ",filename)
  return(df_results)
  
}


