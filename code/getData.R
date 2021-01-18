library(openxlsx)
library(stringr)
library(readr)

getData <- function(country){
  if(country=='Switzerland'){
    return(BAGdata())
  }
  # else if(country=="BelgiumHosp"){
  #   return(BELHospData())
  # }
  # else if(country=="Belgium"){
  #   return(BELCaseData())
  #}
  else{
    worldData <- read.csv('worldData.csv')
    C <- worldData[which(worldData$Country==country),"New_cases"]
    start <- min(which(C >= 10))
    C <- C[start:length(C)]
    date <- worldData[which(worldData$Country==country),1]
    date <- date[start:length(date)]
    N <- length(C)
    return(list(C = C,N = N, date=date))
  }
}

BAGdata <- function(){
  BAGdata         <- read.xlsx(sep=",",startRow = 8, detectDates = TRUE,
                               "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-datengrundlage-lagebericht.xlsx.download.xlsx/200325_Datengrundlage_Grafiken_COVID-19-Bericht.xlsx")
  
  names(BAGdata) <- c("date","cases","casesCumul","hospitalized","hospitalizedCumul",
                      "deaths","deathsCumul")
  C <- round(BAGdata$cases[-length(BAGdata$cases)]) #remove last day (0 or incomplete data)
  N <- length(C)
  date <-BAGdata$date[-length(BAGdata$date)]
  return(list(C = C,N = N, date=date))
}

BELHospData <- function(){
  BELHosp <- read.table("https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv",
                        header=TRUE,sep=",")
  BELHosp$PROVINCE[BELHosp$PROVINCE == "LiÃ¨ge"] = "Liege"
  
  hosp <- data.frame(as.Date(unique(BELHosp$DATE)))
  names(hosp) <- c("date")
  for (province in unique(BELHosp$PROVINCE)){
    hosp <- cbind(hosp,BELHosp$NEW_IN[BELHosp$PROVINCE == province])
    names(hosp)[length(names(hosp))] <- province
  }
  C <- rowSums(hosp[-1])
  N <- length(C)
  return(list(C = C,N = N, date=hosp$date))
}

BELCaseData <- function(){
  BELdata <- read.table("https://epistat.sciensano.be/Data/COVID19BE_tests.csv",
                        header=TRUE,sep=",")
  data <- data.frame(as.Date(unique(BELdata$DATE)))
  names(data) <- c("date")
  for (province in unique(BELdata$PROVINCE)){
    if(!(is.na(province))){
      data <- cbind(data,BELdata$TESTS_ALL_POS[BELdata$PROVINCE == province])
      names(data)[length(names(data))] <- province
    }
    else{
      data <- cbind(data,BELdata$TESTS_ALL_POS[is.na(BELdata$PROVINCE)])
      names(data)[length(names(data))] <- "NA"
    }
  }
  C <- rowSums(data[-1],na.rm=TRUE)
  N <- length(C)
  return(list(C = C[1:(N-3)],N = (N-3), date=data$date[1:(N-3)]))
}
