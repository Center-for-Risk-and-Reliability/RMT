# Load the Shiny library needed to run the app
library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("Life-Stress Modeling Toolkit"),

  sidebarLayout(
    sidebarPanel(
      helpText(
        "The life-stress modeling toolkit takes data and fits it to a
        specified life-stress model with respect to a given life distribution.
        To start, prepare your data in a CSV file with the data
        in the first column, the censored status in the second column (right
        -censored data = ''0'' and non-censored = ''1'', and any stress values
        in the third column and onward.  The censored column has to be included."
      ),
      # Input panel for Data
      fileInput(
        "datafile", "Choose CSV File", multiple = TRUE,
        accept = c("text/csv","text/comma-separated-values,text/plain",
                   ".csv")
      ),
      # Input panel for Life Label
      textInput(
        "xlabel", h4("X Label"),
        value = "X"
      ),
      # Input panel for Stress(es) Label
      textInput(
        "slabel", h4("Stress(es) Label"),
        value = "Y"
      ),
      helpText(
        "Enter the units for life and stress here.  For multi-stress cases, separate each stress by
        a comma.  Example: Enter ''K,RH%'' for dual stresses of temperature and relative humidity."),
      # Input panel for the plotting position options
      selectInput(
        "plotposit", h4("Plotting Position"),
        choices = list(
          "Blom" = 1, "Mean" = 2, "Median" = 3, "Midpoint" = 4,
          "Jenkinson (Beard)" = 5, "Benard and Bos-Levenbach" = 6,
          "Tukey" = 7, "Gringorten" = 8, "Kaplan-Meier" = 9,
          "Nelson-Aalen" = 10),
        selected = 1
      ),
      # Input panel for Probability Plot/Life Distribution
      selectInput(
        "pp", h4("Life Distribution"),
        choices = list("Weibull" = 1,  "3P Weibull" = 2, "Lognormal" = 3,
                       "Normal" = 4, "Exponential" = 5,
                       "2P Exponential" = 6),
        selected = 1
      ),
      # Life-Stress Model
      selectInput(
        "ls", h4("Life-Stress Model"),
        choices = list("Linear" = 1, "Exponential" = 2,  "Arrhenius" = 3, "Eyring" = 4,
                       "Eyring Version 2" = 5, "Power" = 6, "Inverse Power" = 7,
                       "Logarithmic" = 8, "Temperature-Humidity (2-Stress)" = 9,
                       "Temperature-Nonthermal (2-Stress)" = 10, "Eyring Version 3 (2-Stress)" = 11,
                       "Multi-Stress (2 (or greater)-Stress" = 12),
        selected = 1
      ),
      # Input panel for confidence
      numericInput(
        "conf", h4("Confidence Bounds (%)"), min = 1, max = 99.99, value = 95
      ),
      # Input panel for Use Stress
      textInput(
        "usestress", h4("Use Stress"), value = "5"
      ),
      helpText(
        "Insert use stress here.  If use-stress has multiple variables, enter
        the stresses as ''stress1, stress2, stress3, ...''."),
      textInput(
        "stressmin", h4("Minimum Stress"), value = "1"
      ),
      textInput(
        "stressmax", h4("Maximum Stress"), value = "10"
      ),
      helpText(
        "Insert stress range (minimum and maximum) here.  If stress is made up of multiple
        variables, then enter stresses as ''Stress1, Stress2, Stress3, ...''.")
    ),
    mainPanel(
      # Output results in tabular form:
      # Tab 1) Unreliability or reliability plot
      # Tab 2) Table data: x, F, R, F confidence, and R confidence
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Input data",
          tableOutput("rawdata")
        ),
        tabPanel(
          "Probabaility Plot",
          plotOutput("plot1"),
          verbatimTextOutput("paramssum1"),
          verbatimTextOutput("paramsLSQ")
          #downloadButton("downloadprobplot", "Download")
        ),
        tabPanel(
          "MLE Confidence Intervals",
          tableOutput("paramsMLE")
          #downloadButton("downloadprobplot2", "Download")
        ),
        tabPanel(
          "Acceleration Factors",
          tableOutput("AFtable"),
          verbatimTextOutput("paramssum2")
          #downloadButton("downloadprobplot2", "Download")
        ),
        tabPanel(
          "Relationship Plot",
          plotOutput("plot2"),
          verbatimTextOutput("uselife")
          #downloadButton("downloadprobplot2", "Download")
        )
      )
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  # Read the data table
  rawdat <- reactive({
    # Read the data from the CSV
    req(input$datafile)
    data1 <- read.csv(input$datafile$datapath)
  })

  # Identify the Plotting Position that goes into the label
  plotpositlabel <- reactive({
    if (input$plotposit == 1){
      plotlab1 = "Blom"
    }
    if (input$plotposit == 2){
      plotlab1 = "Mean"
    }
    if (input$plotposit == 3){
      plotlab1 = "Median"
    }
    if (input$plotposit == 4){
      plotlab1 = "Midpoint"
    }
    if (input$plotposit == 5){
      plotlab1 = "Beard"
    }
    if (input$plotposit == 6){
      plotlab1 = "BernardBosLevenbach"
    }
    if (input$plotposit == 7){
      plotlab1 = "Tukey"
    }
    if (input$plotposit == 8){
      plotlab1 = "Gringorten"
    }
    if (input$plotposit == 9){
      plotlab1 = "KaplanMeier"
    }
    if (input$plotposit == 10){
      plotlab1 = "NelsonAalen"
    }
    return(plotlab1)
  })

  # Process the plot and LSQ parameters based on Probability Plot type
  plotparamout1 <- reactive({
    if (input$pp == 1){
      plotandparams<-probplot.wbl(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp == 2){
      plotandparams<-probplot.wbl3P(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp == 3){
      plotandparams<-probplot.logn(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp == 4){
      plotandparams<-probplot.nor(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp == 5){
      plotandparams<-probplot.exp(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp == 6){
      plotandparams<-probplot.exp2P(rawdat(),plotpositlabel(),input$xlabel)
    }
    return(plotandparams)
  })

  multistresslist <- reactive({
    databystress<-checkstress(rawdat())
    if (is.null(dim(databystress))){
      stresslist<-colnames(rawdat(), do.NULL = TRUE, prefix = "col")
      stresslist<-stresslist[3:length(stresslist)]
      return(stresslist)
    }
  })

  outputlist1 <- reactive({
    databystress<-checkstress(rawdat())
    if (!is.null(dim(databystress))){
      # Single Stress (this only works for multi-stress so include an error message)
      if (input$pp == 1){
        cat("alpha = ",plotparamout1()[[1]][[1]],", beta = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp == 2){
        cat("alpha = ",plotparamout1()[[1]][[1]],", beta = ",plotparamout1()[[1]][[2]],", gamma = ",plotparamout1()[[1]][[3]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp == 3){
        cat("log-mu = ",plotparamout1()[[1]][[1]],", log-sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp == 4){
        cat("mu = ",plotparamout1()[[1]][[1]],", sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp == 5){
        cat("lambda = ",plotparamout1()[[1]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp == 6){
        cat("theta = ",plotparamout1()[[1]][[1]],", sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
    } else {
      # Multi Stress
      for(i in 1:length(databystress)){
        stressrow<-rep(1,2*length(multistresslist()))

        for(j in 1:length(plotparamout1()[[i*3-2]])){
          stressrow[j*2-1]<-multistresslist()[j]
          stressrow[j*2]<-plotparamout1()[[i*3-2]][[j]]
        }
        if (input$pp == 1){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout1()[[i*3-1]][1],", beta = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp == 2){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout1()[[i*3-1]][1],", beta = ",plotparamout1()[[i*3-1]][2],", gamma = ",plotparamout1()[[i*3-1]][3],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp == 3){
          cat("Stress ",i,stressrow,":log-mu = ",plotparamout1()[[i*3-1]][1],", log-sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp == 4){
          cat("Stress ",i,stressrow,":mu = ",plotparamout1()[[i*3-1]][1],", sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp == 5){
          cat("Stress ",i,stressrow,":lambda = ",plotparamout1()[[i*3-1]],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp == 6){
          cat("Stress ",i,stressrow,":theta = ",plotparamout1()[[i*3-1]][1],", sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
      }
    }
  })

  # Labeling the life distribution and the life-stress models and their columns
  distlslab <- reactive({
    # Life Distribution label
    if (input$pp == 1){
      distlab <- "Weibull"
      distparam <- c("beta")
    }
    if (input$pp == 2){
      distlab <- "3PWeibull"
      distparam <- c("beta")
    }
    if (input$pp == 3){
      distlab <- "Lognormal"
      distparam <- c("sigma")
    }
    if (input$pp == 4){
      distlab <- "Normal"
      distparam <- c("sigma")
    }
    if (input$pp == 5){
      distlab <- "Exponential"
      distparam <- NULL
    }
    if (input$pp == 6){
      distlab <- "2PExponential"
      distparam <- c("sigma")
    }

    # Life-Stress Label
    if (input$ls == 1){
      lslab <- "Linear"
      lsparam <- c("a","b")
    }
    if (input$ls == 2){
      lslab <- "Exponential"
      lsparam <- c("a","b")
    }
    if (input$ls == 3){
      lslab <- "Arrhenius"
      lsparam <- c("Ea","b")
    }
    if (input$ls == 4){
      lslab <- "Eyring"
      lsparam <- c("a","b")
    }
    if (input$ls == 5){
      lslab <- "Eyring2"
      lsparam <- c("a","b")
    }
    if (input$ls == 6){
      lslab <- "Power"
      lsparam <- c("a","b")
    }
    if (input$ls == 7){
      lslab <- "InversePower"
      lsparam <- c("a","b")
    }
    if (input$ls == 8){
      lslab <- "Logarithmic"
      lsparam <- c("a","b")
    }
    if (input$ls == 9){
      lslab <- "TempHumidity"
      lsparam <- c("A","a","b")
    }
    if (input$ls == 10){
      lslab <- "TempNonthermal"
      lsparam <- c("a","b","c")
    }
    if (input$ls == 11){
      lslab <- "Eyring3"
      lsparam <- c("a","b","c","d")
    }
    if (input$ls == 12){
      lslab <- "MultiStress"
      lsparam <- paste("a",0:3)
    }
    return(list(distlab,lslab,c(distparam,lsparam)))
  })

  # Compute the LSQ and MLE estimates
  LSQMLE <- reactive({
    xircSxiSrc <- sort.xircstressdata(rawdat())
    LSQoutput <- lifestress.LSQest(ls=distlslab()[[2]],dist=distlslab()[[1]],pp=plotparamout1())
    MLEoutput <- lifestress.MLEest(LSQoutput[[3]],ls=distlslab()[[2]],dist=distlslab()[[1]],xircSxiSrc[[1]],xircSxiSrc[[3]],xircSxiSrc[[2]],xircSxiSrc[[4]],input$conf/100)
    return(list(LSQoutput,MLEoutput))
  })

  # List of the LSQ estimates vs. MLE estimates
  outputlist2 <- reactive({
    if (input$pp == 1){
      LSQparamresult <- c("beta = ",LSQMLE()[[1]][[3]][1])
      MLEparamresult <- c("beta = ",LSQMLE()[[2]][[1]][1])
    }
    if (input$pp == 2){
      LSQparamresult <- c("beta = ",LSQMLE()[[1]][[3]][1])
      MLEparamresult <- c("beta = ",LSQMLE()[[2]][[1]][1])
    }
    if (input$pp == 3){
      LSQparamresult <- c("sigma = ",LSQMLE()[[1]][[3]][1])
      MLEparamresult <- c("sigma = ",LSQMLE()[[2]][[1]][1])
    }
    if (input$pp == 4){
      LSQparamresult <- c("sigma = ",LSQMLE()[[1]][[3]][1])
      MLEparamresult <- c("sigma = ",LSQMLE()[[2]][[1]][1])
    }
    if (input$pp == 6){
      LSQparamresult <- c("sigma = ",LSQMLE()[[1]][[3]][1])
      MLEparamresult <- c("sigma = ",LSQMLE()[[2]][[1]][1])
    }

    # Life-Stress Label
    if (input$ls == 1){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])
      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 2){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 3){
      if (input$pp == 5){
        LSQparamLSresult <- c("Ea = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("Ea = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", Ea = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", Ea = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 4){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 5){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 6){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 7){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 8){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])
      }
    }
    if (input$ls == 9){
      if (input$pp == 5){
        LSQparamLSresult <- c("A = ",LSQMLE()[[1]][[3]][1],", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c("A = ",LSQMLE()[[2]][[1]][1],", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3])

      } else{
        LSQparamLSresult <- c(", A = ",LSQMLE()[[1]][[3]][2],", a = ",LSQMLE()[[1]][[3]][3],", b = ",LSQMLE()[[1]][[3]][4])
        MLEparamLSresult <- c(", A = ",LSQMLE()[[2]][[1]][2],", a = ",LSQMLE()[[2]][[1]][3],", b = ",LSQMLE()[[2]][[1]][4])
      }
    }
    if (input$ls == 10){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2],", c = ",LSQMLE()[[1]][[3]][3])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2],", c = ",LSQMLE()[[2]][[1]][3])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3],", c = ",LSQMLE()[[1]][[3]][4])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3],", c = ",LSQMLE()[[2]][[1]][4])
      }
    }
    if (input$ls == 11){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2],", c = ",LSQMLE()[[1]][[3]][3],", d = ",LSQMLE()[[1]][[3]][4])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[3]][1],", b = ",LSQMLE()[[2]][[1]][2],", c = ",LSQMLE()[[2]][[1]][3],", d = ",LSQMLE()[[2]][[1]][4])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3],", c = ",LSQMLE()[[1]][[3]][4],", d = ",LSQMLE()[[1]][[3]][5])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3],", c = ",LSQMLE()[[2]][[1]][4],", d = ",LSQMLE()[[2]][[1]][5])
      }
    }
    if (input$ls == 12){
      if (input$pp == 5){
        LSQparamLSresult <- c("a = ",LSQMLE()[[1]][[3]][1],", b = ",LSQMLE()[[1]][[3]][2],", c = ",LSQMLE()[[1]][[3]][3],", d = ",LSQMLE()[[1]][[3]][4])
        MLEparamLSresult <- c("a = ",LSQMLE()[[2]][[1]][1],", b = ",LSQMLE()[[2]][[1]][2],", c = ",LSQMLE()[[2]][[1]][3],", d = ",LSQMLE()[[2]][[1]][4])

      } else{
        LSQparamLSresult <- c(", a = ",LSQMLE()[[1]][[3]][2],", b = ",LSQMLE()[[1]][[3]][3],", c = ",LSQMLE()[[1]][[3]][4],", d = ",LSQMLE()[[1]][[3]][5])
        MLEparamLSresult <- c(", a = ",LSQMLE()[[2]][[1]][2],", b = ",LSQMLE()[[2]][[1]][3],", c = ",LSQMLE()[[2]][[1]][4],", d = ",LSQMLE()[[2]][[1]][5])
      }
    }

    if (input$pp == 5){
      cat("LSQ: ",LSQparamLSresult,"\n")
      cat("MLE: ",MLEparamLSresult,"\n")
    } else {
      cat("LSQ: ",LSQparamresult,LSQparamLSresult,"\n")
      cat("MLE: ",MLEparamresult,MLEparamLSresult,"\n")
    }
  })

  # MLE estimates
  outputtable1 <- reactive({
    pointestvec<-LSQMLE()[[2]][[1]]
    lowvec<-rep(0,length(pointestvec))
    highvec<-rep(0,length(pointestvec))

    for(i in 1:length(pointestvec)){
      lowvec[i]<-LSQMLE()[[2]][[2]][[i]][1]
      highvec[i]<-LSQMLE()[[2]][[2]][[i]][2]
    }
    MLEtable<-matrix(c(distlslab()[[3]],pointestvec,lowvec,highvec),nrow=length(pointestvec),ncol=4,byrow=FALSE)
    colnames(MLEtable)<-c("Parameter","Point Estimate","Low Bound","Upper Bound")
    return(MLEtable)
  })

  # List of AFs
  outputlist3 <- reactive({
    Suse<-as.numeric(unlist(strsplit(input$usestress,",")))
    if (input$pp == 5){
      AFLSQ <- accelfactor(LSQMLE()[[1]][[3]],distlslab()[[2]],sort(LSQMLE()[[1]][[1]]),Suse)
      AFMLE <- accelfactor(LSQMLE()[[2]][[1]],distlslab()[[2]],sort(LSQMLE()[[1]][[1]]),Suse)
    } else {
      AFLSQ <- accelfactor(LSQMLE()[[1]][[3]][2:length(LSQMLE()[[1]])],distlslab()[[2]],sort(LSQMLE()[[1]][[1]]),Suse)
      AFMLE <- accelfactor(LSQMLE()[[2]][[1]][2:length(LSQMLE()[[1]])],distlslab()[[2]],sort(LSQMLE()[[1]][[1]]),Suse)
    }
    AFtable<-matrix(c(sort(LSQMLE()[[1]][[1]]),AFLSQ,AFMLE),nrow=length(LSQMLE()[[1]][[1]]),ncol=3,byrow=FALSE)
    #AFtable<-matrix(c(LSQMLE()[[1]][[1]],AFLSQ),nrow=length(LSQMLE()[[1]][[1]]),ncol=2,byrow=FALSE)
    colnames(AFtable)<-c("Stress","AF (LSQ)","AF (MLE)")
    #colnames(AFtable)<-c("Stress","AF (LSQ)")
    return(AFtable)
  })

  # Relationship Plot
  plotLSrelout1 <- reactive({
    # Prep the stresses (min, max, and use)
    minSinput <- as.vector(as.numeric(strsplit(input$stressmin,",")[[1]]))
    maxSinput <- as.vector(as.numeric(strsplit(input$stressmax,",")[[1]]))
    useSinput <- as.vector(as.numeric(strsplit(input$usestress,",")[[1]]))
    # Compute the plot and use-life (MLE and)
    LSplotandUselife<-lifestress.relationplot(rawdat(),ls=distlslab()[[2]],dist=distlslab()[[1]],pp=plotparamout1(),minSinput,maxSinput,useSinput,1,input$conf*0.01,input$xlabel,input$slabel)
    return(LSplotandUselife)
  })

  outputlist4 <- reactive({
    cat("Characteristic Life at Use Stress in ",input$xlabel," (LSQ):",plotLSrelout1()[[1]],"\n")
    cat("Characteristic Life at Use Stress in ",input$xlabel," (MLE):",plotLSrelout1()[[2]],"\n")
  })

  # OUTPUT 1: Restate the input data in table form
  output$rawdata <- renderTable({
    rawdat()
  })

  # OUTPUT 2: Probability plot which may be either Weibull, Lognormal,
  # Normal, or Exponential
  output$plot1 <- renderPlot({
    plotparamout1()
  })

  # OUTPUT 3: Summary of LSQ Plotting Parameters
  output$paramssum1 <- renderPrint({
    outputlist1()
  })

  # OUTPUT 4: Summary of Life-Stress LSQ and MLE parameters
  output$paramsLSQ <- renderPrint({
    outputlist2()
  })
  output$paramsMLE <- renderTable({
    outputtable1()
  })
  # OUTPUT 5: Acceleration Factor List
  output$AFtable <- renderTable({
    outputlist3()
  })
  # OUTPUT 6: Relationship Plot
  output$plot2 <- renderPlot({
    plotLSrelout1()
  })
  # OUTPUT 7: Summary of LSQ and MLE Use Life
  output$uselife <- renderPrint({
    outputlist4()
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
