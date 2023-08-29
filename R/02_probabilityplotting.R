# Load the Shiny library needed to run the app
library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("Probability Plotting Toolkit"),

  sidebarLayout(
    sidebarPanel(
      helpText(
        "This tool is used for probability plotting of non-parametric data.
        The data may be either single-stress or multi-stress depending on
        user input.  To start, prepare your data in a CSV file with the data
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
      # Input panel for X-Label
      textInput(
        "xlabel", h4("X Label"),
        value = "X"
      ),
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
      # Input panel for Probability Plot
      selectInput(
        "pp1", h4("Probability Plot"),
        choices = list("Weibull" = 1,  "3P Weibull" = 2, "Lognormal" = 3,
                       "Normal" = 4, "Exponential" = 5,
                       "2P Exponential" = 6),
        selected = 1
      ),
      # Compare plotting position
      selectInput(
        "pp2", h4("Compare with other Probability Plot"),
        choices = list("None" = 1, "Weibull" = 2,  "3P Weibull" = 3, "Lognormal" = 4,
                       "Normal" = 5, "Exponential" = 6,
                       "2P Exponential" = 7),
        selected = 1
      )

    ),
    mainPanel(
      # Output results in tabular form:
      # Tab 1) Unreliability or reliability plot
      # Tab 2) Table data: x, F, R, F confidence, and R confidence
      tabsetPanel(
        type = "tabs",
        tabPanel("Input data", tableOutput("rawdata")),
        tabPanel(
          "Probabaility Plot 1",
          plotOutput("plot1"),
          verbatimTextOutput("paramssum1"),
          downloadButton("downloadprobplot", "Download")
          ),
        tabPanel(
          "Probabaility Plot 2",
          plotOutput("plot2"),
          verbatimTextOutput("paramssum2"),
          downloadButton("downloadprobplot2", "Download")
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
    if (input$pp1 == 1){
      plotandparams<-probplot.wbl(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp1 == 2){
      plotandparams<-probplot.wbl3P(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp1 == 3){
      plotandparams<-probplot.logn(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp1 == 4){
      plotandparams<-probplot.nor(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp1 == 5){
      plotandparams<-probplot.exp(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp1 == 6){
      plotandparams<-probplot.exp2P(rawdat(),plotpositlabel(),input$xlabel)
    }
    return(plotandparams)
  })

  # Process the other plot and LSQ parameters based on Probability Plot type (if available)
  plotparamout2 <- reactive({
    if (input$pp2 == 1){
      plotandparams<-1
    }
    if (input$pp2 == 2){
      plotandparams<-probplot.wbl(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp2 == 3){
      plotandparams<-probplot.wbl3P(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp2 == 4){
      plotandparams<-probplot.logn(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp2 == 5){
      plotandparams<-probplot.nor(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp2 == 6){
      plotandparams<-probplot.exp(rawdat(),plotpositlabel(),input$xlabel)
    }
    if (input$pp2 == 7){
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
      # Single Stress
      if (input$pp1 == 1){
        cat("alpha = ",plotparamout1()[[1]][[1]],", beta = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp1 == 2){
        cat("alpha = ",plotparamout1()[[1]][[1]],", beta = ",plotparamout1()[[1]][[2]],", gamma = ",plotparamout1()[[1]][[3]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp1 == 3){
        cat("log-mu = ",plotparamout1()[[1]][[1]],", log-sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
        }
      if (input$pp1 == 4){
        cat("mu = ",plotparamout1()[[1]][[1]],", sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout1()[[2]],"\n")
        }
      if (input$pp1 == 5){
        cat("lambda = ",plotparamout1()[[1]],", R^2 = ",plotparamout1()[[2]],"\n")
      }
      if (input$pp1 == 6){
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
        if (input$pp1 == 1){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout1()[[i*3-1]][1],", beta = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp1 == 2){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout1()[[i*3-1]][1],", beta = ",plotparamout1()[[i*3-1]][2],", gamma = ",plotparamout1()[[i*3-1]][3],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp1 == 3){
          cat("Stress ",i,stressrow,":log-mu = ",plotparamout1()[[i*3-1]][1],", log-sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp1 == 4){
          cat("Stress ",i,stressrow,":mu = ",plotparamout1()[[i*3-1]][1],", sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp1 == 5){
          cat("Stress ",i,stressrow,":lambda = ",plotparamout1()[[i*3-1]],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
        if (input$pp1 == 6){
          cat("Stress ",i,stressrow,":theta = ",plotparamout1()[[i*3-1]][1],", sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
      }
    }
  })

  outputlist2 <- reactive({
    databystress<-checkstress(rawdat())
    if (!is.null(dim(databystress))){
      # Single Stress
      if (input$pp2 == 2){
        cat("alpha = ",plotparamout2()[[1]][[1]],", beta = ",plotparamout2()[[1]][[2]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
      if (input$pp2 == 3){
        cat("alpha = ",plotparamout2()[[1]][[1]],", beta = ",plotparamout2()[[1]][[2]],", gamma = ",plotparamout1()[[1]][[3]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
      if (input$pp2 == 4){
        cat("log-mu = ",plotparamout2()[[1]][[1]],", log-sigma = ",plotparamout2()[[1]][[2]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
      if (input$pp2 == 5){
        cat("mu = ",plotparamout2()[[1]][[1]],", sigma = ",plotparamout2()[[1]][[2]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
      if (input$pp2 == 6){
        cat("lambda = ",plotparamout2()[[1]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
      if (input$pp2 == 7){
        cat("theta = ",plotparamout1()[[1]][[1]],", sigma = ",plotparamout1()[[1]][[2]],", R^2 = ",plotparamout2()[[2]],"\n")
      }
    } else {
      # Multi Stress
      for(i in 1:length(databystress)){
        stressrow<-rep(1,2*length(multistresslist()))

        for(j in 1:length(plotparamout2()[[i*3-2]])){
          stressrow[j*2-1]<-multistresslist()[j]
          stressrow[j*2]<-plotparamout2()[[i*3-2]][[j]]
        }
        if (input$pp2 == 2){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout2()[[i*3-1]][1],", beta = ",plotparamout2()[[i*3-1]][2],", R^2 = ",plotparamout2()[[i*3]],"\n")
        }
        if (input$pp2 == 3){
          cat("Stress ",i,stressrow,":alpha = ",plotparamout2()[[i*3-1]][1],", beta = ",plotparamout1()[[i*3-1]][2],", gamma = ",plotparamout2()[[i*3-1]][3],", R^2 = ",plotparamout2()[[i*3]],"\n")
        }
        if (input$pp2 == 4){
          cat("Stress ",i,stressrow,":log-mu = ",plotparamout2()[[i*3-1]][1],", log-sigma = ",plotparamout2()[[i*3-1]][2],", R^2 = ",plotparamout2()[[i*3]],"\n")
        }
        if (input$pp2 == 5){
          cat("Stress ",i,stressrow,":mu = ",plotparamout2()[[i*3-1]][1],", sigma = ",plotparamout2()[[i*3-1]][2],", R^2 = ",plotparamout2()[[i*3]],"\n")
        }
        if (input$pp2 == 6){
          cat("Stress ",i,stressrow,":lambda = ",plotparamout2()[[i*3-1]],", R^2 = ",plotparamout2()[[i*3]],"\n")
        }
        if (input$pp2 == 7){
          cat("Stress ",i,stressrow,":theta = ",plotparamout1()[[i*3-1]][1],", sigma = ",plotparamout1()[[i*3-1]][2],", R^2 = ",plotparamout1()[[i*3]],"\n")
        }
      }
    }
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

  # OUTPUT 3: Second probability plot which may be either Weibull,
  # Lognormal, Normal, or Exponential
  output$plot2 <- renderPlot({
    plotparamout2()
  })

  # OUTPUT 4: Summary of LSQ Parameters
  output$paramssum1 <- renderPrint({
    outputlist1()
  })


  # OUTPUT 5: Table of LSQ Parameters
  output$paramssum2 <- renderPrint({
    outputlist2()
  })

  # DOWNLOAD 1: The plots
  output$downloadprobplot <- downloadHandler(
    filename = "probabilityplot1.png",
    content = function(file) {
      png(file)
      plotparamout1()
      dev.off()
    }
  )

  output$downloadprobplot2 <- downloadHandler(
    filename = "probabilityplot2.png",
    content = function(file) {
      png(file)
      plotparamout2()
      dev.off()
    }
  )

  # DOWNLOAD: The LSQ Parameter Estimates
  output$downloadLSQparams1 <- downloadHandler(
    filename = "LSQ1.txt",
    content = function(file) {
      txt(file)
      plotparamout1()
      dev.off()
    }
  )

  output$downloadLSQparams2 <- downloadHandler(
    filename = "probabilityplot2.png",
    content = function(file) {
      png(file)
      plotparamout2()
      dev.off()
    }
  )

}

# Run the app ----
shinyApp(ui = ui, server = server)
