# Load the Shiny library needed to run the app
library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("Nonparametric Modeling Toolkit"),

  sidebarLayout(
    sidebarPanel(
      helpText(
        "This tool is used for the purpose of generating several
        nonparametric representations of reliability functions
        from data including right-censored data.  To start,
        prepare your data in a CSV file with the data
        in the first column and the censored status (if any)
        in the second column.  Represent right-censored data
        as ''0'' and non-censored as ''1''."
        ),
      # Input panel for Data
      fileInput(
        "datafile", "Choose CSV File", multiple = TRUE,
        accept = c("text/csv","text/comma-separated-values,text/plain",
                   ".csv")
      ),
      # Turn off censored data box
      checkboxInput(
        "censoroff", "Omit Censored Data", value = FALSE
        ),
      # Input panel for X-Label
      textInput(
        "xlabel", h4("X Label"),
        value = "X"
      ),
      # Input panel for Y-Label
      selectInput(
        "relfcn", h4("Reliability Function"),
        choices = list("Unreliability" = 1,  "Reliability" = 2,
                       "Hazard" = 3, "Cumulative Hazard" = 4),
        selected = 1
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
      # Compare plotting position
      selectInput(
        "plotposit2", h4("Compare with other Plotting Position"),
        choices = list(
          "Blom" = 1, "Mean" = 2, "Median" = 3, "Midpoint" = 4,
          "Jenkinson (Beard)" = 5, "Benard and Bos-Levenbach" = 6,
          "Tukey" = 7, "Gringorten" = 8, "Kaplan-Meier" = 9,
          "Nelson-Aalen" = 10, "None" = 11),
        selected = 11
      ),
      # Input panel for confidence
      numericInput(
        "conf", h4("Confidence Bounds (%)"), min = 1, max = 99.99, value = 95
      )

    ),
    mainPanel(
      # Output results in tabular form:
      # Tab 1) Unreliability or reliability plot
      # Tab 2) Table data: x, F, R, F confidence, and R confidence
      tabsetPanel(
        type = "tabs",
        tabPanel("Input data", tableOutput("rawdata")),
        tabPanel("Nonparametric plot", plotOutput("nonparamplot"), downloadButton("downloadnonparamplot", "Download")),
        tabPanel("Nonparametric output data", tableOutput("nonparamdata"), downloadButton("downloadnonparamdata", "Download")),
        tabPanel("Nonparametric output data (compare)", tableOutput("nonparamdata2"), downloadButton("downloadnonparamdata2", "Download"))
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


  # Pulls the rank, xi, and rc data from the table
  xirc <- reactive({
    # Check if there is right censored data or not
    if(is.null(dim(rawdat()))) {
      xi <- rawdat()[,1]
      i <- rankcalc(xi)
      xitab<-list(i, xi, rep(1, length(xi)))
    }
    else {

      xi <- rawdat()[,1][which(rawdat()[,2] == 1)]
      rc <- rawdat()[,1][which(rawdat()[,2] == 0)]
      i <- rankcalc(xi,rc)
      xitab<-list(i,xi,rc)
      if(input$censoroff==FALSE){
        xi <- rawdat()[,1][which(rawdat()[,2] == 1)]
        rc <- rawdat()[,1][which(rawdat()[,2] == 0)]
        i <- rankcalc(xi,rc)
        xitab<-list(i,xi,rc)
      } else {
        xi <- rawdat()[,1][which(rawdat()[,2] == 1)]
        rc <- rawdat()[,2][which(rawdat()[,2] == 1)]
        i <- rankcalc(xi,rc)
        xitab<-list(i,xi,rc)
      }
    }
    return(xitab)
  })

  # Precompute the plotting position output for each plotting position option
  plotpositoutput <- reactive({
    xitab<-xirc()
    i<-xitab[[1]]
    xi<-xitab[[2]]
    rc<-xitab[[3]]
    if (input$plotposit == 1){
      pp<-plotposit.blom(i, xi, rc)
      plotlab1 = "Blom"
    }
    if (input$plotposit == 2){
      pp<-plotposit.mean(i, xi, rc)
      plotlab1 = "Mean"
    }
    if (input$plotposit == 3){
      pp<-plotposit.median(i, xi, rc)
      plotlab1 = "Median"
    }
    if (input$plotposit == 4){
      pp<-plotposit.midpt(i, xi, rc)
      plotlab1 = "Midpoint"
    }
    if (input$plotposit == 5){
      pp<-plotposit.beard(i, xi, rc)
      plotlab1 = "Jenkinson (Beard)"
    }
    if (input$plotposit == 6){
      pp<-plotposit.bernbos(i, xi, rc)
      plotlab1 = "Benard and Bos-Levenbach"
    }
    if (input$plotposit == 7){
      pp<-plotposit.tukey(i, xi, rc)
      plotlab1 = "Tukey"
    }
    if (input$plotposit == 8){
      pp<-plotposit.grigorten(i, xi, rc)
      plotlab1 = "Gringorten"
    }
    if (input$plotposit == 9){
      pp<-plotposit.kaplanmeier(xi, rc)
      plotlab1 = "Kaplan-Meier"
    }
    if (input$plotposit == 10){
      pp<-plotposit.nelsonaalen(xi, rc)
      plotlab1 = "Nelson-Aalen"
    }
    #colnames(pp)<-c(input$xlabel,"Unreliability","Reliability","Hazard","Cumulative Hazard")
    return(list(pp,plotlab1))
  })

  # Precompute the plotting position output for each plotting position option
  plotpositoutput2 <- reactive({
    xitab<-xirc()
    i<-xitab[[1]]
    xi<-xitab[[2]]
    rc<-xitab[[3]]
    if (input$plotposit2 == 1){
      pp<-plotposit.blom(i, xi, rc)
      plotlab1 = "Blom"
    }
    if (input$plotposit2 == 2){
      pp<-plotposit.mean(i, xi, rc)
      plotlab1 = "Mean"
    }
    if (input$plotposit2 == 3){
      pp<-plotposit.median(i, xi, rc)
      plotlab1 = "Median"
    }
    if (input$plotposit2 == 4){
      pp<-plotposit.midpt(i, xi, rc)
      plotlab1 = "Midpoint"
    }
    if (input$plotposit2 == 5){
      pp<-plotposit.beard(i, xi, rc)
      plotlab1 = "Jenkinson (Beard)"
    }
    if (input$plotposit2 == 6){
      pp<-plotposit.bernbos(i, xi, rc)
      plotlab1 = "Benard and Bos-Levenbach"
    }
    if (input$plotposit2 == 7){
      pp<-plotposit.tukey(i, xi, rc)
      plotlab1 = "Tukey"
    }
    if (input$plotposit2 == 8){
      pp<-plotposit.grigorten(i, xi, rc)
      plotlab1 = "Gringorten"
    }
    if (input$plotposit2 == 9){
      pp<-plotposit.kaplanmeier(xi, rc)
      plotlab1 = "Kaplan-Meier"
    }
    if (input$plotposit2 == 10){
      pp<-plotposit.nelsonaalen(xi, rc)
      plotlab1 = "Nelson-Aalen"
    }
    if (input$plotposit2 == 11){
      pp<-xi
      plotlab1 = "None"
    }
    return(list(pp,plotlab1))
  })

  # Precompute confidence bounds
  plotconfoutput <- reactive({
    if(input$censoroff==FALSE) {
      n<-length(rawdat()[,1])
      #if(input$plotposit == 9){
      #  conlim<-conf.greenwood((100-input$conf)*0.01,n,plotpositoutput()[,3],count.failcen(rawdat(),n))
      #} else {
      #  conlim<-conf.binomial((100-input$conf)*0.01,n,plotpositoutput()[,3])
      #}
      conlim<-conf.greenwood((100-input$conf)*0.01,n,plotpositoutput()[[1]][,3],count.failcen(rawdat(),n))
    } else{
      n<-sum(rawdat()[,2])
      trimdat<-matrix(c(rawdat()[,1][which(rawdat()[,2]==1)],rep(1,sum(rawdat()[,2]))),nrow=sum(rawdat()[,2]),ncol=2)
      conlim<-conf.greenwood((100-input$conf)*0.01,n,plotpositoutput()[[1]][,3],count.failcen(trimdat,n))
    }
    return(conlim)
  })

  # Precompute confidence bounds 2
  plotconfoutput2 <- reactive({
    if(input$censoroff==FALSE) {
      n<-length(rawdat()[,1])
      conlim<-conf.greenwood((100-input$conf)*0.01,n,plotpositoutput2()[[1]][,3],count.failcen(rawdat(),n))
    } else {
      n<-sum(rawdat()[,2])
      trimdat<-matrix(c(rawdat()[,1][which(rawdat()[,2]==1)],rep(1,sum(rawdat()[,2]))),nrow=sum(rawdat()[,2]),ncol=2)
      conlim<-conf.greenwood((100-input$conf)*0.01,n,plotpositoutput2()[[1]][,3],count.failcen(trimdat,n))
    }
    return(conlim)
  })

  nonparamtable1 <- reactive({
    set1<-plotpositoutput()[[1]]
    set2<-plotconfoutput()
    nonparamset<-matrix(
      c(set1[,1],set1[,2],set2[[3]],set2[[4]],set1[,3],
        set2[[1]],set2[[2]],set1[,4],set2[[7]],set2[[8]],set1[,5],set2[[5]],set2[[6]]),
      nrow = length(set1[,1]), ncol = 13)
    colnames(nonparamset)<-c(input$xlabel,"Unreliability","Low bound","High bound","Reliability","Low bound","High bound","Hazard","Low bound","High bound","Cumulative Hazard","Low bound","High bound")
    return(nonparamset)
  })

  nonparamtable2 <- reactive({
    set1<-plotpositoutput2()[[1]]
    set2<-plotconfoutput2()
    nonparamset<-matrix(
      c(set1[,1],set1[,2],set2[[3]],set2[[4]],set1[,3],
        set2[[1]],set2[[2]],set1[,4],set2[[7]],set2[[8]],set1[,5],set2[[5]],set2[[6]]),
      nrow = length(set1[,1]), ncol = 13)
    colnames(nonparamset)<-c(input$xlabel,"Unreliability","Low bound","High bound","Reliability","Low bound","High bound","Hazard","Low bound","High bound","Cumulative Hazard","Low bound","High bound")
    return(nonparamset)
  })

  ylabelout <- reactive({
    if (input$relfcn == 1){
      ylabel1 = "Unreliability"
    }
    if (input$relfcn == 2){
      ylabel1 = "Reliability"
    }
    if (input$relfcn == 3){
      ylabel1 = "Hazard Function"
    }
    if (input$relfcn == 4){
      ylabel1 = "Cumulative Hazard Function"
    }
    return(ylabel1)
  })

  fullplot <- function(){
    # Plot for unreliability
    if (input$relfcn ==1){
      medplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotpositoutput()[[1]][,2])
      lowplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotconfoutput()[[3]])
      hiplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotconfoutput()[[4]])
      ylims<-c(0,1)
      # Secondary comparison plot
      if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
          input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
          input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
          input$plotposit2 == 10){
        medplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotpositoutput2()[[1]][,2])
        lowplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotconfoutput2()[[3]])
        hiplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotconfoutput2()[[4]])
      }
    }
    # Plot for reliability
    if (input$relfcn ==2){
      medplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotpositoutput()[[1]][,3])
      lowplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotconfoutput()[[1]])
      hiplot<-plot.nonparam(plotpositoutput()[[1]][,1], plotconfoutput()[[2]])
      ylims<-c(0,1)
      # Secondary comparison plot
      if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
          input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
          input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
          input$plotposit2 == 10){
        medplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotpositoutput2()[[1]][,3])
        lowplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotconfoutput2()[[1]])
        hiplot2<-plot.nonparam(plotpositoutput2()[[1]][,1], plotconfoutput2()[[2]])
      }
    }
    # Plot for hazard
    if (input$relfcn ==3){
      medplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotpositoutput()[[1]][,4])
      lowplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotconfoutput()[[7]])
      hiplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotconfoutput()[[8]])
      if (is.finite(max(hiplot[[2]]))){
        ylims<-c(0,max(hiplot[[2]]))
      } else {
        ylims<-c(0,hiplot[[2]][length(hiplot[[2]])-2])
      }
      # Secondary comparison plot
      if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
          input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
          input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
          input$plotposit2 == 10){
        medplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotpositoutput2()[[1]][,4])
        lowplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotconfoutput2()[[7]])
        hiplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotconfoutput2()[[8]])
        if (is.finite(max(c(max(hiplot[[2]]),max(hiplot2[[2]]))))) {
          ylims<-c(0,max(c(max(hiplot[[2]]),max(hiplot2[[2]]))))
        } else {
          ylims<-c(0,max(c(hiplot[[2]][length(hiplot[[2]])-2],hiplot2[[2]][length(hiplot2[[2]])-2])))
        }
      }
    }
    # Plot for cumulative hazard
    if (input$relfcn ==4){
      medplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotpositoutput()[[1]][,5])
      lowplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotconfoutput()[[5]])
      hiplot<-plot.nonparamhaz(plotpositoutput()[[1]][,1], plotconfoutput()[[6]])
      if (is.finite(max(hiplot[[2]]))){
        ylims<-c(0,max(hiplot[[2]]))
      } else {
        ylims<-c(0,hiplot[[2]][length(hiplot[[2]])-2])
      }
      # Secondary comparison plot
      if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
          input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
          input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
          input$plotposit2 == 10){
        medplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotpositoutput2()[[1]][,5])
        lowplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotconfoutput2()[[5]])
        hiplot2<-plot.nonparamhaz(plotpositoutput2()[[1]][,1], plotconfoutput2()[[6]])
        if (is.finite(max(c(max(hiplot[[2]]),max(hiplot2[[2]]))))) {
          ylims<-c(0,max(c(max(hiplot[[2]]),max(hiplot2[[2]]))))
        } else {
          ylims<-c(0,max(c(hiplot[[2]][length(hiplot[[2]])-2],hiplot2[[2]][length(hiplot2[[2]])-2])))
        }
        #ylims<-c(0,max(c(max(hiplot[[2]]),max(hiplot2[[2]]))))
      }
    }
    plot(medplot[[1]], medplot[[2]], type="l", xlab=input$xlabel, ylab=ylabelout(),ylim=ylims , col="blue")
    lines(medplot[[1]],lowplot[[2]], type="l", lty=2, col="blue")
    lines(medplot[[1]],hiplot[[2]], type="l", lty=2, col="blue")
    if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
        input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
        input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
        input$plotposit2 == 10){
      lines(medplot2[[1]],medplot2[[2]], type="l", col="green")
      lines(medplot2[[1]],lowplot2[[2]], type="l", lty=2, col="green")
      lines(medplot2[[1]],hiplot2[[2]], type="l", lty=2, col="green")
      if(input$relfcn ==2){
        legend( "bottomleft", c( plotpositoutput()[[2]], plotpositoutput2()[[2]]),text.col=c("blue", "green") )
      } else {
        legend( "topleft", c( plotpositoutput()[[2]], plotpositoutput2()[[2]]),text.col=c("blue", "green") )
      }
    }
  }

  # OUTPUT 1: Restate the input data in table form
  output$rawdata <- renderTable({
    rawdat()
  })

  # OUTPUT 2: Non-parametric reliability function based on
  # (1) plotting position input (and comparison if any), (2)
  # confidence interval, and (3) reliability function
  output$nonparamplot <- renderPlot({
    fullplot()
  })

  # OUTPUT 3: Table of non-parametric data output
  output$nonparamdata <- renderTable({
    nonparamtable1()
  })

  # OUTPUT 4: Table of non-parametric data output
  output$nonparamdata2 <- renderTable({
    if (input$plotposit2 == 1 | input$plotposit2 == 2 | input$plotposit2 == 3|
        input$plotposit2 == 4|input$plotposit2 == 5|input$plotposit2 == 6|
        input$plotposit2 == 7|input$plotposit2 == 8|input$plotposit2 == 9|
        input$plotposit2 == 10){
      nonparamtable2()
    }
  })

  # DOWNLOAD 1: The non-parametric output data
  output$downloadnonparamdata <- downloadHandler(
    filename = function() {
      paste("outputdata1", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(nonparamtable1(), file, row.names = FALSE)
    }
  )

  output$downloadnonparamdata2 <- downloadHandler(
    filename = function() {
      paste("outputdata2", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(nonparamtable2(), file, row.names = FALSE)
    }
  )
  # DOWNLOAD 2: The plot
  output$downloadnonparamplot <- downloadHandler(
    filename = "outputplot.png",
    content = function(file) {
      png(file)
      fullplot()
      dev.off()
    })

}

# Run the app ----
shinyApp(ui = ui, server = server)
