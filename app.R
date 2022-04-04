## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)
library(tidyverse)
library(ggrepel)


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Assignment 7"),
  HTML("<p>To use this application, download the CSV <code>deseq_res.csv</code> from the data directory of this app's repository.</p>"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fileinput", "Load differential expression results",
              accept = ".csv", placeholder = "deseq_res.csv"),
    
      HTML(paste(rep("<p>A volcano plot can be generated with <b>'log<sub>2</sub> fold-change'</b> on the x-axis and <b>'p-adjusted'</b> on the y-axis.</p>"),
                 collapse = "")),
      
      radioButtons("radio_1", "Choose the column for the x-axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                   selected = "log2FoldChange"),
      radioButtons("radio_2", "Choose the column for the y-axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                   selected = "padj"),
      
      colourInput("col_1", label = "Base point color" ,"#22577A"),
      colourInput("col_2", label = "Highlight point color" ,"#FFCF56"),
      
      sliderInput(inputId = "slider_padj", min = -300, max = 0,
                  label = "Select the magnitude of the p adjusted coloring:", 
                  value = -150, step = 1),
      
      submitButton("Plot", width = '100%')
    ), # sidebarPanel
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("volcano_plot")),
        tabPanel("Table", tableOutput("data_table"))
        )
      ) # end of mainPanel
    
    ) # end of sidebarLayout
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
      req(input$fileinput)
      file <- input$fileinput
      if (is.null(file)) {
        return(NULL)
      } 
      else {
          data <- read.csv(file$datapath, header= TRUE, sep=",")
      }
      data <- data %>% rename(Gene = X)
      return(data)
      })
    
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    #' 
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
          plot <- ggplot(data = dataf, aes(x = !!sym(x_name), y=-log10(!!sym(y_name)))) +
            geom_point(aes(color = padj< 1*10^(slider))) +
            labs( color = str_glue('{y_name} 1 x 10^ {slider}')) +
            theme(legend.position = "bottom")
          return(plot)
        }
    
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
      #print(slider)
      dataf <- arrange(dataf, pvalue)
      dataf <- filter(dataf, between(padj, 0, 10 ^slider ))
      dataf <- mutate(dataf, pvalue = formatC(dataf$pvalue, digits = 2, format = "E" ),
               padj = formatC(dataf$padj, digits = 2, format = "E"))
      return(dataf)
    }
    
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    I.radio_1 <- reactive({
      if (is.null(input$radio_1)) return("log2FoldChange")
      input$radio_1
    })
    I.radio_2 <- reactive({
      if (is.null(input$radio_2)) return("padj")
      input$radio_2
    })
    I.slider_padj <- reactive({
      if (is.null(input$slider_padj)) return(-210)
      input$slider_padj
    })
    I.col_1 <- reactive({
      if (is.null(input$col_1)) return("blue")
      input$col_1
    })
    I.col_2 <- reactive({
      if (is.null(input$col_2)) return("green")
      input$col_2
    })
    
    
    output$volcano_plot <- renderPlot({
      volcano_plot(dataf = load_data(), slider = I.slider_padj(), 
                   x_name=I.radio_1(), y_name=I.radio_2(), 
                   color1=I.col_1(), color2=I.col_2())
      }) # replace this NULL
    # Same here, just return the table as you want to see it in the web page
    output$data_table <- renderTable({
      draw_table(dataf = load_data(), slider = I.slider_padj())
      }) # replace this NULL
}


# Run the application
shinyApp(ui = ui, server = server)

