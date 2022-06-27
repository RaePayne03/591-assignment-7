## Author: Rachel Payne
## rtpayne@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(tidyverse)
library(colourpicker)


# Define UI
ui <- fluidPage(
  titlePanel('BF591 Assignment 7'),
  p("Rachel Payne"),
  sidebarLayout(
    sidebarPanel(
      p('Load differential expression results'),
      fileInput('file1', 'Choose file to upload',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  '.csv'),
                placeholder = "deseq_results.csv"),
      radioButtons('xaxis', "choose column for the x axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue", "padj"),
                   selected = "log2FoldChange"), 
      
      radioButtons('yaxis', "choose column for the y axis",
                   choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue", "padj"),
                   selected = "pvalue"), 
      
      colourInput('accent_col',"Select threshold accent color","#FF5500"), 
      
      colourInput('base_col',"Select base color","#0C5E9C"),
      
      sliderInput('update_slider', "select the magnitude of p-adjusted coloring:",
                  value = -150,
                  min = -300, 
                  max = 0), 
      submitButton(text = "Apply Changes", width = "100%")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
        tabPanel("Table", tableOutput("results_table")))
    )
    
  )
)

# Define server logic
server <- function(input, output, session) {
  
  #' load_Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is 
  #' still a "function", but it will take no arguments. The `reactive({})` bit 
  #' says "if any of my inputs (as in, input$...) are changed, run me again". 
  #' This is useful when a user clicks a new button or loads a new file. In 
  #' our case, look for the uploaded file's datapath argument and load it with 
  #' read.csv. Return this data frame in the normal return() style.
  
  load_data<- reactive({
    inFile <- read.delim(file = input$file1$datapath, sep = ",", header = TRUE) %>%
      rename(ENSEMBL_ID = X)
    if (is.null(inFile))
      return(NULL)
    return(inFile)
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
  
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    volc_plot <- dataf %>% 
      select(c(x_name,y_name)) %>%
      mutate(volcano_status = ifelse(.[,y_name] <10^slider, 'TRUE', 'FALSE')) %>%
      ggplot()+ 
      geom_point(aes(x= !!sym(x_name),
                     y= -log10(!!sym(y_name)),
                     color = factor(volcano_status, levels = c("TRUE", "FALSE")))) +
      scale_color_manual(values = c(color1,color2)) + 
      theme_minimal()+ 
      theme(legend.position = "bottom") + 
      labs(x = x_name,
           y = paste0('-log10(',y_name,')'),
           color = paste0(y_name, ' < 1 x 10^',slider))
    
    return(volc_plot)
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
    new_table <- dataf %>%
      arrange(pvalue) %>%
      filter(padj < 10^slider) %>% 
      mutate(pvalue = formatC(.$pvalue, digits = 2, format = 'e'), 
             padj = formatC(.$padj, digits =2, format = 'e'))
    
    return(new_table)
  }
  
  #output results
  output$results_table<-renderTable({
    req(input$file1$datapath)
    
    draw_table(dataf = load_data(),
               slider = input$update_slider)
  })
  
  output$volcano <- renderPlot({
    req(input$file1$datapath)
    
    volcano_plot(dataf = load_data(),
                 x_name = input$xaxis,
                 y_name = input$yaxis,
                 slider = input$update_slider,
                 color1 = input$accent_col,
                 color2 = input$base_col)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


