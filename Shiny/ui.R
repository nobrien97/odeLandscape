#############################################################################
# ODE Fitness Landscape App
# Author: Nick O'Brien
# Last changed: 17/7/2022
#############################################################################
library(shiny)
library(shinyjs)
library(latex2exp)

molTraits <- list(
  "aZ",
  "bZ",
  "KZ",
  "KXZ"
)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Molecular trait fitness landscape visualiser"),  
  fluidRow(
    column(width = 4,
      fileInput("dataInput", "Upload data file", accept = ".csv")),
    column(width = 4,
      textInput("seedInput", "Replicate")),
    column(width = 4,
      textInput("modelindexInput", "Model Index")
    )
  ),
  fluidRow(
    column(width = 6,
      numericInput("widthInput", "Selection strength", 0.005, min = 0.0, 
                   max = NA, step = 0.001)
    ),
    column(width = 6,
      numericInput("optInput", "Optimum phenotype", 0, min = 0, 
                   max = NA, step = 1)
    ),
  ),
  fluidRow(
    column(width = 12,
      actionButton("goButton", "Go!", icon = icon("play"))       
    )
  ),
  fluidRow(
    hr(),
    column(width = 4,
           radioButtons("xaxisSelector", "X-axis molecular trait", 
                        choices = molTraits, selected = "aZ")),
    column(width = 4,
           radioButtons("yaxisSelector", "Y-axis molecular trait",
                        choices = molTraits, selected = "bZ")),
    column(width = 4,
           sliderInput("missing_value1", "aZ value", min = 0, 
                       max = 50, value = 1.0, step = 0.1),
           sliderInput("missing_value2", "bZ value", min = 0, 
                       max = 50, value = 1.0, step = 0.1)
           ),
  ),
  fluidRow(
    column(width = 12, align = "center",
      plotOutput(outputId = "main_plot", height = "800px", width = "600px"))
  )      
)

