#############################################################################
# ODE Fitness Landscape App
# Author: Nick O'Brien
# Last changed: 17/7/2022
#############################################################################


library(tidyverse)
library(shiny)
library(shinyjs)
library(rvest)

# Helper functions

molTraits <- list(
  "aZ",
  "bZ",
  "KZ",
  "KXZ"
)


# Lock buttons when we hit go
buttonLocker <- function(buttons) {
  for (button in buttons)
    toggleState(button)
} 
buttons <- c("dataInput", "seedInput", "modelindexInput", "widthInput", 
             "optInput", "goButton", "xaxisSelector", "yaxisSelector",
             "missingAxis1", "missingAxis2")


# Run landscaper program
runLandscaper <- function(df_path, width, optimum, threads) {
  system(sprintf("ODELandscaper -i %s -o ./landscape.csv -w %f -p %f -t %i",
                 df_path, width, optimum, threads))
  return(read_csv("./landscape.csv", col_names = F, col_types = "double"))
}

# Generate input file for landscaper
genInputFile <- function(aZ, bZ, KZ, KXZ) {
  # Generate a table of molecular trait combinations
  samples <- expand.grid(aZ, bZ, KZ, KXZ)
  write.table(samples, "./samples.csv", sep = ",", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  return()
}

# Determine which input is pointing to which molecular trait
getMolTraitRanges <- function(xaxis, yaxis, missing1, missing2) {
  result <- list()
  result[[xaxis]] <- seq(0, 50, by = 0.05)
  result[[yaxis]] <- seq(0, 50, by = 0.05)
  inputs <- list(xaxis, yaxis)
  missing <- setdiff(molTraits, inputs)
  
  result[[missing[1]]] <- missing1
  result[[missing[2]]] <- missing2

  return(result)
}

# Fetch simulation data with the matching seed and model index
fetchSimData <- function(file, seed, modelindex) {
  return(read_csv(system(sprintf("awk -F , \'$2 == \"%s\" && $3 == %i\' %s",
                 seed, modelindex, file))))
}

# Update missing trait names for our sliders
updateMissingValues <- function(input, session) {
  # Our current missing_value inputs (mol traits we are keeping track of)
  
  # Updated list of missing_values after changing x and y axes
  newNames <- getMissingNames(input$xaxisSelector, input$yaxisSelector)

  index_x <- match(input$xaxisSelector, session$userData$missing_values)
  index_y <- match(input$yaxisSelector, session$userData$missing_values)
    
  if (!is.na(index_x)) {
    updateSliderInput(session, paste0("missing_value", index_x), label = paste(newNames[[index_x]], "value"))
    session$userData$missing_values[index_x] <- newNames[[index_x]] 
  } else {
    session$userData$missing_values[1] <- input$xaxisSelector
    }
  
  if (!is.na(index_y)) {
    updateSliderInput(session, paste0("missing_value", index_y), label = paste(newNames[[index_y]], "value"))
    session$userData$missing_values[index_y] <- newNames[[index_y]] 
  } else {
    session$userData$missing_values[2] <- input$yaxisSelector
  }

}

getMissingNames <- function(xaxis, yaxis) {
  inputs <- list(xaxis, yaxis)
  return(setdiff(molTraits, inputs))
  
}

nextTrait <- function(curTrait) {
  return(ifelse(curTrait != molTraits[[length(molTraits)]], 
         molTraits[[match(curTrait, molTraits) + 1]], molTraits[1]))
}



server <- function(input, output, session) {
  # user variables for storing current values of input$missing_value1 and 2
  session$userData$missing_values <- c("aZ", "bZ")

  # Create temp directory for storing simulation info and landscape csv
  session$userData$rootTmpPath <- tempdir()
  session$userData$rootTmpPath <- file.path(session$userData$rootTmpPath, paste0("landscape", sample(1:.Machine$integer.max, 1)))
  if (!dir.exists(session$userData$rootTmpPath))
    dir.create(session$userData$rootTmpPath)
  
  observe({
    enable(selector = "[name=yaxisSelector]")
    disable(selector = sprintf("[name=yaxisSelector][value=%s]", input$xaxisSelector))
    if (input$xaxisSelector == input$yaxisSelector) {
      updateRadioButtons(session, "yaxisSelector", 
                         selected = nextTrait(input$yaxisSelector))
    }
    updateMissingValues(input, session)
  }) 

  observeEvent(input$goButton, {
    # When we press the button, we first need to lock all inputs
    buttonLocker(buttons)
    # Generate an input to the landscaper
    inputs <- 
    session$userData$landscape_input <- genInputFile(
      input
    )
  }, ignoreNULL = TRUE)
}
