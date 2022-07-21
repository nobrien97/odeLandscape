#############################################################################
# ODE Fitness Landscape App
# Author: Nick O'Brien
# Last changed: 17/7/2022
#############################################################################


library(shiny)
library(shinyjs)
library(shinyFiles)
library(tidyverse)
library(latex2exp)
library(cowplot)
library(rvest)
library(metR)

# Helper functions

molTraits <- list(
  "aZ",
  "bZ",
  "KZ",
  "KXZ"
)

molTraitsFigLab <- list(
  "aZ" = TeX("$\\alpha_Z$ value"),
  "bZ" = TeX("$\\beta_Z$ value"),
  "KZ" = TeX("$K_Z$ value"),
  "KXZ" = TeX("$K_{XZ}$ value")
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
  result <- read_csv("./landscape.csv", col_names = F, col_types = "d")
  names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  return(result)
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
  result[[xaxis]] <- seq(0, 10, by = 0.05)
  result[[yaxis]] <- seq(0, 10, by = 0.05)
  missing <- getMissingNames(xaxis, yaxis)
  
  result[[missing[[1]]]] <- missing1
  result[[missing[[2]]]] <- missing2

  return(result)
}

# Update missing trait names for our sliders
updateMissingValues <- function(input, session) {
  # Our current missing_value inputs (mol traits we are keeping track of)
  
  # Updated list of missing_values after changing x and y axes
  newNames <- getMissingNames(input$xaxisSelector, input$yaxisSelector)
  session$userData$missing_values <- unlist(newNames)

  updateSliderInput(session, "missing_value1", 
                    label = paste(newNames[[1]], "value"))
  updateSliderInput(session, "missing_value2", 
                    label = paste(newNames[[2]], "value"))

}

getMissingNames <- function(xaxis, yaxis) {
  inputs <- list(xaxis, yaxis)
  return(setdiff(molTraits, inputs))
  
}

nextTrait <- function(curTrait) {
  return(ifelse(curTrait != molTraits[[length(molTraits)]], 
         molTraits[[match(curTrait, molTraits) + 1]], molTraits[1]))
}

# Fetch simulation data with the matching seed and model index
fetchSimData <- function(file, seed, modelindex) {
  if (file.exists(file)) {
    # If we want a random seed, find one
    if (seed == "random")
      seed <- sample(system(sprintf("tail -n +2 %s | awk -F',' '{print $2}' | sort -u", 
                                    file), intern = T), 1)
    
    header <- "gen,seed,modelindex,mutType,mutID,position,constraint,originGen,value,chi,Freq,mutCount,fixGen,meanH,VA,phenomean,phenovar,dist,mean_w,deltaPheno,deltaW,Q1_AUC,Q1_aZ,Q1_bZ,Q1_KZ,Q1_KXZ,Q2_AUC,Q2_aZ,Q2_bZ,Q2_KZ,Q2_KXZ,Q3_AUC,Q3_aZ,Q3_bZ,Q3_KZ,Q3_KXZ"
    system(sprintf("awk -F , \'$2 == \"%s\" && $3 == %i\' %s > ./sim.csv",
                                   seed, modelindex, file))
    # Stick a header on
    system(sprintf("sed -i -e '1i%s' ./sim.csv", header))
    return(read_csv("./sim.csv", 
                    col_select = c(1:3, 16, 19, 28:31), 
                    col_types = "iffdddddd") %>% distinct())
  }
}

getRandomSeed <- function(file) {
  sample(system(sprintf("tail -n +2 %s | awk -F',' '{print $2}' | sort -u", 
                                file), intern = T), 1)
}

plotElements <- function(xaxis, yaxis, breaks) {
  list(
    labs(x = molTraitsFigLab[[xaxis]], y = molTraitsFigLab[[yaxis]]),
    geom_contour_filled(breaks = breaks),
    scale_x_continuous(limits = c(0, 10)),
    scale_y_continuous(limits = c(0, 10)),
    theme_bw(),
    theme(text = element_text(size = 16), legend.position = "bottom"),
    guides(fill = guide_legend(nrow=3, byrow = T))
    )
}

# Plot data
plotFigure <- function(df_landscape, input) {
  fitnesscc <- colorRampPalette(c("#D42B4F", "#2BD4B0"))(11)
  phenocc <- colorRampPalette(c("#D42B4F", "#2BD4B0"))(11)
  xaxis <- input$xaxisSelector
  yaxis <- input$yaxisSelector
  phenoPlot <- ggplot(df_landscape, aes_string(x = xaxis, 
                                    y = yaxis,
                                    z = "pheno")) +
    plotElements(xaxis, yaxis, breaks = seq(0, max(df_landscape$pheno), length.out = 11)) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = phenocc) +
    labs(fill = "Phenotypic\nvalue (Z)")
    
  
  fitnessPlot <- ggplot(df_landscape, aes_string(x = xaxis, 
                                          y = yaxis,
                                          z = "fitness")) +
    plotElements(xaxis, yaxis, breaks = seq(0, 1, length.out = 11)) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = fitnesscc) +
    labs(fill = "Fitness (w)")

  return(list(phenoPlot, fitnessPlot))
}

# Overlay simulation: returns a list of lists to stick onto the numeric solution
overlaySimData <- function(plt, df_sim, xaxis, yaxis) {
  xaxis <- paste0("Q2_", xaxis)
  yaxis <- paste0("Q2_", yaxis)
  list(list(
  geom_point(mapping = aes_string(x = xaxis, 
                                  y = yaxis,
                                  size = "gen"),
                       data = df_sim, inherit.aes = F),
    scale_size_continuous(range = c(0, 4), guide = "none")),
  list(
  geom_point(mapping = aes_string(x = xaxis, 
                                  y = yaxis, 
                                  size = "gen"),
                       data = df_sim, inherit.aes = F),
    scale_size_continuous(range = c(0, 4), guide = "none")
  ))
}



server <- function(input, output, session) {
  # volumes <- c(getVolumes()())
  # shinyFileChoose(input, "dataInput", roots = volumes, session = session)
  
  # Storage for current missing_value names
  session$userData$missing_values <- c("KZ", "KXZ")
    
  # Create temp directory for storing simulation info and landscape csv
  session$userData$rootTmpPath <- tempdir()
  session$userData$rootTmpPath <- file.path(session$userData$rootTmpPath, paste0("landscape", sample(1:.Machine$integer.max, 1)))
  if (!dir.exists(session$userData$rootTmpPath))
    dir.create(session$userData$rootTmpPath)
  setwd(session$userData$rootTmpPath)
  
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
    molTraitValues <- getMolTraitRanges(input$xaxisSelector, input$yaxisSelector,
                                        input$missing_value1, input$missing_value2)
    genInputFile(
      molTraitValues[["aZ"]], molTraitValues[["bZ"]], 
      molTraitValues[["KZ"]], molTraitValues[["KXZ"]]
    )
    # Calculate landscape
    session$userData$landscape <- runLandscaper("./samples.csv", 
                                                input$widthInput, 
                                                input$optInput, 
                                                input$threadsInput)
    
    session$userData$plotObject <- plotFigure(req(session$userData$landscape), input)

    buttonLocker(buttons)
  }, ignoreNULL = TRUE)
  
  observeEvent(input$randomButton, {
    updateTextInput(session, "seedInput", value = getRandomSeed(input$dataInput))
  }, ignoreNULL = TRUE)
  
  simData <- reactive({
    if (!isTruthy(input$dataInput) | !isTruthy(input$seedInput) | !isTruthy(input$modelindexInput)) {return()}
    #file_selected <- parseFilePaths(volumes, input$dataInput)$datapath
    file_selected <- input$dataInput
    fetchSimData(file_selected, input$seedInput, input$modelindexInput)
  })
  
  output$main_plot <- renderPlot({
    input$goButton
    if (!is.null(simData()))
      session$userData$overlayObject <- overlaySimData(req(session$userData$plotObject), req(simData()),
                                                    input$xaxisSelector, input$yaxisSelector)

    # Only plot if we have a plotObject, only stick on overlays if they exist
    if (!is.null(session$userData$plotObject)) {
      if (!is.null(session$userData$overlayObject)) {
        plot_grid(session$userData$plotObject[[1]] + session$userData$overlayObject[[1]], 
                  NULL, 
                  session$userData$plotObject[[2]] + session$userData$overlayObject[[2]], 
                  rel_widths = c(1, 0.1, 1), labels = c("A", "", "B"), nrow = 1)
      } else {
        plot_grid(session$userData$plotObject[[1]], 
                  NULL, 
                  session$userData$plotObject[[2]], 
                  rel_widths = c(1, 0.1, 1), labels = c("A", "", "B"), nrow = 1)
      }
    }
  })
  
}
