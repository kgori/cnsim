#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel("Copynumber simulation"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            # Samplesize
            numericInput("samplesize",
                         "Number of SNPs", 10000,
                         min = 0, max = 1000000, step = 100),
            # Depth
            numericInput("depth",
                         "Sequencing Depth", 100,
                         min = 0, max = 1000, step = 1),
            # T_a
            fluidRow(
                column(6, numericInput("ta", "Tumour A-allele copy number", 1,
                                       min = 0, max = 10, step = 1)),
                # T_b
                column(6, numericInput("tb", "Tumour B-allele copy number", 1,
                                       min = 0, max = 10, step = 1))
            ),
            # H_a
            fluidRow(
                column(6, numericInput("ha", "Host A-allele copy number", 1,
                                       min = 0, max = 10, step = 1)),
            # H_b
                column(6, numericInput("hb", "Host B-allele copy number", 1,
                                       min = 0, max = 10, step = 1))
            ),
            # Purity
            sliderInput("purity",
                        "Sample Purity",
                        min = 0, max = 1, value = 0.7,
                        animate = animationOptions(interval = 300, loop = TRUE),
                        step = 0.01, round = 2),
            
            # Ploidy
            sliderInput("ploidy",
                        "Average tumour ploidy",
                        min = 1, max = 4, value = 2.0, 
                        animate = animationOptions(interval = 300, loop = TRUE), 
                        step = 0.05, round = 2),
            
            # (Balance)
            numericInput("balance",
                         "VAF symmetry",
                         min = 0, max = 1, value = 0.5, 
                         step = 0.01),
            
            # (Host ploidy)
            
            # Histogram bins
            numericInput("bins",
                         "Bins",
                         min = 0, max = 200, value = 80, 
                         step = 1)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(plotOutput("vafPlot"), plotOutput("logRPlot"))
        )
    )
))
