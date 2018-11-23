#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

replace_na <- function(val, default) {
    ifelse(is.na(val), default, val)
}

simulated_vaf <- function(samplesize, depth, ta, tb, ha, hb, purity, tploidy = 2, hploidy = 2, balance = 0) {
    stopifnot(purity >= 0 & purity <= 1)
    
    # Adjust the sequencing depth of this region from the average `depth`
    # based on this region's total copy number and the tumour and host average ploidies
    region_depth <- (((ta+tb) / tploidy) * purity + (1 - purity) * ((ha+hb) / hploidy)) * depth
    
    # Work out the expected number of reads based on this region's depth
    expected_ta <- replace_na(region_depth * purity * ta / (ta+tb), 0)
    expected_tb <- replace_na(region_depth * purity * tb / (ta+tb), 0)
    expected_ha <- replace_na(region_depth * (1-purity) * ha / (ha+hb), 0)
    expected_hb <- replace_na(region_depth * (1-purity) * hb / (ha+hb), 0)
    
    # Sample read counts from a Poisson distribution (smooth out zeros by adding a small weight)
    nta <- rpois(samplesize, expected_ta+0.01)
    ntb <- rpois(samplesize, expected_tb+0.01)
    nha <- rpois(samplesize, expected_ha+0.01)
    nhb <- rpois(samplesize, expected_hb+0.01)
    
    vaf <- ((ntb+nhb) / (nta+ntb+nha+nhb))
    expected_vaf <- ((expected_tb+expected_hb) / (expected_ta+expected_tb+expected_ha+expected_hb))
    
    # Rebalance the SNPs onto major and minor alleles
    ix <- seq(1, balance * length(vaf))
    vaf[ix] <- 1 - vaf[ix]
    list(vaf = vaf, expected_vaf = expected_vaf)
}

simulated_logr <- function(samplesize, depth, tn, hn, purity, tploidy = 2, hploidy = 2) {
    stopifnot(purity >= 0 & purity <= 1)
    
    # Work out the expected number of reads in this segment in a host sample
    # sequenced to equal depth
    expected_h <- (hn / hploidy) * depth
    
    # Work out the expected number of reads mapping to this segment in the tumour sample,
    # given the average depth, ploidy and purity
    expected_t <- ((tn / tploidy) * purity + (1 - purity) * (hn / hploidy)) * depth
    
    # Sample read counts from a Poisson distribution
    nt <- rpois(samplesize, expected_t)
    nh <- rpois(samplesize, expected_h)
    logr <- log2(nt/nh)
    expected_logr <- log2(expected_t / expected_h)
    list(logr = logr, expected_logr = expected_logr)
}


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$vafPlot <- renderPlot({
        
        # generate bins based on input$bins from ui.R
        vaf  <- simulated_vaf(input$samplesize,
                              input$depth,
                              input$ta,
                              input$tb,
                              input$ha,
                              input$hb,
                              input$purity,
                              tploidy = input$ploidy,
                              balance = input$balance)
        bins <- seq(0, 1, length.out = input$bins + 1)
        
        # draw the histogram with the specified number of bins
        hist(vaf$vaf, breaks = bins, col = "steelblue1", border = "steelblue1", probability = TRUE,
             main = "Histogram of VAF", xlab = "VAF")
        lines(density(vaf$vaf), col = "steelblue", lwd = 2)
        abline(v = c(vaf$expected_vaf, 1-vaf$expected_vaf), col = "darkorange", lty = 2, lwd = 2)
    })
    
    output$logRPlot <- renderPlot({
        
        # generate bins based on input$bins from ui.R
        logr <- simulated_logr(input$samplesize,
                               input$depth,
                               input$ta + input$tb,
                               input$ha + input$hb,
                               input$purity,
                               input$ploidy)
        bins <- seq(min(-2, min(logr$logr)), max(4, max(logr$logr)), length.out = input$bins + 1)
        
        # draw the histogram with the specified number of bins
        hist(logr$logr, breaks = bins, col = 'firebrick1', border = 'firebrick1', probability = TRUE,
             main = "Histogram of logR", xlab = "logR")
        lines(density(logr$logr), col = "firebrick", lwd = 2)
        abline(v = logr$expected_logr, col = "green3", lty = 2, lwd = 2)
    })
})
