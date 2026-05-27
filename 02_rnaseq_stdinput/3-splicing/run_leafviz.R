# ==============================================================================
# Description:  A utility script to launch the LeafViz interactive Shiny application 
#               for exploring the differential splicing results locally.
# Input:        ./results/*.Rdata (Compiled LeafViz objects from prepare_results)
# Output:       Interactive web application deployed to local browser.
# ==============================================================================

library(dplyr)
library(leafviz)


# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
# Set the default browser for launching the local Shiny server.
# External users: Change this to "chrome", "safari", etc., as appropriate for 
# your local machine environment.

options(browser = "firefox")

# ------------------------------------------------------------------------------
# 2. Launch Visualization
# ------------------------------------------------------------------------------
# Launch the interactive Shiny app for specific pairwise comparisons.
# Note: These commands are blocking (they start a local server). Run them one at 
# a time in your R console, interact with the UI, and press Esc to stop the server 
# before launching the next one.

# View Alternative Splicing: TLR7 vs Baseline
leafviz("./results/TLR7.Rdata")

# View Alternative Splicing: BCR vs Baseline
leafviz("./results/BCR.Rdata")

# View Alternative Splicing: DN2 vs Baseline
leafviz("./results/unstday0vs.DN2.Rdata")
