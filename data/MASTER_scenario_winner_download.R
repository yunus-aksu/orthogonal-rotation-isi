# Install the googledrive package if it is not already installed
# install.packages("googledrive")
library(googledrive)

# Download the file from Google Drive to your local data folder (1.3 GB)
drive_download(
  "https://drive.google.com/file/d/1uIPgbcNrZ--cJhbrf1_7WnY6wnA54G44/view?usp=sharing", 
  path = "data/MASTER_scenario_winner.rds", 
  overwrite = TRUE
)

# Load the file into R
model_results <- readRDS("data/MASTER_scenario_winner.rds")
