# Load necessary library
library(ggplot2)

# Set working directory inside the tests folder
setwd("~/Documents/PACS-challenges/challenge_3/tests")

# Read the data from the CSV file
data <- read.csv("plot.csv")

# Function to generate plots for each process group
generate_plots <- function(data) {
  # Create a directory to store the plots
  dir.create("./plot", showWarnings = FALSE)
  
  # Group the data by "Processes"
  grouped_data <- split(data, data$Processes)
  
  # Loop through each group and create an XY plot
  for (process_group in names(grouped_data)) {
    print(paste("Generating plot for Process Group:", process_group))
    
    # Filter the data for the current process group
    filtered_data <- grouped_data[[process_group]]
    
    # Create the XY plot connecting the points
    p <- ggplot(filtered_data, aes(x = Size, y = Time_ms)) +
      geom_line() + 
      geom_point() +# Use points instead of lines for individual data points
      labs(title = paste("Process Group:", process_group),
           x = "Size",
           y = "Total Time (ms)") +
      scale_y_log10() +
      scale_x_continuous(trans = 'log2')
      theme_minimal()
    
    # Save the plot with with processes and size into the name
    ggsave(paste("./plot/plot_size_time_", process_group, ".png", sep = ""), plot = p)
    
    # Create the XY plot connecting the points in log scale

    p <- ggplot(filtered_data, aes(x = Size, y = Error)) +
      geom_line() + 
      geom_point() +# Use points instead of lines for individual data points
      labs(title = paste("Process Group:", process_group),
           x = "Size",
           y = "Error") +
      scale_y_log10() +
      scale_x_continuous(trans = 'log2')
      theme_minimal()
    
    # Save the plot with with processes and size into the name
    ggsave(paste("./plot/plot_size_error_", process_group, ".png", sep = ""), plot = p)
    
    # Print the plot
    print(p)
  }
}

# Call the function with the loaded data
generate_plots(data)
