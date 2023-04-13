std_error <- function(x) { # Function for calculating standard error of the mean (SEM)
	sd(x) / sqrt(length(x)) 
}