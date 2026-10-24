#' Calculate Peak RAM Used
#'
#' This function monitors the total and peak RAM used by any number of
#'  R expressions or functions.
#'
#' When working with big datasets, RAM conservation is critically
#'  important. However, it is not always enough to just monitor the
#'  size of the objects created. So-called "copy-on-modify" behavior,
#'  characteristic of R, means that some expressions or functions may
#'  require an unexpectedly large amount of RAM overhead. For example,
#'  replacing a single value in a matrix (e.g., with \code{'[<-'})
#'  duplicates that matrix in the backend, making this task
#'  require twice as much RAM as that used by the matrix itself.
#'  The \code{peakRAM} package makes it easy to monitor the total
#'  and peak RAM used so that developers can quickly identify and
#'  eliminate RAM hungry code.
#'
#' @param ... R expressions or function calls. Anonymous functions
#'  (e.g., \code{function() 1:1e7}) also accepted.
#' @return A \code{data.frame} tallying total and peak RAM use.
#'
#' @examples
#' peakRAM(function() 1:1e7,
#'         1:1e7,
#'         1:1e7 + 1:1e7,
#'         1:1e7 * 2)
#' @export
peakRAM <- function(...){
  
  # Capture R expressions or function calls
  args <- c(as.list(match.call(expand.dots = FALSE)$`...`), NULL)
  Function <- sapply(args, function(e) paste0(deparse(e), collapse = ""))
  Function <- sapply(Function, function(e) gsub(" ", "", e))
  
  # Initialize containers for output
  numfunc <- length(args)
  totaltime <- vector("numeric", numfunc)
  RAMused <- vector("numeric", numfunc)
  RAMpeak <- vector("numeric", numfunc)
  
  i <- 1
  for(arg in args){
    
    # Reset garbage collector and save baseline
    start <- gc(verbose = FALSE, reset = TRUE)
    
    # Evaluate regular and anonymous functions
    evalTime <- system.time(result <- eval.parent(arg))
    if(inherits(result, "function")){
      
      evalTime <- system.time(output <- result())
      rm(result)
      
    }else{
      
      output <- result
      rm(result)
    }
    
    # Call garbage collector and save post-eval
    end <- gc(verbose = FALSE, reset = FALSE)
    rm(output)
    
    # Calculate total and peak RAM used
    totaltime[i] <- as.numeric(evalTime["elapsed"])
    get_column_mb = function(object, before_column) {
      ind = which(colnames(object) == before_column)
      stopifnot(colnames(object)[ind+1] == "(Mb)")
      if (length(ind) == 0) {
        stop("No value found for that column")
      }
      # need to ref index and not colname because duplicated colnames for (Mb)
      object[2, ind + 1]
    }
    
    RAMused[i] <- get_column_mb(end, "used")- get_column_mb(start, "used")
    RAMpeak[i] <- get_column_mb(end, "max used") - get_column_mb(start, "max used")
    
    i <- i + 1
  }
  
  data.frame("Function_Call" = Function,
             "Elapsed_Time_sec" = totaltime,
             "Total_RAM_Used_MiB" = RAMused,
             "Peak_RAM_Used_MiB" = RAMpeak,
             row.names = 1:numfunc)
}

convert_to_bytes <- function(x) {
  
  # Remove spaces for safety
  x <- trimws(x)
  
  if (grepl("MB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("MB$", "", x, ignore.case = TRUE))
    return(value * 1024^2)  # MB → bytes
  } else if (grepl("GB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("GB$", "", x, ignore.case = TRUE))
    return(value * 1024^3)  # GB → bytes
  } else if (grepl("KB$", x, ignore.case = TRUE)) {
    value <- as.numeric(sub("KB$", "", x, ignore.case = TRUE))
    return(value * 1024)  # KB → bytes
  } 
  
  return(x)
}

transform_breaks_time <- function(x, breaks){
  breaks <- c(-1,breaks)
  x_ind <- as.integer(cut(x, breaks = breaks))
  diff <- breaks[-1] - breaks[-length(breaks)]
  low <- breaks[x_ind]
  high <- breaks[x_ind + 1]
  x_norm <- ((x - low)/(high - low)) + (x_ind - 1)
  return(x_norm)
}

transform_breaks_mem <- function(x){
  breaks <- c(1,1024, 1024^2, 1024^3, 4*1024^3)
  x_ind <- as.integer(cut(x, breaks = breaks))
  diff <- breaks[-1] - breaks[-length(breaks)]
  low <- breaks[x_ind]
  high <- breaks[x_ind + 1]
  # x_norm <- ((x - low)/(high - low)) + (x_ind - 1)
  x_norm <- (x)/(1024^(x_ind)) + (x_ind - 1)
  return(x_norm)
}