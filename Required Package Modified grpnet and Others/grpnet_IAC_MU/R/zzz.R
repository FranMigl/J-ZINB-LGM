# Some code (for startup message) is re-purposed from the R package
# "mclust" (Scrucca et al., 2016) https://cran.r-project.org/package=mclust

grpnetStartupMessage <- 
  function(){
    msg <- c(paste0("                                    _   
         __ _ _ __ _ __  _ __   ___| |_ 
        / _` | '__| '_ \\| '_ \\ / _ \\ __|
       | (_| | |  | |_) | | | |  __/ |_ 
        \\__, |_|  | .__/|_| |_|\\___|\\__|
        |___/     |_|      version ", 
                    packageVersion("grpnet"), "\n"),
             "\nType 'citation(\"grpnet\")' to cite this package.\n")
    return(msg)
  }

.onAttach <- 
  function(lib, pkg){
    msg <- grpnetStartupMessage()
    if(!interactive()) msg[1] <- paste("Package 'grpnet' version", packageVersion("grpnet"))
    packageStartupMessage(msg)      
    invisible()
  }




