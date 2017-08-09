setupProject <- function(projectName, wd){
  
  dir.test1 <- dir.exists(paste(wd, "//", projectName,sep = ""))
  project.folder <- gsub("\\\\", "/", file.path(wd, projectName))
  
  if(!dir.test1){
    
    selected.folder <- choose.dir()
    dir.test2 <- dir.exists(paste(selected.folder, "//", projectName,sep = ""))
    project.folder <- gsub("\\\\", "/", file.path(selected.folder, projectName))
  
    if(!dir.test2){
    
      data.files <- choose.files(default = "", caption = "Select *All* Experimental Data",
                                  multi = TRUE, filters = "txt",
                                  index = nrow(Filters))
      
      bed.file <- choose.files(default = "", caption = "Select bed file",
                               multi = FALSE, filters = "bed",
                               index = nrow(Filters))
      

      dir.create(project.folder, showWarnings = FALSE)
      dir.create(file.path(project.folder,'bed'), showWarnings = FALSE)
      dir.create(file.path(project.folder,'save'), showWarnings = FALSE)
      dir.create(file.path(project.folder,'data'), showWarnings = FALSE)
      dir.create(file.path(project.folder,'output'), showWarnings = FALSE)
      
      file.copy(from=data.files, to=file.path(project.folder,'data'), 
                overwrite = TRUE, recursive = FALSE, 
                copy.mode = TRUE)
      
      file.copy(from=bed.file, to=file.path(project.folder,'bed'), 
                overwrite = TRUE, recursive = FALSE, 
                copy.mode = TRUE)
    }
  }
    
  return(project.folder)

}