##### ABOUT #####

# RStudio v.2024.12.0
# R v. 4.4.1

##### Packages #####

install.packages(c("data.table", "devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "reticulate", "fastqcr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "Biobase", "tximport", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79", 
                       "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))

devtools::install_github("uc-bd2k/GREP2")
library(GREP2)

packages <- c("data.table", "devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "reticulate", "fastqcr")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}))
fastqc_install()

##### Variables #####

args <- commandArgs(trailingOnly=TRUE)
parfilepath <- args[1]
parfilepath <- "/Users/vborges/RNAseq/RNAseq.par"
if (!file.exists(parfilepath)) {stop("Error: The file does not exist at the specified path: ", parfilepath)
  } else {print("Reading parameter file")}
parfile <- Filter(function(x) !grepl("^#", x) && nzchar(x), trimws(readLines(parfilepath, warn = FALSE)))

##### FastQC #####

java_path <- dirname(dirname(parfile[3]))
Sys.setenv(JAVA_HOME=java_path)
system("java -version")
fastqc(fq.dir = (parfile[1]))

##### Summarizing QC Results #####

if (parfile[2] %in% c("YES", "yes", "Yes")) { 
  
#multiqc - Checking and installing
  venv_path <- "multiqc_env"
  if (!file.exists(venv_path)) {
    virtualenv_create(venv_path)}
  virtualenv_install(venv_path, "multiqc")
  use_virtualenv(venv_path, required = TRUE)
  multiqc_path <- file.path(path.expand("~"), ".virtualenvs", venv_path, "bin", "multiqc")
  system(paste(multiqc_path, "--version"))

#fastq - Checking and Running
  if (!file.exists(parfile[1])) { stop("Error: The folder does not exist at the specified path: ", parfile[1]) }
  fastq_files <- list.files(parfile[1], pattern = "\\.fastq(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(fastq_files) == 0) {stop("Error: No .fastq files found in the folder: ", parfile[1])
  } else { print("Valid folder with .fastq files found. Proceeding with analysis.")
  system(paste(multiqc_path, parfile[1],"-o", parfile[1]))}
}


