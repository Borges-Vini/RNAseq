##### ABOUT #####
# RStudio v.2024.12.0
# R v. 4.4.1

##### Future Changes #####

#change the parfile input
#Create warning and error message for each step


##### Packages #####

install.packages(c("data.table", "devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "reticulate", "fastqcr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "Biobase", "tximport", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79", 
                       "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))


packages <- c("data.table", "devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "reticulate", "fastqcr")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}))
fastqc_install()

#fastp - Checking and installing
  Sys.getenv("PATH")
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/bin", sep = ":"))
  Sys.setenv(RETICULATE_MINICONDA_PATH = "~/anaconda3")
  reticulate::use_condaenv("base", required = TRUE)
  conda_create(envname = "fastp_env", packages = "python=3.9", channel = "conda-forge")
  conda_install(envname = "fastp_env", packages = "fastp", channel = "bioconda")
  use_condaenv("fastp_env", required = TRUE)
  system("conda run -n fastp_env fastp --version")

  
##### Variables #####

args <- commandArgs(trailingOnly=TRUE)
parfilepath <- args[1]
parfilepath <- "/Users/vborges/RNAseq/RNAseq.par"
if (!file.exists(parfilepath)) {stop("Error: The file does not exist at the specified path: ", parfilepath)
  } else {print("Reading parameter file")}
parfile <- Filter(function(x) !grepl("^#", x) && nzchar(x), trimws(readLines(parfilepath, warn = FALSE)))

##### Output Directory #####

if (!file.exists(parfile[2])) {system(paste("mkdir -p", parfile[2]))
} else {
  system(paste("rm -rf", parfile[2])) 
  system(paste("mkdir -p", parfile[2])) 
  print("Folder already exists. Overwriting.")
}

system(paste("mkdir -p", file.path(parfile[2], "Raw_fastq")))
system(paste("mkdir -p", file.path(parfile[2], "QC_results")))
system(paste("mkdir -p", file.path(parfile[2], "Trimmed")))
system(paste("cp -f -r", file.path(parfile[1], "**/*.fastq*"), file.path(parfile[2], "Raw_fastq")))


##### FastQC #####

java_path <- dirname(dirname(parfile[4]))
Sys.setenv(JAVA_HOME=java_path)
system(paste("java -version"))
fastqc(fq.dir = file.path(parfile[2], "Raw_fastq"), qc.dir = file.path(parfile[2], "QC_results"), threads = 12)

##### Summarizing QC Results #####

if (parfile[3] %in% c("YES", "yes", "Yes")) { 
  
#multiqc - Checking and installing
  venv_path <- "multiqc_env"
  if (!file.exists(venv_path)) {
    virtualenv_create(venv_path)}
  virtualenv_install(venv_path, "multiqc")
  use_virtualenv(venv_path, required = TRUE)
  multiqc_path <- file.path(path.expand("~"), ".virtualenvs", venv_path, "bin", "multiqc")
  system(paste(multiqc_path, "--version"))

#fastq - Checking and Running
  if (!file.exists(file.path(parfile[2], "Raw_fastq"))) { stop("Error: The folder does not exist at the specified path: ", file.path(parfile[2], "Raw_fastq")) }
  fastq_files <- list.files(file.path(parfile[2], "Raw_fastq"), pattern = "\\.fastq(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(fastq_files) == 0) {stop("Error: No .fastq files found in the folder: ", file.path(parfile[2], "Raw_fastq"))
  } else { print("Valid folder with .fastq files found. Proceeding with analysis.")
  system(paste(multiqc_path, file.path(parfile[2], "QC_results"),"-o", file.path(parfile[2], "QC_results")))}
}

##### Identifying Single/Paired-End #####

fastq_files <- list.files(file.path(parfile[2], "Raw_fastq"), pattern = "\\.fastq(\\.gz)?$", full.names = TRUE, recursive = TRUE)
fastq_pairs <- list()
for (file in fastq_files) {
  sample_name <- sub("_R[12].fastq(\\.gz)?$", "", basename(file))
  if (!sample_name %in% names(fastq_pairs)) {
    fastq_pairs[[sample_name]] <- list(R1 = NULL, R2 = NULL)
  }
  if (grepl("R1", file)) {
    fastq_pairs[[sample_name]]$R1 <- file
  } else if (grepl("R2", file)) {
    fastq_pairs[[sample_name]]$R2 <- file
  }
}

##### Trimming #####

for (sample in names(fastq_pairs)) {
  if (!is.null(fastq_pairs[[sample]]$R1) & !is.null(fastq_pairs[[sample]]$R2)) {
    message("Sample ", sample, " has both R1 and R2. Trimming adaptors now.")
    
    trim_fastq(fastq1 = file.path(fastq_pairs[[sample]]$R1), fastq2 = file.path(fastq_pairs[[sample]]$R2), adapter1 = NULL, adapter2 = NULL,
               illumina = FALSE, nextera = FALSE, small_rna = FALSE,
               minlength = 20, minqual = 20, trimN = TRUE,
               retainUnpaired = TRUE, retain1length = 35, retain2length = 35,
               clipR1 = NULL, clipR2 = NULL, clip3primeR1 = NULL,
               clip3primeR2 = NULL, robust_check = FALSE, dest.dir = NULL,
               threads = NULL, trimgalore = "trim_galore")
    
  } else if (!is.null(fastq_pairs[[sample]]$R1)) {message("Sample ", sample, " has only R1. Trimming adaptors now.")
  
    trim_fastq
    
    } else {message("Sample ", sample, " does not have R1 or R2 (Unexpected case).")
  }}


  



