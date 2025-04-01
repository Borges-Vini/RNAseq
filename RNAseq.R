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

system(paste("mkdir -p", file.path(parfile[2], "Reference_files")))
system(paste("mkdir -p", file.path(parfile[2], "Raw_fastq")))
system(paste("mkdir -p", file.path(parfile[2], "QC_results")))
system(paste("mkdir -p", file.path(parfile[2], "Trimmed")))
system(paste("mkdir -p", file.path(parfile[2], "Trimmed", "Unpaired")))
system(paste("mkdir -p", file.path(parfile[2], "Aligned")))
system(paste("cp -f -r", file.path(parfile[1], "**/*.fastq*"), file.path(parfile[2], "Raw_fastq")))

##### Reference_files #####

#TruSeq
truseq_content <- c(
  ">PrefixPE/1",
  "TACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  ">PrefixPE/2",
  "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
  ">PE1",
  "TACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  ">PE1_rc",
  "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
  ">PE2",
  "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
  ">PE2_rc",
  "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
)
writeLines(truseq_content, file.path(parfile[2], "Reference_files", "TruSeq3-PE-2.fa"))

#HISAT_index




##### Conda Env #####

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/bin", sep = ":"))
Sys.setenv(RETICULATE_MINICONDA_PATH = "~/anaconda3")
reticulate::use_condaenv("base", required = TRUE)
conda_create(envname = "RNAseq_env", packages = "python=3.9", channel = "conda-forge")
conda_install(envname = "RNAseq_env", packages = "fastp", channel = "bioconda")
conda_install(envname = "RNAseq_env", packages = "hisat2", channel = "bioconda")
conda_install(envname = "RNAseq_env", packages = "samtools", channel = "bioconda")
use_condaenv("RNAseq_env", required = TRUE)
system("conda run -n RNAseq_env samtools --version")

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

##### Trimming reads #####

for (sample in names(fastq_pairs)) {
  if (!is.null(fastq_pairs[[sample]]$R1) & !is.null(fastq_pairs[[sample]]$R2)) {
    message("Sample ", sample, " has both R1 and R2. Trimming adaptors now.")
    
    system(paste("conda run -n RNAseq_env fastp",
                 "--in1", file.path(fastq_pairs[[sample]]$R1), 
                 "--in2", file.path(fastq_pairs[[sample]]$R2), 
                 "--out1", file.path(parfile[2], "Trimmed", paste0(sample, "_R1.fastq.gz")),
                 "--out2", file.path(parfile[2], "Trimmed", paste0(sample, "_R2.fastq.gz")),
                 "--unpaired1", file.path(parfile[2], "Trimmed", "Unpaired", paste0(sample, "_R1.fastq.gz")),
                 "--unpaired2", file.path(parfile[2], "Trimmed", "Unpaired", paste0(sample, "_R2.fastq.gz")),
                 "--failed_out", file.path(parfile[2], "Trimmed", paste0(sample, "_failed.out")), 
     #            "--phred64",
                 "--adapter_fasta", file.path(parfile[2], "Reference_files", "TruSeq3-PE-2.fa"),
                 "--cut_front --cut_tail",
                 "--cut_window_size 4 --cut_mean_quality 20",
                 "--trim_front1 12 --trim_front2 12",
    #            "--max_len1 37 --max_len2 37",
                 "--length_required 25",
                 "--thread 12",
                 "--html ", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.html")), 
                 "--json", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.json"))
                 ))
  
    
  } else if (!is.null(fastq_pairs[[sample]]$R1)) {
    message("Sample ", sample, " has only R1. Trimming adaptors now.")
  
    system(paste("conda run -n RNAseq_env fastp",
                 "--in1", file.path(fastq_pairs[[sample]]$R1), 
                 "--out1", file.path(parfile[2], "Trimmed", paste0(sample, "_R1.fastq.gz")),
                 "--unpaired1", file.path(parfile[2], "Trimmed", "Unpaired", paste0(sample, "_R1.fastq.gz")),
                 "--failed_out", file.path(parfile[2], "Trimmed", paste0(sample, "_failed.out")), 
    #            "--phred64",
                 "--adapter_fasta", file.path(parfile[2], "Reference_files", "TruSeq3-PE-2.fa"),
                 "--cut_front --cut_tail",
                 "--cut_window_size 4 --cut_mean_quality 20",
                 "--trim_front1 12 --trim_front2 12",
    #            "--max_len1 37 --max_len2 37",
                 "--length_required 25",
                 "--thread 12",
                 "--html ", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.html")), 
                 "--json", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.json"))
                 ))
    
  } else {message("Sample ", sample, " does not have R1 or R2 (Unexpected case).")
  }}


##### Aligning reads #####  
  
    system(paste("conda run -n RNAseq_env hisat2",
                  "-1", file.path(parfile[2], "Trimmed", paste0(sample, "_R1.fastq.gz")),
                  "-2", file.path(parfile[2], "Trimmed", paste0(sample, "_R2.fastq.gz")),
                  "-x", 
                  "-S",  
   
                 
                 
                  ))
  
  
  
  center="Marshall_Genomics_Core"
  platform="Illumina"
  model="NextSeq2000"
  genome_file="/Genomes/HISAT_indexes/genome"
  
  samples=( `ls ${trimmed_dir}/*_R1.fastq.gz | \
           (while read f; do basename ${f} _R1.fastq.gz ; done) ` )
  
  if [ -z ${maxAligners} ] ; then
  maxAligners=12
  fi
  
  numSamples=${#samples[@]}
    
    start=0
    while [ $start -lt $numSamples ] ; do
    for s in ${samples[@]:start:maxAligners}; do (
      f1=${trimmed_dir}/${s}_R1.fastq.gz
      f2=${trimmed_dir}/${s}_R2.fastq.gz
      ${hisat2} --phred33 --no-mixed \
      --rg-id ${s} \
      --rg SM:${s} \
      --rg CN:${center} \
      --rg PL:${platform} \
      --rg PM:${model}  \
      -x ${genome_file} \
      -1 ${f1} \
      -2 ${f2} \
      2>${log_dir}/align_${s}.log \
      | ${samtools} sort -n -m 6G - 2>${log_dir}/name_sort_${s}.log \
      | samtools fixmate -m - -  2>${log_dir}/fixmate_${s}.log \
      | samtools sort -m 6G - 2>${log_dir}/pos_sort_${s}.log \
      | samtools markdup -r - ${aligned_dir}/${s}.bam 2>${log_dir}/markdup_${s}.log && \
      samtools index ${aligned_dir}/${s}.bam 2>${log_dir}/index_${s}.log
    )&
      done
    wait
    
    start=$[ $start + $maxAligners ]
    done
    wait
    

