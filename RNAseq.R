# About -------------------------------------------------------------------
# RStudio v.2024.12.0
# R v. 4.4.1

# Future Changes ----------------------------------------------------------

#change the parfile input
#Create warning and error message for each step

# Packages ----------------------------------------------------------------

#General
required_packages <- c("data.table", "devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "reticulate", "fastqcr")
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)
if (length(missing_packages) > 0) {install.packages(missing_packages)}
invisible(lapply(required_packages, library, character.only = TRUE))
fastqc_install()

#BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
bioc_packages <- c("GEOquery", "Biobase", "tximport", "EnsDb.Hsapiens.v86","EnsDb.Mmusculus.v79", "AnnotationDbi","org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db")
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {BiocManager::install(missing_bioc)}
invisible(lapply(bioc_packages, library, character.only = TRUE))

# PAR file ----------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
parfilepath <- args[1]
parfilepath <- "/Users/vborges/RNAseq/RNAseq.par"
if (!file.exists(parfilepath)) {stop("Error: The file does not exist at the specified path: ", parfilepath)
  } else {print("Reading parameter file")}
parfile <- Filter(function(x) !grepl("^#", x) && nzchar(x), trimws(readLines(parfilepath, warn = FALSE)))

# Output Directory --------------------------------------------------------

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
system(paste("mkdir -p", file.path(parfile[2], "Aligned", "Logs")))
system(paste("cp -f -r", file.path(parfile[1], "**/*.fastq*"), file.path(parfile[2], "Raw_fastq")))

# Reference_files ---------------------------------------------------------

if (grepl("(?i)truseq3-pe", parfile[5], perl = TRUE)) {
  
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
  adapter_path <- file.path(parfile[2], "Reference_files", "TruSeq3-PE-2.fa")
  
  }else {
    adapter_path <- parfile[5]
    }

if (grepl("(?i)grcm38", parfile[6], perl = TRUE)) {
  options(timeout = 600)
  download.file("https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download", 
              destfile = file.path(parfile[2], "Reference_files", "grcm38.tar.gz"), mode = "wb")
  hisatindexgenome_path <- file.path(parfile[2], "Reference_files", "grcm38.tar.gz")
}else {
  hisatindexgenome_path <- parfile[6]
}

# Conda Env ---------------------------------------------------------------

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/anaconda3/bin", sep = ":"))
Sys.setenv(RETICULATE_MINICONDA_PATH = "~/anaconda3")
reticulate::use_condaenv("base", required = TRUE)
conda_create(envname = "RNAseq_env",packages = c("python=3.10", "samtools=1.21", "fastp","hisat2"), channel = c("conda-forge", "bioconda"))
use_condaenv("RNAseq_env", required = TRUE)
system("conda run -n RNAseq_env samtools --version")

# FastQC ------------------------------------------------------------------
java_path <- dirname(dirname(parfile[4]))
Sys.setenv(JAVA_HOME=java_path)
system(paste("java -version"))
fastqc(fq.dir = file.path(parfile[2], "Raw_fastq"), qc.dir = file.path(parfile[2], "QC_results"), threads = as.integer(parallel::detectCores() * 2 / 3))

# Summarizing QC Results --------------------------------------------------

if (grepl("(?i)y", parfile[3], perl = TRUE)) {
    
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

# Identifying Single/Paired-End -------------------------------------------
fastq_files <- list.files(file.path(parfile[2], "Raw_fastq"), pattern = "\\.fastq(\\.gz)?$", full.names = TRUE, recursive = TRUE)
samples_fastq_files <- unique(sub("_R[12]\\.fastq\\.gz$", "", basename(fastq_files)))

fastq_pairs <- list()
for (sample in samples_fastq_files) {
  fastq_pairs[[sample]] <- list(
    R1 = file.path(file.path(parfile[2], "Raw_fastq"), paste0(sample, "_R1.fastq.gz")),
    R2 = file.path(file.path(parfile[2], "Raw_fastq"), paste0(sample, "_R2.fastq.gz"))
  )
}
# Trimming reads ----------------------------------------------------------
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
                 "--adapter_fasta", adapter_path,
                 "--cut_front --cut_tail",
                 "--cut_window_size 4 --cut_mean_quality 15",
                 "--trim_front1 12 --trim_front2 12",
                 "--max_len1 37 --max_len2 37",
                 "--length_required 25",
                 "--thread", as.integer(parallel::detectCores() * 2 / 3),
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
                 "--adapter_fasta", adapter_path,
                 "--cut_front --cut_tail",
                 "--cut_window_size 4 --cut_mean_quality 15",
                 "--trim_front1 12 --trim_front2 12",
                 "--max_len1 37 --max_len2 37",
                 "--length_required 25",
                 "--thread", as.integer(parallel::detectCores() * 2 / 3),
                 "--html ", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.html")), 
                 "--json", file.path(parfile[2], "Trimmed", paste0(sample, "_fastp_report.json"))
                 ))
    
  } else {message("Sample ", sample, " does not have R1 or R2 (Unexpected case).")
  }}


# Identifying Trimmed -----------------------------------------------------

trimmed_files <- list.files(file.path(parfile[2], "Trimmed"), pattern = "_R[12]\\.fastq\\.gz$", full.names = TRUE)
samples_trimmed_files <- unique(sub("_R[12]\\.fastq\\.gz$", "", basename(trimmed_files)))

trimmed_pairs <- list()
for (sample in samples_trimmed_files) {
  trimmed_pairs[[sample]] <- list(
    R1 = file.path(file.path(parfile[2], "Trimmed"), paste0(sample, "_R1.fastq.gz")),
    R2 = file.path(file.path(parfile[2], "Trimmed"), paste0(sample, "_R2.fastq.gz"))
  )
}
# Aligning reads ----------------------------------------------------------

for (sample in names(trimmed_pairs)[27:length(trimmed_pairs)]) {
  
  if (!is.null(trimmed_pairs[[sample]]$R1) & !is.null(trimmed_pairs[[sample]]$R2)) {
    message("Aligning reads now:", sample)
    
    maxAligners <- as.integer(parallel::detectCores() * 2 / 3)
    
      pair <- trimmed_pairs[[sample]]
      
      if (!is.null(pair$R1) && !is.null(pair$R2)) {
        
        if (!file.exists(pair$R1) || !file.exists(pair$R2)) {
          message("Missing FASTQ for sample: ", sample)
          next
        }

        cmd <- paste(
          "conda run -n RNAseq_env bash -c", shQuote(paste(
            "set -euo pipefail;",
            
            # Step 1: Align to SAM
            "hisat2 --phred33 --no-mixed",
            "-p", maxAligners,
            "--rg-id", sample,
            "--rg", paste0("SM:", sample),
            "--rg", "CN:Marshall_Genomics_Core",
            "--rg", "PL:Illumina",
            "--rg", "PM:NextSeq2000",
            "-x", shQuote(hisatindexgenome_path),
            "-1", shQuote(pair$R1),
            "-2", shQuote(pair$R2),
            "-S", shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_align.log"))), "&&",
            
            # Step 2: Sort by name (required for fixmate)
            "samtools sort -n -@", maxAligners,
            "-o", shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_name_sort.log"))), "&&",
            
            # Step 3: Fixmate (add ms and MC tags needed by markdup)
            "samtools fixmate -m",
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_fixmate.log"))), "&&",
            
            # Step 4: Coordinate sort (required for markdup)
            "samtools sort -@", maxAligners,
            "-o", shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))),
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_pos_sort.log"))), "&&",
            
            # Step 4.5: quickcheck to validate sorted BAM
            "samtools quickcheck -v",
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))), "&&",
            
            # Step 5: Mark duplicates (write to final dedup BAM)
            "samtools markdup -r",
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))),
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".bam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_markdup.log"))), "&&",
            
            # Step 6: Index the final deduplicated BAM
            "samtools index",
            shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".bam"))),
            "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_index.log")))
          ))
        )
        
        # Run command
        system(cmd)
        system(paste("rm", 
               paste0(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
               paste0(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
               paste0(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
               paste0(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam")))))
               
      }
  } else if (!is.null(trimmed_pairs[[sample]]$R1)) {
    message("Aligning reads now:", sample)
    
    maxAligners <- as.integer(parallel::detectCores() * 2 / 3)

    pair <- trimmed_pairs[[sample]]
    
    if (!is.null(pair$R1)) {
      
      if (!file.exists(pair$R1)) {
        message("Missing FASTQ for sample: ", sample)
        next
      }
      
      cmd <- paste(
        "conda run -n RNAseq_env bash -c", shQuote(paste(
          "set -euo pipefail;",
          
          # Step 1: Align to SAM
          "hisat2 --phred33 --no-mixed",
          "-p", maxAligners,
          "--rg-id", sample,
          "--rg", paste0("SM:", sample),
          "--rg", "CN:Marshall_Genomics_Core",
          "--rg", "PL:Illumina",
          "--rg", "PM:NextSeq2000",
          "-x", shQuote(hisatindexgenome_path),
          "-1", shQuote(pair$R1),
          "-2", shQuote(pair$R2),
          "-S", shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_align.log"))), "&&",
          
          # Step 2: Sort by name (required for fixmate)
          "samtools sort -n -@", maxAligners,
          "-o", shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_name_sort.log"))), "&&",
          
          # Step 3: Fixmate (add ms and MC tags needed by markdup)
          "samtools fixmate -m",
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_fixmate.log"))), "&&",
          
          # Step 4: Coordinate sort (required for markdup)
          "samtools sort -@", maxAligners,
          "-o", shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))),
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_pos_sort.log"))), "&&",
          
          # Step 4.5: quickcheck to validate sorted BAM
          "samtools quickcheck -v",
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))), "&&",
          
          # Step 5: Mark duplicates (write to final dedup BAM)
          "samtools markdup -r",
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam"))),
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".bam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_markdup.log"))), "&&",
          
          # Step 6: Index the final deduplicated BAM
          "samtools index",
          shQuote(file.path(parfile[2], "Aligned", paste0(sample, ".bam"))),
          "2>", shQuote(file.path(parfile[2], "Aligned", "Logs", paste0(sample, "_index.log")))
        ))
      )
      
      # Run command
      system(cmd)
      system(paste("rm", 
                   paste0(file.path(parfile[2], "Aligned", paste0(sample, ".sam"))),
                   paste0(file.path(parfile[2], "Aligned", paste0(sample, "_name_sorted.bam"))),
                   paste0(file.path(parfile[2], "Aligned", paste0(sample, "_fixmate.bam"))),
                   paste0(file.path(parfile[2], "Aligned", paste0(sample, "_coord_sorted.bam")))))
      
    }
  } else {message("Sample ", sample, " does not have R1 or R2 (Unexpected case).")}
  
  }

    


# Alignment rate -------------------------------------------------------------------

log_aligned_bam <- list.files(file.path(parfile[2], "Aligned", "Logs"), pattern = "_merged_align.log$", full.names = TRUE)
samples <- character(length(log_aligned_bam))
rates <- numeric(length(log_aligned_bam))

for (i in seq_along(log_aligned_bam)) {
  lines <- readLines(log_aligned_bam[i], warn = FALSE)
  rates[i] <- if (length(lines) >= 9) {
    as.numeric(sub("([0-9.]+)%.*", "\\1", lines[9]))
  } else {
    NA
  }
  samples[i] <- sub(".*/([^/]+)_merged_align\\.log$", "\\1", log_aligned_bam[i])
}

overall_align_rate <- data.frame(sample = samples, rate = rates)
overall_align_rate <- overall_align_rate[order(overall_align_rate$rate), ]

mean_rate <- mean(overall_align_rate$rate, na.rm = TRUE)
print(mean_rate)
