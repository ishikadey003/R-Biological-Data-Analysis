# ============================================================
# BRCA1 SEQUENCE ANALYSIS
# ============================================================
# Gene: BRCA1
# NCBI Accession: NM_007294
# Organism: Homo sapiens

# Load required packages
library(rentrez)
library(Biostrings)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)

# This script answers Q1-Q7 of the assignment.
# ============================================================
# ------------------------------------------------------------
# DOWNLOAD BRCA1 SEQUENCE FROM NCBI
# ------------------------------------------------------------

brca1_fasta <- entrez_fetch(
  db = "nuccore",
  id = "NM_007294",
  rettype = "fasta",
  retmode = "text"
)

# Save the FASTA file
writeLines(brca1_fasta, "brca1.fasta")
# ------------------------------------------------------------
# READ FASTA SEQUENCE
# ------------------------------------------------------------

brca1_data <- readDNAStringSet("brca1.fasta")

# Extract the individual DNA sequence
brca1_seq <- brca1_data[[1]]
# ============================================================
# Q1. HOW MANY BASES LONG IS THE SEQUENCE?
# ============================================================

brca1_length <- width(brca1_data)

brca1_length
7088

# ============================================================
# Q2. WHAT ARE THE FIRST 30 BASES?
# ============================================================

brca1_first30 <- subseq(
  brca1_seq,
  start = 1,
  end = 30
)

brca1_first30
as.character(brca1_first30)
"GCTGAGACTTCCTGGACGGGGGACAGGCTG"

# ============================================================
# Q3. HOW MANY TIMES DOES "ATG" APPEAR IN THE SEQUENCE?
# ============================================================

brca1_atg_count <- countPattern(
  "ATG",
  brca1_seq
)

brca1_atg_count
140

# ============================================================
# Q4. REVERSE COMPLEMENT OF THE FIRST 30 BASES
# ============================================================

brca1_reverse_complement <- reverseComplement(brca1_first30)

brca1_reverse_complement
"CAGCCTGTCCCCCGTCCAGGAAGTCTCAGC"

# ============================================================
# Q5. CREATE A GRANGES OBJECT FOR BRCA1
# ============================================================

brca1_gr <- GRanges(
  seqnames = "chr17",
  ranges = IRanges(
    start = 43044295,
    end = 43125483
  ),
  strand = "-"
)

brca1_gr

# Calculate genomic range width
brca1_gr_width <- width(brca1_gr)
brca1_gr_width

# Calculate end - start
brca1_end_minus_start <- end(brca1_gr) - start(brca1_gr)
brca1_end_minus_start

# ============================================================
# Q6. BRCA1 ANNOTATION
# ============================================================

brca1_annotation <- select(
  org.Hs.eg.db,
  keys = "BRCA1",
  keytype = "SYMBOL",
  columns = c("ENTREZID", "GENENAME")
)

brca1_annotation

brca1_entrez_id <- brca1_annotation$ENTREZID[1]
brca1_gene_name <- brca1_annotation$GENENAME[1]

brca1_entrez_id
brca1_gene_name

# ============================================================
# Q7. SUMMARY
# ============================================================

# BRCA1 is a 7,088-base nucleotide sequence in the NCBI record
# analyzed, with 140 occurrences of the ATG pattern. The gene
# is located on chromosome 17 (chr17), spanning positions
# 43,044,295–43,125,483 on the negative strand.


