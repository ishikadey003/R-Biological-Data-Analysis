# ============================================================
# EGFR SEQUENCE ANALYSIS
# ============================================================
# Gene: EGFR
# NCBI Accession: NM_005228
# Organism: Homo sapiens
# ============================================================

# ============================================================
# Q1. DOWNLOAD AND READ THE EGFR SEQUENCE
# ============================================================

egfr_fasta <- entrez_fetch(
  db = "nuccore",
  id = "NM_005228",
  rettype = "fasta",
  retmode = "text"
)

writeLines(egfr_fasta, "egfr.fasta")

egfr_data <- readDNAStringSet("egfr.fasta")

egfr_seq <- egfr_data[[1]]


# Determine sequence length
egfr_length <- width(egfr_data)

egfr_length

# ============================================================
# Q2. WHAT ARE THE FIRST 30 BASES?
# ============================================================

egfr_first30 <- subseq(
  egfr_seq,
  start = 1,
  end = 30
)

as.character(egfr_first30)

# ============================================================
# Q3. HOW MANY TIMES DOES "ATG" APPEAR?
# ============================================================

egfr_atg_count <- countPattern(
  "ATG",
  egfr_seq
)

egfr_atg_count

# ============================================================
# Q4. REVERSE COMPLEMENT OF THE FIRST 30 BASES
# ============================================================

egfr_reverse_complement <- reverseComplement(egfr_first30)

as.character(egfr_reverse_complement)

# ============================================================
# Q5. CREATE A GRANGES OBJECT FOR EGFR
# ============================================================

egfr_gr <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(
    start = 55019017,
    end = 55211628
  ),
  strand = "+"
)

# Width requested in the assignment: end - start
egfr_end_minus_start <- end(egfr_gr) - start(egfr_gr)

egfr_end_minus_start

# Actual width of the GRanges interval
egfr_gr_width <- width(egfr_gr)

egfr_gr_width

# ============================================================
# Q6. EGFR ANNOTATION
# ============================================================

egfr_annotation <- select(
  org.Hs.eg.db,
  keys = "EGFR",
  keytype = "SYMBOL",
  columns = c("ENTREZID", "GENENAME")
)

egfr_annotation

egfr_entrez_id <- egfr_annotation$ENTREZID
egfr_gene_name <- egfr_annotation$GENENAME

egfr_entrez_id
egfr_gene_name
egfr_gr

# ============================================================
# Q7. SUMMARY
# ============================================================

# EGFR is a 9,905-base nucleotide sequence in the NCBI record
# analyzed, with 160 occurrences of the ATG pattern. The EGFR
# gene is located on chromosome 7 (chr7), spanning positions
# 55,019,017–55,211,628 on the positive strand, and its Entrez
# ID is 1956.