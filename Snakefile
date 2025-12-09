# Snakefile – UPGRADED VERSION 8.4-lite (no AMR / no Plasmid modules)
# - Polypolish v0.5.x (no subcommand)
# - Name-sorted BAM for Polypolish + coord-sorted BAM for general use
# - BWA-MEM uses -a to keep all alignments for Polypolish
# - fastp outputs go to QC/fastp/
# - BUSCO HTML globbing in final report

import os, glob
from datetime import datetime

configfile: "config.yaml"

# --- Configuration ---
ISOLATE       = config["isolate"]
BASE          = config.get("base_dir", ".")
THREADS       = int(config.get("threads", 8))
KEEP_PCT      = int(config.get("keep_percent", 95))
RACON_ROUNDS  = int(config.get("racon_rounds", 2))
BUSCO_LINEAGE = config.get("busco_lineage", "bacteria_odb10")
MEDAKA_MODEL  = config.get("medaka_model", "")
TARGET_TAXID  = config.get("target_taxid", None)
TOTAL_STEPS   = 24  # AMR + Platon removed

# --- Paths ---
ISOD = os.path.join(BASE, ISOLATE)
DB   = os.path.join(BASE, "db")
RAW  = os.path.join(ISOD, "raw")
ILL  = os.path.join(RAW, "illumina")
ONT  = os.path.join(RAW, "nanopore")
TRIM = os.path.join(ISOD, "trimmed")
WORK = os.path.join(ISOD, "work")
ASM  = os.path.join(ISOD, "asm")
QC   = os.path.join(ISOD, "reports")
ANN  = os.path.join(ISOD, "annotation")
LOGS = os.path.join(ISOD, "logs")

KRAKEN2_DB_PATH = os.path.join(DB, "kraken2_std_db")
KRAKEN2_DIR    = os.path.join(QC, "kraken2")
FASTQC_PRE_DIR = os.path.join(QC, "fastqc_pre")
NP_RAW_DIR     = os.path.join(QC, "nanoplot_raw")
NP_FILT_DIR    = os.path.join(QC, "nanoplot_filt")
SEQKIT_DIR     = os.path.join(QC, "seqkit")
COV_DIR        = os.path.join(QC, "coverage")
QUAST_DIR      = os.path.join(QC, "quast")
BUSCO_DIR      = os.path.join(QC, "busco")
MQC_DIR        = os.path.join(QC, "multiqc")

MERGE_R1       = os.path.join(WORK, "illumina_R1.merged.fq.gz")
MERGE_R2       = os.path.join(WORK, "illumina_R2.merged.fq.gz")
CLEAN_R1       = os.path.join(WORK, "illumina_R1.clean.fq.gz")
CLEAN_R2       = os.path.join(WORK, "illumina_R2.clean.fq.gz")

ONT_MERGED     = os.path.join(WORK, "ont_merged.fastq.gz")
ONT_CLEAN      = os.path.join(WORK, "ont_clean.fastq.gz")

TRIM_R1        = os.path.join(TRIM, "R1.trimmed.fq.gz")
TRIM_R2        = os.path.join(TRIM, "R2.trimmed.fq.gz")
ONT_FILTERED   = os.path.join(WORK, "ont_filtered.fastq.gz")

UNI_DIR        = os.path.join(ASM, "unicycler")
UNI_ASM        = os.path.join(UNI_DIR, "assembly.fasta")
ALL_CONTIGS    = os.path.join(WORK, "assembly.all_contigs.fasta")
ROUND0         = os.path.join(WORK, "assembly.racon0.fasta")
RACON_LAST     = os.path.join(WORK, f"assembly.racon{RACON_ROUNDS}.fasta")

MEDAKA_DIR     = os.path.join(WORK, "medaka")
MEDAKA_OUT     = os.path.join(MEDAKA_DIR, "consensus.fasta")

NAME_BAM       = os.path.join(WORK, "illumina.namesort.bam")

POLISHED       = os.path.join(WORK, "assembly.polished.fasta")
FINAL          = os.path.join(WORK, "assembly.final.fasta")

PROKKA_DIR     = ANN
PROKKA_PREFIX  = ISOLATE
PROKKA_FAA     = os.path.join(PROKKA_DIR, f"{PROKKA_PREFIX}.faa")

# --- Helpers ---
def illumina_r1():
    pats = ["*_R1*.fastq.gz","*_R1*.fq.gz","*_R1*.fastq","*_R1*.fq"]
    return sorted(sum((glob.glob(os.path.join(ILL,p)) for p in pats), []))

def illumina_r2():
    pats = ["*_R2*.fastq.gz","*_R2*.fq.gz","*_R2*.fastq","*_R2*.fq"]
    return sorted(sum((glob.glob(os.path.join(ILL,p)) for p in pats), []))

def nanopore_files():
    pats = ["*.fastq.gz","*_fq.gz","*.fastq","*.fq"]
    return sorted(sum((glob.glob(os.path.join(ONT,p)) for p in pats), []))

def get_log_message(step_num, total_steps, message_text):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"---\n[{now}] --- Step {step_num} of {total_steps}: {message_text} ---\n"


# --- rule all ---
rule all:
    input:
        FINAL,
        os.path.join(QUAST_DIR, "report.txt"),
        os.path.join(BUSCO_DIR, f"busco_{ISOLATE}", "short_summary.txt"),
        os.path.join(MQC_DIR, "multiqc_report.html"),
        os.path.join(PROKKA_DIR, f"{PROKKA_PREFIX}.gff"),
        os.path.join(KRAKEN2_DIR, "illumina.report.txt"),
        os.path.join(KRAKEN2_DIR, "nanopore.report.txt"),
        FASTQC_PRE_DIR,
        NP_RAW_DIR,
        NP_FILT_DIR,
        os.path.join(SEQKIT_DIR, "illumina_pre.tsv"),
        os.path.join(SEQKIT_DIR, "ont.tsv"),
        os.path.join(SEQKIT_DIR, "assembly.tsv"),
        os.path.join(COV_DIR, "ONT_coverage.checked"),

# --- Merge & QC (raw) ---
rule merge_illumina:
    input: R1=illumina_r1(), R2=illumina_r2()
    output: R1=MERGE_R1, R2=MERGE_R2
    log: os.path.join(LOGS, "merge_illumina.log")
    conda: "envs/assembly.yml"
    threads: 1
    message: get_log_message(1, TOTAL_STEPS, "Merging Illumina lanes (R1/R2)")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{WORK}";
        if [ -z "{input.R1}" ] || [ -z "{input.R2}" ]; then
          echo "No Illumina reads found in {ILL}" >&2;
          : > "{output.R1}"; : > "{output.R2}";
          exit 0;
        fi;
        tmp1="{WORK}/.r1.tmp.fq"; tmp2="{WORK}/.r2.tmp.fq";
        >"$tmp1"; >"$tmp2";
        for f in {input.R1}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp1"; done;
        for f in {input.R2}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp2"; done;
        pigz -n -c "$tmp1" > "{output.R1}";
        pigz -n -c "$tmp2" > "{output.R2}";
        rm -f "$tmp1" "$tmp2";
        ) 2>&1 | tee "{log}"
        """

rule fastqc_pre:
    input: r1=MERGE_R1, r2=MERGE_R2
    output: directory(FASTQC_PRE_DIR)
    log: os.path.join(LOGS, "fastqc_pre.log")
    conda: "envs/qc.yml"
    threads: 2
    message: get_log_message(2, TOTAL_STEPS, "FastQC on raw merged Illumina")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{output}";
        if [ -s "{input.r1}" ]; then fastqc -t {threads} -o "{output}" "{input.r1}"; fi;
        if [ -s "{input.r2}" ]; then fastqc -t {threads} -o "{output}" "{input.r2}"; fi;
        ) 2>&1 | tee "{log}"
        """

rule kraken2_illumina:
    input: r1=MERGE_R1, r2=MERGE_R2
    output: kraken=os.path.join(KRAKEN2_DIR, "illumina.kraken2.txt"), report=os.path.join(KRAKEN2_DIR, "illumina.report.txt")
    params: db=KRAKEN2_DB_PATH
    log: os.path.join(LOGS, "kraken2_illumina.log")
    resources: mem_mb=32000
    conda: "envs/qc.yml"
    threads: THREADS
    message: get_log_message(3, TOTAL_STEPS, "Kraken2 contamination (raw Illumina)")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{KRAKEN2_DIR}";
        if [ -s "{input.r1}" ] && [ -s "{input.r2}" ]; then
            kraken2 --db "{params.db}" --threads {threads} --paired \
                --report "{output.report}" "{input.r1}" "{input.r2}" > "{output.kraken}";
        else
            echo "Skipping Kraken2 on Illumina (no inputs).";
            : > "{output.kraken}"; : > "{output.report}";
        fi
        ) 2>&1 | tee "{log}"
        """

rule merge_ont:
    input: files=nanopore_files()
    output: ONT_MERGED
    log: os.path.join(LOGS, "merge_ont.log")
    conda: "envs/assembly.yml"
    threads: 1
    message: get_log_message(4, TOTAL_STEPS, "Merging ONT FASTQ files")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{WORK}";
        if [ -z "{input.files}" ]; then echo "No ONT files in {ONT}" >&2; exit 1; fi;
        tmp="{WORK}/.ont_tmp.fastq"; > "$tmp";
        for f in {input.files}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp"; done;
        pigz -n -c "$tmp" > "{output}";
        rm -f "$tmp";
        ) 2>&1 | tee "{log}"
        """

rule nanoplot_raw:
    input: ONT_MERGED
    output: directory(NP_RAW_DIR)
    log: os.path.join(LOGS, "nanoplot_raw.log")
    conda: "envs/qc.yml"
    threads: 2
    message: get_log_message(5, TOTAL_STEPS, "NanoPlot on raw ONT")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{output}";
        NanoPlot --fastq "{input}" -o "{output}" --threads {threads} --tsv_stats;
        ) 2>&1 | tee "{log}"
        """

rule kraken2_nanopore:
    input: ont=ONT_MERGED
    output: kraken=os.path.join(KRAKEN2_DIR, "nanopore.kraken2.txt"), report=os.path.join(KRAKEN2_DIR, "nanopore.report.txt")
    params: db=KRAKEN2_DB_PATH
    log: os.path.join(LOGS, "kraken2_nanopore.log")
    resources: mem_mb=32000
    conda: "envs/qc.yml"
    threads: THREADS
    message: get_log_message(6, TOTAL_STEPS, "Kraken2 contamination (raw ONT)")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{KRAKEN2_DIR}";
        if [ -s "{input.ont}" ]; then
            kraken2 --db "{params.db}" --threads {threads} \
                --report "{output.report}" "{input.ont}" > "{output.kraken}";
        else
            echo "Skipping Kraken2 on ONT (no inputs).";
            : > "{output.kraken}"; : > "{output.report}";
        fi
        ) 2>&1 | tee "{log}"
        """

# --- Whitelist decontamination (KrakenTools) ---
rule clean_illumina_krakentools:
    input:
        r1 = MERGE_R1,
        r2 = MERGE_R2,
        kr = os.path.join(KRAKEN2_DIR, "illumina.kraken2.txt"),
        report = os.path.join(KRAKEN2_DIR, "illumina.report.txt")
    output:
        r1 = CLEAN_R1,
        r2 = CLEAN_R2
    params:
        taxid = lambda wc: str(TARGET_TAXID) if TARGET_TAXID else ""
    log: os.path.join(LOGS, "clean_illumina.log")
    conda: "envs/qc.yml"
    threads: 4
    message: get_log_message(7, TOTAL_STEPS, "Whitelist Illumina reads (KrakenTools)")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{WORK}";
        if [ -z "{params.taxid}" ]; then
          echo "[WARN] target_taxid not set; passing Illumina reads through." >&2;
          ln -sf "$(realpath "{input.r1}")" "{output.r1}" 2>/dev/null || cp -a "{input.r1}" "{output.r1}";
          ln -sf "$(realpath "{input.r2}")" "{output.r2}" 2>/dev/null || cp -a "{input.r2}" "{output.r2}";
          exit 0;
        fi

        tmp1="{WORK}/.illumina_R1.clean.fq"
        tmp2="{WORK}/.illumina_R2.clean.fq"

        extract_kraken_reads.py \
          -k "{input.kr}" \
          -r "{input.report}" \
          -s1 "{input.r1}" -s2 "{input.r2}" \
          --taxid "{params.taxid}" --include-children \
          -o "$tmp1" -o2 "$tmp2" --fastq-output

        pigz -n -f "$tmp1"; pigz -n -f "$tmp2"
        mv "$tmp1.gz" "{output.r1}"
        mv "$tmp2.gz" "{output.r2}"

        # sanity check: ensure valid gzip
        zcat -t "{output.r1}" >/dev/null
        zcat -t "{output.r2}" >/dev/null
        ) 2>&1 | tee "{log}"
        """

rule clean_ont_krakentools:
    input:
        ont = ONT_MERGED,
        kr  = os.path.join(KRAKEN2_DIR, "nanopore.kraken2.txt"),
        report = os.path.join(KRAKEN2_DIR, "nanopore.report.txt")
    output:
        ONT_CLEAN
    params:
        taxid = lambda wc: str(TARGET_TAXID) if TARGET_TAXID else ""
    log: os.path.join(LOGS, "clean_ont.log")
    conda: "envs/qc.yml"
    threads: 4
    message: get_log_message(8, TOTAL_STEPS, "Whitelist ONT reads (KrakenTools)")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{WORK}";
        if [ -z "{params.taxid}" ]; then
          echo "[WARN] target_taxid not set; passing ONT reads through." >&2;
          ln -sf "$(realpath "{input.ont}")" "{output}" 2>/dev/null || cp -a "{input.ont}" "{output}";
          exit 0;
        fi

        tmp="{WORK}/.ont_clean.fastq"

        extract_kraken_reads.py \
          -k "{input.kr}" \
          -r "{input.report}" \
          -s "{input.ont}" \
          --taxid "{params.taxid}" --include-children \
          -o "$tmp" --fastq-output

        pigz -n -f "$tmp"
        mv "$tmp.gz" "{output}"

        # sanity check: ensure valid gzip
        zcat -t "{output}" >/dev/null
        ) 2>&1 | tee "{log}"
        """

# --- Pre-process ---
rule fastp_trim:
    input: R1=CLEAN_R1, R2=CLEAN_R2
    output: R1=TRIM_R1, R2=TRIM_R2, html=os.path.join(QC,"fastqc","fastp.html"), json=os.path.join(QC,"fastqc","fastp.json")
    log: os.path.join(LOGS, "fastp_trim.log")
    conda: "envs/assembly.yml"
    threads: THREADS
    message: get_log_message(9, TOTAL_STEPS, "Illumina trimming/QC with fastp")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{TRIM}" "$(dirname "{output.html}")";
        if [ -s "{input.R1}" ] && [ -s "{input.R2}" ]; then
          fastp -i "{input.R1}" -I "{input.R2}" -o "{output.R1}" -O "{output.R2}" --thread {threads} --html "{output.html}" --json "{output.json}";
        else
          : > "{output.R1}"; : > "{output.R2}"; : > "{output.html}"; : > "{output.json}";
        fi
        ) 2>&1 | tee "{log}"
        """

rule filtlong_ont:
    input: ONT_CLEAN
    output: ONT_FILTERED
    log: os.path.join(LOGS, "filtlong_ont.log")
    conda: "envs/assembly.yml"
    threads: THREADS
    params: keep=KEEP_PCT
    message: get_log_message(10, TOTAL_STEPS, "ONT filtering (Filtlong)")
    shell:
        r"""
        (
        set -euo pipefail;
        filtlong --min_length 1000 --keep_percent {params.keep} "{input}" | pigz -n > "{output}";
        ) 2>&1 | tee "{log}"
        """

rule nanoplot_filt:
    input: ONT_FILTERED
    output: directory(NP_FILT_DIR)
    log: os.path.join(LOGS, "nanoplot_filt.log")
    conda: "envs/qc.yml"
    threads: 2
    message: get_log_message(11, TOTAL_STEPS, "NanoPlot on filtered ONT")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{output}";
        NanoPlot --fastq "{input}" -o "{output}" --threads {threads} --tsv_stats;
        ) 2>&1 | tee "{log}"
        """

# --- Assembly ---
rule unicycler:
    input: ont=ONT_FILTERED, r1=TRIM_R1, r2=TRIM_R2
    output: dir=directory(UNI_DIR), asm=UNI_ASM
    log: os.path.join(LOGS, "unicycler.log")
    resources: mem_mb=32000
    conda: "envs/assembly.yml"
    threads: THREADS
    message: get_log_message(12, TOTAL_STEPS, "Hybrid assembly with Unicycler")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{output.dir}" "{LOGS}";
        if [ -s "{input.r1}" ] && [ -s "{input.r2}" ]; then
          unicycler -l "{input.ont}" -1 "{input.r1}" -2 "{input.r2}" -o "{output.dir}" -t {threads} --min_fasta_length 1000;
        else
          unicycler -l "{input.ont}" -o "{output.dir}" -t {threads} --min_fasta_length 1000;
        fi;
        [ -s "{output.asm}" ]
        ) 2>&1 | tee "{log}"
        """

# --- Polishing ---
rule collect_unicycler_contigs:
    input: UNI_ASM
    output: ALL_CONTIGS
    threads: 1
    message: get_log_message(13, TOTAL_STEPS, "Collecting assembled contigs")
    shell: "set -euo pipefail; cp '{input}' '{output}';"

rule racon_round0:
    input: ALL_CONTIGS
    output: ROUND0
    threads: 1
    message: get_log_message(14, TOTAL_STEPS, "Seed Racon polishing (round 0)")
    shell: "set -euo pipefail; cp '{input}' '{output}';"

rule racon_map:
    input: asm=lambda wc: os.path.join(WORK, f"assembly.racon{int(wc.i)-1}.fasta"), ont=ONT_FILTERED
    output: paf=os.path.join(WORK, "ont_vs_racon{i}.paf")
    log: os.path.join(LOGS, "racon_map_{i}.log")
    conda: "envs/assembly.yml"
    threads: THREADS
    message: "Minimap2 mapping for Racon round {wildcards.i}"
    shell: "set -euo pipefail; (minimap2 -t {threads} -x map-ont '{input.asm}' '{input.ont}' > '{output.paf}') 2>&1 | tee '{log}';"

rule racon_consensus:
    input: ont=ONT_FILTERED, paf=lambda wc: os.path.join(WORK, f"ont_vs_racon{wc.i}.paf"), asm=lambda wc: os.path.join(WORK, f"assembly.racon{int(wc.i)-1}.fasta")
    output: out=os.path.join(WORK, "assembly.racon{i}.fasta")
    log: os.path.join(LOGS, "racon_consensus_{i}.log")
    resources: mem_mb=32000
    conda: "envs/assembly.yml"
    threads: THREADS
    message: "Racon consensus for round {wildcards.i}"
    shell: "set -euo pipefail; (racon -t {threads} '{input.ont}' '{input.paf}' '{input.asm}' > '{output.out}') 2>&1 | tee '{log}';"

rule racon_rounds:
    input: expand(os.path.join(WORK, "assembly.racon{i}.fasta"), i=range(1, RACON_ROUNDS + 1))

rule medaka:
    input: RACON_LAST, ONT_FILTERED
    output: MEDAKA_OUT
    log: os.path.join(LOGS, "medaka.log")
    resources: mem_mb=32000
    conda: "envs/medaka.yml"
    threads: THREADS
    params: model=MEDAKA_MODEL, outdir=MEDAKA_DIR
    message: get_log_message(15, TOTAL_STEPS, "Medaka polishing")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{params.outdir}";
        if [ -n "{params.model}" ]; then
          medaka_consensus -t {threads} -i "{input[1]}" -d "{input[0]}" -o "{params.outdir}" -m "{params.model}";
        else
          medaka_consensus -t {threads} -i "{input[1]}" -d "{input[0]}" -o "{params.outdir}";
        fi
        ) 2>&1 | tee "{log}"
        """

# Illumina mapping -> produce coord-sorted BAM, then derive name-sorted BAM
rule illumina_map:
    input: ref=MEDAKA_OUT, r1=TRIM_R1, r2=TRIM_R2
    output:
        bam=os.path.join(WORK, "illumina.bam"),
        bai=os.path.join(WORK, "illumina.bam.bai"),
        namesort=NAME_BAM
    log: os.path.join(LOGS, "illumina_map.log")
    conda: "envs/polishing.yml"
    threads: THREADS
    message: get_log_message(16, TOTAL_STEPS, "Illumina mapping (BWA-MEM; all alignments) + name-sort")
    shell:
        r"""
        (
        set -euo pipefail
        if [ -s "{input.r1}" ] && [ -s "{input.r2}" ]; then
          bwa index "{input.ref}"
          # -a keeps all alignments that Polypolish requires
          bwa mem -t {threads} -a -x intractg "{input.ref}" "{input.r1}" "{input.r2}" \
            | samtools sort -@ {threads} -O BAM -o "{output.bam}"
          samtools index "{output.bam}"
          samtools sort -n -@ {threads} -O BAM "{output.bam}" -o "{output.namesort}"
        else
          samtools faidx "{input.ref}"
          samtools view -b -T "{input.ref}" -h -o "{output.bam}" /dev/null
          samtools index "{output.bam}"
          samtools view -b -T "{input.ref}" -h -o "{output.namesort}" /dev/null
        fi
        ) 2>&1 | tee "{log}"
        """

# Polypolish -> use name-sorted BAM; v0.5.0 syntax has no subcommand
rule polypolish:
    input:
        ref=MEDAKA_OUT,
        bam=NAME_BAM  # name-sorted (SO:queryname)
    output: POLISHED
    log: os.path.join(LOGS, "polypolish.log")
    resources: mem_mb=32000
    conda: "envs/polishing.yml"
    threads: 1
    message: get_log_message(17, TOTAL_STEPS, "Illumina-based polishing (Polypolish)")
    shell:
        r"""
        (
        set -euo pipefail

        # ensure name-sorted; fix on the fly if not
        if ! samtools view -H "{input.bam}" | grep -q '^@HD.*SO:queryname'; then
          echo "[WARN] {input.bam} not name-sorted; resorting on the fly." >&2
          samtools sort -n -@1 -O BAM "{input.bam}" -o "{input.bam}.tmp"
          mv "{input.bam}.tmp" "{input.bam}"
        fi

        tmp_sam="{WORK}/.polypolish.sam"
        tmp_out="{output}.tmp"
        trap 'rm -f "$tmp_sam" "$tmp_out"' EXIT

        # Write SAM with headers to a temp file (PP 0.5.x needs a real file; not '-')
        samtools view -h "{input.bam}" > "$tmp_sam"

        # Count alignments (reads) in SAM
        read_count=$(samtools view -c "$tmp_sam" || echo 0)

        if [ "$read_count" -gt 0 ]; then
          polypolish "{input.ref}" "$tmp_sam" > "$tmp_out"
          mv "$tmp_out" "{output}"
        else
          echo "[WARN] No alignments; copying Medaka consensus." >&2
          cp "{input.ref}" "{output}"
        fi
        ) 2>&1 | tee "{log}"
        """

rule final_symlink:
    input: POLISHED
    output: FINAL
    threads: 1
    message: get_log_message(18, TOTAL_STEPS, "Creating final polished assembly")
    shell: "set -euo pipefail; cp '{input}' '{output}';"

# --- QC & summaries ---
rule quast:
    input: FINAL
    output: report=os.path.join(QUAST_DIR, "report.txt")
    log: os.path.join(LOGS, "quast.log")
    conda: "envs/qc.yml"
    threads: THREADS
    params: outdir=QUAST_DIR
    message: get_log_message(19, TOTAL_STEPS, "QUAST: assembly quality")
    shell: "set -euo pipefail; (quast.py -o '{params.outdir}' -t {threads} '{input}') 2>&1 | tee '{log}';"

rule busco:
    input: fasta=os.path.abspath(FINAL)
    output: summary=os.path.join(BUSCO_DIR, f"busco_{ISOLATE}", "short_summary.txt")
    log: os.path.join(LOGS, "busco.log")
    resources: mem_mb=32000
    conda: "envs/qc.yml"
    threads: THREADS
    params: lineage=BUSCO_LINEAGE, outdir=BUSCO_DIR, outname=f"busco_{ISOLATE}"
    message: get_log_message(20, TOTAL_STEPS, "BUSCO: biological completeness")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{LOGS}";
        busco -i "{input.fasta}" -o "{params.outname}" --out_path "{params.outdir}" -m genome -l "{params.lineage}" -c {threads} --force;
        rundir="{params.outdir}/{params.outname}";
        sfile=$(ls -1 "$rundir"/short_summary*.txt 2>/dev/null | head -n1 || true);
        if [ -n "$sfile" ] ; then
          mkdir -p "$(dirname "{output.summary}")";
          cp -f "$sfile" "{output.summary}";
        else
          mkdir -p "$(dirname "{output.summary}")" && : > "{output.summary}";
        fi
        ) 2>&1 | tee "{log}"
        """

rule seqkit_stats:
    input: i1=MERGE_R1, i2=MERGE_R2, ont_raw=ONT_MERGED, ont_filt=ONT_FILTERED, asm=FINAL
    output: illumina_pre=os.path.join(SEQKIT_DIR, "illumina_pre.tsv"), ont_tbl=os.path.join(SEQKIT_DIR, "ont.tsv"), asm_tbl=os.path.join(SEQKIT_DIR, "assembly.tsv")
    log: os.path.join(LOGS, "seqkit_stats.log")
    conda: "envs/qc.yml"
    threads: 1
    message: get_log_message(21, TOTAL_STEPS, "seqkit stats summaries")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{SEQKIT_DIR}";
        if [ -s "{input.i1}" ] || [ -s "{input.i2}" ]; then
          seqkit stats -a "{input.i1}" "{input.i2}" > "{output.illumina_pre}";
        else
          : > "{output.illumina_pre}";
        fi;
        seqkit stats -a "{input.ont_raw}" "{input.ont_filt}" > "{output.ont_tbl}";
        seqkit stats -a "{input.asm}" > "{output.asm_tbl}";
        ) 2>&1 | tee "{log}"
        """

rule coverage_warn:
    input: ont=ONT_FILTERED, asm=FINAL
    output: touch(os.path.join(COV_DIR,"ONT_coverage.checked"))
    conda: "envs/qc.yml"
    threads: 1
    message: get_log_message(22, TOTAL_STEPS, "Coverage check (warn-only)")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{COV_DIR}"

        # Extract total bp from the 2nd row of seqkit stats table
        ont_bp=$(seqkit stats -a -T "{input.ont}" | awk 'NR==2 {{gsub(/,/,"",$5); print $5}}' || echo 0)
        asm_bp=$(seqkit stats -a -T "{input.asm}" | awk 'NR==2 {{gsub(/,/,"",$5); print $5}}' || echo 0)

        [ -n "$ont_bp" ] || ont_bp=0
        [ -n "$asm_bp" ] || asm_bp=0

        if [ "$asm_bp" -gt 0 ]; then
          cov=$(( ont_bp / asm_bp ))
        else
          cov=0
        fi

        printf "ONT bp: %s  Assembly bp: %s  ONT_coverage(~x): %s\n" "$ont_bp" "$asm_bp" "$cov" \
          | tee "{COV_DIR}/ONT_coverage.txt"

        : > "{output}"
        """

rule multiqc:
    input:
        quast_report=os.path.join(QUAST_DIR, "report.txt"),
        busco=os.path.join(BUSCO_DIR, f"busco_{ISOLATE}", "short_summary.txt"),
        kraken_illumina_report=os.path.join(KRAKEN2_DIR, "illumina.report.txt"),
        kraken_nanopore_report=os.path.join(KRAKEN2_DIR, "nanopore.report.txt"),
        prokka_gff=os.path.join(PROKKA_DIR, f"{PROKKA_PREFIX}.gff"),  # <--- NEW: wait for Prokka
    output: os.path.join(MQC_DIR, "multiqc_report.html")
    log: os.path.join(LOGS, "multiqc.log")
    conda: "envs/qc.yml"
    threads: 1
    message: get_log_message(23, TOTAL_STEPS, "MultiQC aggregation")
    shell:
        "set -euo pipefail; "
        "(mkdir -p '{MQC_DIR}'; "
        " multiqc -f -n multiqc_report.html -o '{MQC_DIR}' '{ISOD}') "
        "2>&1 | tee '{log}';"


# --- Annotation ---
rule prokka:
    input: fasta=os.path.abspath(FINAL)
    output: gff=os.path.join(PROKKA_DIR, f"{PROKKA_PREFIX}.gff"), faa=PROKKA_FAA
    log: os.path.join(LOGS, "prokka.log")
    resources: mem_mb=32000
    conda: "envs/annotation.yml"
    threads: THREADS
    params: outdir=PROKKA_DIR, prefix=PROKKA_PREFIX
    message: get_log_message(24, TOTAL_STEPS, "Prokka genome annotation")
    shell:
        r"""
        (
        set -euo pipefail; mkdir -p "{LOGS}";
        prokka "{input.fasta}" --outdir "{params.outdir}" --prefix "{params.prefix}" --cpus {threads} --force
        ) 2>&1 | tee "{log}"
        """


# --- Final Integrated Report ---
report:
    """
    =================================
    Bacterial Assembly Report
    =================================

    **Isolate**: {ISOLATE}

    **Date**: {datetime.now().strftime("%Y-%m-%d")}

    This report summarizes the hybrid assembly, quality control, and functional analysis of the specified bacterial isolate.

    Quick Links
    -----------
    - `Full MultiQC Report <{os.path.relpath(os.path.join(MQC_DIR, 'multiqc_report.html'))}>`_

    Assembly Quality (QUAST)
    ------------------------
    .. csv-table:: QUAST Summary
       :file: {os.path.join(QUAST_DIR, "report.tsv")}
       :header-rows: 1

    Genome Completeness (BUSCO)
    ---------------------------
    .. raw:: html
       :file: {os.path.join(BUSCO_DIR, f"busco_{ISOLATE}", "short_summary.specific.{BUSCO_LINEAGE}.busco_{ISOLATE}.html")}

    """
