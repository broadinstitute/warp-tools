#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "Usage: index_phenotype_bed.sh INPUT MODALITY OUTPUT_DIR THREADS" >&2
  exit 2
fi

input="$1"
modality="$2"
output_dir="$3"
threads="$4"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

if [[ ! -f "$input" ]]; then
  echo "Phenotype input does not exist: $input" >&2
  exit 1
fi
if [[ -z "$modality" ]]; then
  echo "Modality cannot be empty." >&2
  exit 1
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "Threads must be a positive integer." >&2
  exit 1
fi
for command_name in awk bgzip gzip sort tabix; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "Required command is not available: $command_name" >&2
    exit 1
  fi
done

mkdir -p "$output_dir"
temporary_dir="$(mktemp -d "${TMPDIR:-/tmp}/phenotype-index.XXXXXX")"
trap 'rm -rf "$temporary_dir"' EXIT

normalized="$temporary_dir/normalized.tsv"
lookup_tsv="$temporary_dir/phenotype_lookup.tsv"
stats="$temporary_dir/stats.tsv"
indexed="$output_dir/phenotypes.bed.gz"
lookup="$output_dir/phenotype_lookup.tsv.gz"
qc="$output_dir/index_qc.tsv"
started_at="$(date +%s)"
input_bytes="$(wc -c < "$input" | tr -d ' ')"

case "$input" in
  *.gz|*.bgz)
    reader=(gzip -cd -- "$input")
    ;;
  *)
    reader=(cat -- "$input")
    ;;
esac

log "Reading and validating the $modality phenotype file: $input"
"${reader[@]}" | awk -F '\t' -v OFS='\t' -v lookup="$lookup_tsv" -v stats="$stats" '
  NR == 1 {
    expected = NF
    if (expected < 5) {
      print "Phenotype input must contain four metadata columns and at least one sample column." > "/dev/stderr"
      exit 1
    }
    $1 = "#chrom"
    $2 = "start"
    $3 = "end"
    $4 = "phenotype_id"
    print
    print "phenotype_id", "chrom", "start", "end" > lookup
    next
  }
  {
    if (NF != expected) {
      print "Phenotype input has an inconsistent field count at line " NR "." > "/dev/stderr"
      exit 1
    }
    if ($1 == "") {
      print "Phenotype chromosome cannot be empty at line " NR "." > "/dev/stderr"
      exit 1
    }
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || ($3 + 0) <= ($2 + 0)) {
      print "Phenotype coordinates are invalid at line " NR "." > "/dev/stderr"
      exit 1
    }
    if ($4 == "") {
      print "Phenotype ID cannot be empty at line " NR "." > "/dev/stderr"
      exit 1
    }
    if (seen[$4]++) {
      print "Phenotype ID is duplicated: " $4 > "/dev/stderr"
      exit 1
    }
    records++
    print
    print $4, $1, $2, $3 > lookup
  }
  END {
    if (records < 1) {
      print "Phenotype input contains no data rows." > "/dev/stderr"
      exit 1
    }
    print records, expected - 4 > stats
  }
' > "$normalized"

read -r phenotype_rows sample_columns < "$stats"
log "Validated $phenotype_rows phenotype rows and $sample_columns sample columns."

log "Sorting phenotype rows by genomic coordinate."
{
  head -n 1 "$normalized"
  tail -n +2 "$normalized" | LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k3,3n
} | bgzip -@ "$threads" -c > "$indexed"

log "Creating the tabix index."
tabix -f -p bed "$indexed"

log "Compressing the phenotype-coordinate lookup."
gzip -c "$lookup_tsv" > "$lookup"

IFS=$'\t' read -r first_id first_chrom first_start first_end < <(tail -n +2 "$lookup_tsv" | head -n 1)
log "Testing the tabix index with phenotype $first_id."
tabix "$indexed" "${first_chrom}:$((first_start + 1))-${first_end}" | \
  awk -F '\t' -v expected_id="$first_id" '
    $4 == expected_id { found = 1 }
    END { exit(found ? 0 : 1) }
'

elapsed_seconds="$(( $(date +%s) - started_at ))"
bgzf_bytes="$(wc -c < "$indexed" | tr -d ' ')"
lookup_bytes="$(wc -c < "$lookup" | tr -d ' ')"
{
  printf 'modality\tphenotype_rows\tsample_columns\tinput_bytes\tbgzf_bytes\tlookup_bytes\telapsed_seconds\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$modality" "$phenotype_rows" "$sample_columns" "$input_bytes" \
    "$bgzf_bytes" "$lookup_bytes" "$elapsed_seconds"
} > "$qc"

log "Phenotype indexing complete: $indexed"
