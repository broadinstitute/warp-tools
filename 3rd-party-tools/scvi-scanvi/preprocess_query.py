"""Preprocess one query sample for label transfer.

Moved here from the PreprocessFilter task that used to be inline Python in scANVI.wdl; the steps
are unchanged, but the result is returned in memory instead of written to disk:
  GEX   keep STARsolo cell calls (star_IsCell, added as True when missing), drop genes and cells
        with fewer than 3 counts, keep the counts in layers['counts'].
  ATAC  (multiome only) index by gex_barcodes, keep the barcodes shared with GEX, and convert the
        cell-by-bin matrix to a gene-activity matrix for the chosen genome.
  Both  the batch label, the placeholder label 'Unknown', and the modality tag used downstream.
"""
import scanpy as sc

GEX_BATCH = "pd-multiome_sci_gex"
ATAC_BATCH = "pd-multiome_sci_atac"


def _gene_annotation(snap, genome):
    annotations = {
        "hg38": getattr(snap.genome, "hg38", None),
        "mm10": getattr(snap.genome, "mm10", None),
        "mm39": getattr(snap.genome, "mm39", None) or getattr(snap.genome, "GRCm39", None),
    }
    if annotations.get(genome) is None:
        raise SystemExit(f"ERROR: unsupported/unavailable genome '{genome}' (expected hg38|mm10|mm39)")
    return annotations[genome]


def preprocess_query(gex_path, atac_path=None, genome="hg38"):
    """Return (gex, atac_activity) ready for label transfer; atac_activity is None without ATAC."""
    gex = sc.read_h5ad(gex_path)
    print(f"  GEX loaded: {gex.shape}", flush=True)
    if "star_IsCell" not in gex.obs.columns:
        print("  Patching: adding star_IsCell = True (column was missing)", flush=True)
        gex.obs["star_IsCell"] = True
    gex = gex[gex.obs["star_IsCell"] == True]
    print(f"  After star_IsCell filter: {gex.shape}", flush=True)
    sc.pp.filter_genes(gex, min_counts=3)
    sc.pp.filter_cells(gex, min_counts=3)
    print(f"  After min_counts filter: {gex.shape}", flush=True)
    gex.layers["counts"] = gex.X.copy()

    atac_activity = None
    if atac_path:
        import snapatac2 as snap
        atac = snap.read(atac_path, backed=None)
        print(f"  ATAC loaded: {atac.shape}", flush=True)
        if "gex_barcodes" not in atac.obs.columns:
            print("  Patching: adding gex_barcodes from obs index (column was missing)", flush=True)
            atac.obs["gex_barcodes"] = atac.obs.index
        atac.obs = atac.obs.set_index("gex_barcodes")
        shared = atac.obs.index.intersection(gex.obs.index)
        print(f"  Shared barcodes between GEX and ATAC: {len(shared)}", flush=True)
        atac = atac[atac.obs.index.isin(shared)].copy()
        gex = gex[gex.obs.index.isin(shared)].copy()
        atac.obs["batch"] = ATAC_BATCH
        print(f"  Converting ATAC cell-by-bin to gene activity ({genome})...", flush=True)
        atac_activity = snap.pp.make_gene_matrix(atac, gene_anno=_gene_annotation(snap, genome))
        atac_activity.obs["final_annotation"] = "Unknown"
        atac_activity.obs["modality"] = "atac_unannotated"
        print(f"  Gene activity matrix: {atac_activity.shape}", flush=True)
    else:
        gex = gex.copy()

    gex.obs["batch"] = GEX_BATCH
    gex.obs["final_annotation"] = "Unknown"
    gex.obs["modality"] = "rna_unannotated"
    return gex, atac_activity
