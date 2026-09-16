FROM python:3.11-slim

ARG MMC_VERSION=1.7.4

RUN apt-get update && \
    apt-get install -y git build-essential wget curl && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir git+https://github.com/AllenInstitute/cell_type_mapper.git@v${MMC_VERSION}

# Bake in the built-in reference-atlas assets so the WDL never has to fetch them from
# external storage at runtime. Source: Allen Institute's public MapMyCells S3 bucket,
# per https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/input_data_files/running_online_taxonomies_locally.md
ENV MAPMYCELLS_DATA_DIR=/opt/mapmycells/data
RUN mkdir -p ${MAPMYCELLS_DATA_DIR} && \
    # Whole Mouse Brain / ABC atlas taxonomy (CCN20230722)
    curl -fL -o ${MAPMYCELLS_DATA_DIR}/precomputed_stats_ABC_revision_230821.h5 \
        https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/precomputed_stats_ABC_revision_230821.h5 && \
    curl -fL -o ${MAPMYCELLS_DATA_DIR}/mouse_markers_230821.json \
        https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/mouse_markers_230821.json && \
    # Human MTG / SEA-AD taxonomy (CCN20230505) -- small enough that cell_type_mapper
    # finds markers at runtime, so no marker file is needed for this one
    curl -fL -o ${MAPMYCELLS_DATA_DIR}/precomputed_stats.20231120.sea_ad.MTG.h5 \
        https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5

# ponytail: gene-symbol -> Ensembl-ID translation (mmc_gene_mapper) is NOT baked in here.
# Building its db requires `bkbit download-ncbi-taxonomy`, which pulls a 15GB file --
# too large/slow for a routine image build (likely exceeds GH Actions runner disk).
# Left as an explicit follow-up; see MapMyCells.wdl's custom_gene_mapping_db input.

CMD ["/bin/bash"]
