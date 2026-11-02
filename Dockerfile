# Runtime image for annovar-nf: bcftools/htslib, Perl and Python.
# ANNOVAR is licensed and NOT included: pass --annovar_src/--annovar_db,
# Nextflow mounts them into the container automatically.
#
#   docker build -t ghcr.io/assetdaniyarov/annovar-nf:0.1.1 .
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="annovar-nf" \
      org.opencontainers.image.description="Runtime for the annovar-nf Nextflow pipeline (ANNOVAR not included)" \
      org.opencontainers.image.source="https://github.com/assetdaniyarov/annovar-nf" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba install -y -n base -c conda-forge procps-ng \
    && micromamba clean --all --yes

# Put the environment on PATH for non-interactive Nextflow tasks.
ENV PATH="/opt/conda/bin:${PATH}"
