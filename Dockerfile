# Lightweight conda base with mamba
FROM condaforge/miniforge3:24.3.0-0

SHELL ["/bin/bash", "-lc"]

# System basics (optional but handy)
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tini && \
    rm -rf /var/lib/apt/lists/*

# Workdir
WORKDIR /workspace

# Copy env specs first to leverage Docker layer caching
COPY environment.yml /tmp/environment.yml
COPY requirements.txt /tmp/requirements.txt

# Create "maws" env. If environment.yml is missing packages, we top it up explicitly.
RUN mamba env create -n maws -f /tmp/environment.yml || mamba create -n maws -y python=3.11 && \
    mamba install -n maws -y -c conda-forge \
        ambertools \
        openmm \
        jupyterlab \
        ipykernel && \
    source activate maws && \
    pip install -r /tmp/requirements.txt && \
    mamba clean -afy

# Put the env on PATH so we don't need "conda activate" at runtime
ENV PATH="/opt/conda/envs/maws/bin:${PATH}"

# Copy repo (compose also mounts it; this lets the image run standalone too)
COPY . /workspace

# (Optional) Register a Jupyter kernel for this env
RUN python -m ipykernel install --sys-prefix --name maws --display-name "Python (MAWS)"

# Non-root user for Jupyter
RUN useradd -ms /bin/bash -u 1000 jupyter && \
    chown -R jupyter:jupyter /workspace /opt/conda

USER jupyter

EXPOSE 8888
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--NotebookApp.token=", "--notebook-dir=/workspace"]
