# 1. Use stable Ubuntu LTS base
FROM ubuntu:22.04

# 2. Install essential packages, including HTTPS certs for wget
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    wget curl git file sudo build-essential ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# 3. Install Miniconda (not quiet, so we see errors if any)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/miniconda/condabin:/miniconda/bin:${PATH}"

# 4. Optional: Install Linuxbrew (if needed)
# RUN /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 5. Clone your code repo
RUN git clone https://gitlab.igem.org/2023/software-tools/dtu-denmark.git

# 6. Create Conda environment
RUN conda env create --file dtu-denmark/environment.yml

# 7. Install science tools (AutoDock Vina, OpenBabel, GROMACS)
RUN apt-get update && \
    apt-get install -y autodock-vina openbabel cmake && \
    rm -rf /var/lib/apt/lists/*

# 8. Download & build GROMACS
RUN curl -O https://ftp.gromacs.org/pub/gromacs/gromacs-2021.4.tar.gz && \
    tar -xf gromacs-2021.4.tar.gz && \
    mkdir /gromacs-2021.4/build
WORKDIR /gromacs-2021.4/build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_USE_RDTSCP=OFF && \
    make && \
    make install && \
    rm ../../gromacs-2021.4.tar.gz && \
    rm -rf ../../gromacs-2021.4
RUN ln -s /usr/local/gromacs/bin/gmx /usr/bin/gmx

# 9. Setup project folder and create user
WORKDIR /dtu-denmark
RUN useradd -ms /bin/bash jupyter && \
    echo "jupyter ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/jupyter

# 10. Set permissions for Jupyter user
RUN chown -R jupyter:jupyter /miniconda/envs/AptaLoop && \
    chown -R jupyter:jupyter /dtu-denmark && \
    chmod -R 777 /dtu-denmark

# 11. Switch to user
USER jupyter

# 12. Register kernel
RUN /bin/bash -c "source /miniconda/bin/activate AptaLoop && \
    python -m ipykernel install --user --name AptaLoop --display-name 'Python (AptaLoop)'"

# 13. Run JupyterLab
EXPOSE 8888
CMD ["/bin/bash", "-c", "source /miniconda/bin/activate AptaLoop && jupyter lab --ip=0.0.0.0 --port=8888 --allow-root"]
