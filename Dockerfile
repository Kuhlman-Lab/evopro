FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04
ADD evopro/ /home/evopro/evopro/ 
RUN apt-get -y update
RUN apt-get -y install git
RUN git clone https://github.com/Kuhlman-Lab/proteinmpnn.git /home/evopro/proteinmpnn/
RUN git clone https://github.com/Kuhlman-Lab/alphafold.git /home/evopro/alphafold/
RUN apt-get install -y wget
RUN /home/evopro/alphafold/setup/download_alphafold_params.sh /home/evopro/alphafold/setup/

ENV PYTHONPATH /home/evopro/:/home/evopro/evopro/score_funcs/:/home/evopro/alphafold/:/home/evopro/alphafold/run/:/home/evopro/proteinmpnn/:/home/evopro/proteinmpnn/run/

ENV EVOPRO /home/evopro/evopro/
ENV ALPHAFOLD_RUN /home/evopro/alphafold/run/
ENV PROTEIN_MPNN_RUN /home/evopro/proteinmpnn/run/

# Install base utilities
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

ADD setup_conda.yaml /setup_conda.yaml
RUN conda env update -n base --file /setup_conda.yaml
RUN pip3 install --upgrade jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
RUN python3 -m pip install /home/evopro/alphafold/alphafold/
RUN python3 -m pip install numpy scipy matplotlib pandas

ENTRYPOINT ["python", "/home/evopro/evopro/run/run_geneticalg_gpus.py"]
