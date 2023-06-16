FROM ubuntu:latest
ADD alphafold/ /home/evopro/alphafold/
ADD evopro/ /home/evopro/evopro/ 
ADD proteinmpnn/ /home/evopro/proteinmpnn/

ENV PYTHONPATH /home/evopro/:/home/evopro/alphafold/

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
#RUN conda env create -n evopro -f /setup_conda.yaml
RUN conda env update -n base --file /setup_conda.yaml
#SHELL ["conda", "run", "-n", "evopro", "/bin/bash", "-c"]
RUN pip3 install --upgrade jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
RUN python3 -m pip install /home/evopro/alphafold/alphafold/
RUN python3 -m pip install numpy scipy matplotlib pandas
