FROM python:3.11-slim AS poetry
RUN pip install poetry==1.8.3
WORKDIR /app
COPY pyproject.toml  ./
RUN poetry export -f requirements.txt --output requirements.txt --without-hashes

ARG MINIFORGE_NAME=Mambaforge
ARG MINIFORGE_VERSION=22.11.1-4


ENV CONDA_DIR=/opt/conda
ENV CONDA_BIN_DIR=${CONDA_DIR}/condabin
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_BIN_DIR}:${PATH}

RUN apt-get update && apt-get install unzip wget uuid-runtime procps --no-install-recommends -y


# Install mamba from mambaforge
# https://github.com/conda-forge/miniforge#mambaforge
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh
RUN /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR}

# Install mamba dependencies
COPY env/ani.yaml /env/ani.yaml
RUN conda config --add channels bioconda \
    && conda config --add channels anaconda \
    &&  mamba env update -n ani -f=/env/ani.yaml && mamba clean -afy

# Install ani
RUN wget https://github.com/lskatz/ani-m/archive/refs/tags/v0.1.tar.gz && \
 tar xzf v0.1.tar.gz && \
 rm -v v0.1.tar.gz
COPY ani-m_1.pl /app
ENV LC_ALL=C \
 PATH="${PATH}:/ani-m-0.1:/mash-Linux64-v2.3"

#install ncbi datasets tool
RUN wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets && \
 chmod +x datasets && \
 mv -v datasets /usr/local/bin


ENV HOME=/app
WORKDIR $HOME
ENV PATH="/opt/conda/envs/ani/bin:${PATH}"
# Install Process code & its dependencies

ENV PYTHONPATH="${PYTHONPATH}:${HOME}"
RUN python3 -m pip install --upgrade pip
RUN pip3 install --no-deps -r requirements.txt
COPY src/ src/
