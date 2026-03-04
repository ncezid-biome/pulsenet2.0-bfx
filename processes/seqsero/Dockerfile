FROM python:3.11-slim AS poetry
RUN pip install poetry==1.8.3
WORKDIR /app
COPY pyproject.toml  ./
RUN poetry export -f requirements.txt --output requirements.txt --without-hashes \
    && poetry install --no-root


FROM python:3.11-slim
SHELL ["/bin/bash", "-c"]

ARG MINIFORGE_NAME=Mambaforge
ARG MINIFORGE_VERSION=22.11.1-4


ENV CONDA_DIR=/opt/conda
ENV CONDA_BIN_DIR=${CONDA_DIR}/condabin
ENV CONDA_BIN_DIR2=${CONDA_DIR}/bin
ENV CONDA_BIN_DIR3=${CONDA_DIR}/envs/seqsero/bin
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_BIN_DIR}:${PATH}
ENV PATH=/opt/conda/envs/seqsero/bin/:${PATH}
ENV PATH=${CONDA_BIN_DIR2}:${PATH}
ENV PATH=${CONDA_BIN_DIR3}:${PATH}
RUN apt-get update && apt-get install unzip wget procps git --no-install-recommends -y


# Install mamba from mambaforge
# https://github.com/conda-forge/miniforge#mambaforge
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh
RUN /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR}
RUN mamba install libarchive=3.5.2=hb890918_2 -y
# Install mamba dependencies
COPY conda.yaml /env/conda.yaml
RUN mamba create -n seqsero
RUN mamba env update -n seqsero -f /env/conda.yaml

ENV HOME=/app
WORKDIR $HOME

# Install Process code & its dependencies

ARG PIP_REGISTRY_PASSWORD
ARG PIP_REGISTRY_USERNAME
ENV PYTHONPATH="/usr/local/bin/"
ENV PYTHONPATH="${PYTHONPATH}:${HOME}"


RUN python3.11 -m pip install --upgrade pip
COPY --from=poetry /app/requirements.txt .
RUN pip3.11 install --no-deps -r requirements.txt

# edit stringMLST to avoid "read-only file system" error from trying to write to conda bin dir
RUN sed -i 's/log = dbPrefix/log = "stringmlst_log"/' /opt/conda/envs/seqsero/bin/stringMLST.py

COPY src/ src/
