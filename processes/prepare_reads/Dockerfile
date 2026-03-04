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
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_BIN_DIR}:${PATH}

RUN apt-get update && apt-get install unzip wget uuid-runtime procps --no-install-recommends -y


# Install mamba from mambaforge
# https://github.com/conda-forge/miniforge#mambaforge
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh
RUN /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR}

# Install mamba dependencies
COPY conda.yaml /env/cleanreads.yaml
RUN conda config --add channels bioconda \
    && conda config --add channels anaconda \
    &&  mamba env update -n base -f=/env/cleanreads.yaml && mamba clean -afy


ENV HOME=/app
WORKDIR $HOME

# Configure NCBI
RUN mkdir ${HOME}/.ncbi/ 
RUN cp /root/.ncbi/user-settings.mkfg ${HOME}/.ncbi/user-settings.mkfg \
    && printf '/LIBS/GUID = "%s"\n' `uuidgen` >> ${HOME}/.ncbi/user-settings.mkfg \
    && printf '/libs/cloud/accept_aws_charges = "false"\n' >> ${HOME}/.ncbi/user-settings.mkfg \
    && printf '/libs/cloud/report_instance_identity = "true"\n' >> ${HOME}/.ncbi/user-settings.mkfg \
    && printf '/repository/user/main/public/cache-disabled = "true"\n'>> ${HOME}/.ncbi/user-settings.mkfg

# Install Process code & its dependencies
ENV PYTHONPATH="${PYTHONPATH}:${HOME}"
RUN python3 -m pip install --upgrade pip
COPY --from=poetry /app/requirements.txt .
RUN cat requirements.txt
RUN pip3 install --no-deps -r requirements.txt
COPY src/ src/
