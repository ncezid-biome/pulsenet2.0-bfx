FROM python:3.11-slim AS poetry
RUN pip install poetry==1.8.3
WORKDIR /app
COPY pyproject.toml poetry.lock ./
RUN poetry export -f requirements.txt --output requirements.txt --without-hashes

FROM python:3.11-slim
SHELL ["/bin/bash", "-c"]

ARG MINIFORGE_NAME=Mambaforge
ARG MINIFORGE_VERSION=22.11.1-4

ENV CONDA_DIR=/opt/conda
ENV CONDA_BIN_DIR=${CONDA_DIR}/condabin
ENV CONDA_BIN_DIR2=${CONDA_DIR}/bin
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_BIN_DIR}:${PATH}
ENV PATH=${CONDA_BIN_DIR2}:${PATH}

RUN apt-get update \
	&& apt-get install unzip wget procps --no-install-recommends -y

# Install mamba from mambaforge
# https://github.com/conda-forge/miniforge#mambaforge
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh
RUN /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR}

# Install mamba dependencies
COPY conda.yaml /env/conda.yaml
RUN conda config --add channels bioconda \
    && conda config --add channels anaconda \
    &&  mamba env update -n base -f=/env/conda.yaml && mamba clean -afy \
    && find /opt/conda/ -maxdepth 1 -type f -name "*" | grep -P -v "[A-Z]" | while read file; do mv $file ${file%/*}/bin/; done

ENV HOME /app
WORKDIR $HOME

# Install Process code & its dependencies

ARG PIP_REGISTRY_PASSWORD
ARG PIP_REGISTRY_USERNAME
ENV PYTHONPATH "${PYTHONPATH}:${HOME}"
RUN python3.11 -m pip install --upgrade pip
# RUN pip config set global.index-url https://$PIP_REGISTRY_USERNAME:$PIP_REGISTRY_PASSWORD@pdx-platform-224016688692.d.codeartifact.eu-west-1.amazonaws.com/pypi/pdx-python-libs/simple/
COPY --from=poetry /app/requirements.txt .
COPY requirements.txt /app/requirements_module.txt
RUN cat requirements.txt
RUN pip3.11 install --no-deps -r requirements.txt
RUN pip3.11 install --no-deps -r /app/requirements_module.txt
COPY src/ src/
