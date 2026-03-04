FROM python:3.11-slim AS poetry
RUN pip install poetry
WORKDIR /app
COPY pyproject.toml poetry.lock ./
RUN poetry export -f requirements.txt --output requirements.txt --without-hashes


FROM python:3.11-slim
SHELL ["/bin/bash", "-c"]

ARG MINIFORGE_NAME=Mambaforge
ARG MINIFORGE_VERSION=22.11.1-4


ENV CONDA_DIR=/opt/conda/
ENV CONDA_BIN_DIR=${CONDA_DIR}/condabin
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_BIN_DIR}:${PATH}

RUN apt-get update && apt-get install unzip wget uuid-runtime procps --no-install-recommends -y

# Install mamba from mambaforge
# https://github.com/conda-forge/miniforge#mambaforge
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh
RUN /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR}

# Install mamba dependencies
COPY conda.yaml /env/conda.yaml
# RUN conda config --add channels bioconda
RUN mamba create -n midas
RUN mamba env update -n midas -f /env/conda.yaml


ENV HOME=/app
WORKDIR $HOME
ADD run.py /app
ADD lib /app/lib
RUN chmod +x run.py && \
	ln -s /app/run.py /usr/bin/
RUN apt-get install -y curl 


RUN curl -k -L -o azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux && \     
     tar zxvf azcopy.tar.gz --strip-components=1 && \
     rm azcopy.tar.gz && \     
     chmod +x azcopy && \     
     mv azcopy /usr/local/bin/azcopy
# Install Process code & its dependencies
RUN chmod 777 /opt/conda/envs/midas/bin/*
ENV PYTHONPATH="${PYTHONPATH}:${HOME}"
ENV PATH="/opt/conda/envs/midas/bin/:${PATH}"
ENV PYTHONPATH="/opt/conda/envs/midas/lib/python2.7/site-packages/$PYTHONPATH"
RUN update-alternatives --install /usr/local/bin/ python3 /usr/local/bin/python3.11 1


COPY --from=poetry /app/requirements.txt .
RUN cat requirements.txt
RUN pip3.11 install --no-deps -r requirements.txt

COPY src/ src/
