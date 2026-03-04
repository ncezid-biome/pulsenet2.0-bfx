# Use a base image with a compatible Java runtime and essential build tools
ARG NXF_VER="24.10.0"
FROM nextflow/nextflow:${NXF_VER} AS builder

# Use a minimal base image for runtime
FROM adoptopenjdk:11-jre-hotspot

# Copy Nextflow binary from the builder stage
COPY --from=builder /usr/local/bin/nextflow /usr/local/bin/nextflow
RUN chmod 755 /usr/local/bin/nextflow
ENV PATH="/root/.nextflow:${PATH}"
ENV NXF_HOME="/root/.nextflow"
ENV NXF_PLUGINS_DIR=$NXF_HOME/plugins
ENV NXF_VER="24.10.0"
ENV NXF_OFFLINE="true"
ENV NXF_CLOUDCACHE_PATH='az://batch-inputs/nf_cloud_cache/'
ENV NXF_IGNORE_RESUME_HISTORY=true
RUN apt-get update && apt-get install -y curl 
RUN curl -s https://github.com/nextflow-io/nextflow/releases/download/v$NXF_VER/nextflow-$NXF_VER-all
# Set the working directory inside the container
RUN apt-get update && apt-get install -y curl 
RUN curl -k -L -o azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux && \     
     tar zxvf azcopy.tar.gz --strip-components=1 && \
     rm azcopy.tar.gz && \     
     chmod +x azcopy && \     
     mv azcopy /usr/local/bin/azcopy

WORKDIR /app


# Copy the Nextflow script and any necessary files to the container
RUN mkdir -p /root/.nextflow/plugins
COPY . /app
RUN nextflow -version
RUN nextflow plugin install nf-azure@1.10.1
RUN nextflow plugin install nf-weblog@1.1.2
RUN nextflow plugin install nf-cloudcache@0.4.2