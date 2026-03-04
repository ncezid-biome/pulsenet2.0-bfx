FROM python:3.11-slim AS poetry
RUN pip install poetry
WORKDIR /app
COPY pyproject.toml poetry.lock ./
RUN poetry export -f requirements.txt --output requirements.txt --without-hashes

FROM python:3.11-slim
SHELL ["/bin/bash", "-c"]

ENV HOME=/app
WORKDIR $HOME

# Install Process code & its dependencies

ENV PYTHONPATH="${PYTHONPATH}:${HOME}"
RUN python3.11 -m pip install --upgrade pip
COPY --from=poetry /app/requirements.txt .
COPY requirements.txt /app/requirements_ispcr.txt
RUN pip3.11 install --no-deps -r requirements.txt
RUN pip3.11 install --no-deps -r /app/requirements_ispcr.txt
COPY src/ src/
