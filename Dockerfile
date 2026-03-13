FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CDR_OUT_DIR=/data

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    build-essential \
    cargo \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install fpocket via micromamba (conda-forge).
ENV PATH="/opt/fpocket/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/fpocket/lib"
RUN curl -fsSL "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64" \
        -o /usr/local/bin/micromamba \
    && chmod +x /usr/local/bin/micromamba \
    && micromamba install -y --no-rc -c conda-forge fpocket --prefix /opt/fpocket \
    && rm /usr/local/bin/micromamba

# Build cons-capra07 from source (needed for step 02 conservation scoring).
ARG CONS_CAPRA07_REPO=https://github.com/IONLACE/cons-capra07.git
RUN git clone --depth 1 "${CONS_CAPRA07_REPO}" /tmp/cons-capra07 \
    && cargo build --release --manifest-path /tmp/cons-capra07/Cargo.toml \
    && cp /tmp/cons-capra07/target/release/cons-capra07 /usr/local/bin/cons-capra07 \
    && rm -rf /tmp/cons-capra07

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /app/requirements.txt

COPY . /app
RUN mkdir -p /data

ENTRYPOINT ["python", "run_pipeline.py"]
