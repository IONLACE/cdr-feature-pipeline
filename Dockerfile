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
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Build fpocket from source.
ARG FPOCKET_REPO=https://github.com/Discngine/fpocket.git
RUN git clone --depth 1 "${FPOCKET_REPO}" /tmp/fpocket \
    && sed -i 's/char \*residue_string\[M_MAX_CUSTOM_POCKET_LEN\];/char residue_string[M_MAX_CUSTOM_POCKET_LEN];/' /tmp/fpocket/src/fparams.c \
    && sed -i 's/strcpy(&residue_string, pt);/strcpy(residue_string, pt);/' /tmp/fpocket/src/fparams.c \
    && make -C /tmp/fpocket qhull \
    && make -C /tmp/fpocket fpocket \
    && install -m 0755 /tmp/fpocket/bin/fpocket /usr/local/bin/fpocket \
    && rm -rf /tmp/fpocket

# Build cons-capra07 from source (needed for step 02 conservation scoring).
ARG CONS_CAPRA07_REPO=https://github.com/shin-kinos/cons-capra07.git
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
