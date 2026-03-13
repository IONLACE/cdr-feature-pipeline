FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CDR_OUT_DIR=/data

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    build-essential \
    cargo \
    cmake \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Build fpocket from source (no apt package exists).
ARG FPOCKET_REPO=https://github.com/Discngine/fpocket.git
RUN git clone --depth 1 "${FPOCKET_REPO}" /tmp/fpocket \
        && if [ -f /tmp/fpocket/CMakeLists.txt ]; then \
                 cmake -S /tmp/fpocket -B /tmp/fpocket/build \
                 && cmake --build /tmp/fpocket/build --parallel $(nproc) \
                 && cp /tmp/fpocket/build/bin/fpocket /usr/local/bin/fpocket; \
             elif [ -f /tmp/fpocket/src/Makefile ]; then \
                 make -C /tmp/fpocket/src \
                 && cp /tmp/fpocket/bin/fpocket /usr/local/bin/fpocket; \
             else \
                 echo "Unsupported fpocket source layout" && ls -la /tmp/fpocket && exit 1; \
             fi \
    && rm -rf /tmp/fpocket

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
