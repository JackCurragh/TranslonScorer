FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# System deps for pysam/pyBigWig and friends
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    git \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Paths are relative to the repo root, which IS this repo — not a translonscorer/
# subdirectory. The earlier COPY translonscorer/... form was inherited from the
# monorepo layout and fails on any case-sensitive filesystem.
# Dependencies first, for layer caching.
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Then the package. See .dockerignore — data/ is multi-GB and must stay out.
COPY . /app
RUN pip install --no-cache-dir /app

ENTRYPOINT ["translonscorer"]
