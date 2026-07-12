FROM python:3.10-slim

# System deps: build tools for python-igraph / leidenalg / scipy wheels,
# plus git for any VCS installs.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        libgl1 \
        libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Install Python deps first for better layer caching.
COPY requirements.txt /work/requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir "numpy<2" \
    && pip install --no-cache-dir -r requirements.txt

# Copy the project source.
COPY . /work

# Data and results are mounted / generated at runtime, not baked in.
CMD ["bash"]
