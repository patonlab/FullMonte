FROM python:2.7.9-slim

RUN apt-get update && apt-get install -y \
                python-numpy \
                python-rdkit \
    --no-install-recommends && rm -rf /var/lib/apt/lists/*
