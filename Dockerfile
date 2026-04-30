FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt || \
    pip install --no-cache-dir flask flask-cors pyyaml pandas pyarrow jinja2 google-cloud-bigquery google-cloud-storage requests tqdm

COPY . .

RUN mkdir -p data/vcf data/annotated data/enriched data/classified data/phenotype reports eval logs

EXPOSE 8080

CMD ["python3", "dashboard/app.py"]
