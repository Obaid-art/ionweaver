# IonWeaver engine server. Built for Hugging Face Spaces (Docker SDK, port 7860)
# and for any other container host.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends libxrender1 libxext6 libsm6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY engine ./engine
COPY data ./data
COPY web ./web
COPY server ./server

# a non-root user, as Spaces expects
RUN useradd -m -u 1000 user && chown -R user /app
USER user

ENV PORT=7860
EXPOSE 7860
CMD ["sh", "-c", "uvicorn app:app --app-dir server --host 0.0.0.0 --port ${PORT}"]
