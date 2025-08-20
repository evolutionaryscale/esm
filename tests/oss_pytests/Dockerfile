# Dockerfile.sdktest
FROM python:3.12-slim

# Install pip and basic dependencies
RUN apt-get update && apt-get install -y curl build-essential && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /tests

# Copy requirements and install them (assumes esm-oss is one of them)
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy in the tests
COPY . .

# Default command (can be overridden in docker run)
CMD ["pytest", "-v", "test_oss_client.py"]
