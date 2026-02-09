# Beagle 5.5 requires Java 8
# This Docker is pushed to pushed to us.gcr.io/broad-gotc-prod/beagle5:0.0.2
FROM eclipse-temurin:8-jre-jammy

ARG BEAGLE_JAR_URL="https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"

# Install prerequisites + Google Cloud CLI repo, then install google-cloud-cli
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl gnupg apt-transport-https \
    && mkdir -p /etc/apt/keyrings \
    && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /etc/apt/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/etc/apt/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# Fetch Beagle jar to a stable path
RUN mkdir -p /opt/beagle \
    && curl -fsSL "${BEAGLE_JAR_URL}" -o /opt/beagle/beagle.jar

WORKDIR /work

# IMPORTANT for WDL/Cromwell:
# Do NOT set ENTRYPOINT to java -jar, otherwise Cromwell's "bash -c ..." command won't run.
# Leave ENTRYPOINT unset and let the WDL command drive execution.
