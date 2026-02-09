# Beagle 5.5 requires Java 8
# This Docker file is pushed to us.gcr.io/broad-gotc-prod/beagle5:0.1.1
FROM eclipse-temurin:8-jre-jammy

ARG BEAGLE_JAR_URL="https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"

# Install curl + certs, fetch the jar, and make a stable path
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/beagle \
    && curl -fsSL "${BEAGLE_JAR_URL}" -o /opt/beagle/beagle.jar

WORKDIR /work

# Default: run beagle; you can pass normal beagle args at `docker run ...`
ENTRYPOINT ["java","-jar","/opt/beagle/beagle.jar"]

