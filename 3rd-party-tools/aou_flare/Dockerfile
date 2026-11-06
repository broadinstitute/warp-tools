FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
 && add-apt-repository -y universe \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
    openjdk-17-jdk \
    git \
    wget \
    ca-certificates \
    tabix \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/flare

RUN git clone https://github.com/browning-lab/flare.git

RUN javac -cp flare/src/ flare/src/admix/AdmixMain.java \
 && jar cfe /opt/flare/flare.jar admix/AdmixMain -C flare/src/ ./

# Optional: download the prebuilt UW jar as well
RUN wget -O /opt/flare/flare.uw.jar https://faculty.washington.edu/browning/flare.jar

RUN printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -euo pipefail' \
  'exec java -jar /opt/flare/flare.jar "$@"' \
  > /usr/local/bin/flare \
 && chmod +x /usr/local/bin/flare

CMD ["flare", "--help"]
