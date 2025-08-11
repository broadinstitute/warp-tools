FROM hailgenetics/hail:0.2.133

# Set working directory
WORKDIR /opt

# Clone mtSwirl repository (assuming git is available in base image)
RUN git clone https://github.com/jamesemery/mtSwirl.git

# Upgrade pip and install only essential dependencies
RUN pip install --upgrade pip==25.2 && \
    pip install --upgrade aiohttp==3.12.15

# Set the mtSwirl directory in PATH for easier access
ENV PATH="/opt/mtSwirl:${PATH}"

