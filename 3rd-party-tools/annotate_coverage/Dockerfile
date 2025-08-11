FROM hailgenetics/hail:0.2.133

# Set working directory
WORKDIR /opt

# Clone mtSwirl repository (assuming git is available in base image)
RUN git clone https://github.com/jamesemery/mtSwirl.git

# Upgrade pip and install only essential dependencies
RUN pip install --upgrade pip && \
    pip install --upgrade aiohttp

# Set the mtSwirl directory in PATH for easier access
ENV PATH="/opt/mtSwirl:${PATH}"

# Set working directory back to root for task execution
WORKDIR /
