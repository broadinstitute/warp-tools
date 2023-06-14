FROM alpine:3.14

ARG SAMTOOLS_VERSION=1.6 \
    PICARD_PUBLIC_VERSION=2.26.10

ENV TERM=xterm-256color \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    PICARD_URL=https://github.com/broadinstitute/picard/releases/download/${PICARD_PUBLIC_VERSION}/picard.jar

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        SAMTOOLS_VERSION=${SAMTOOLS_VERSION} \
	PICARD_PUBLIC_VERSION=2.26.10

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
	apk add --no-cache \
            autoconf \
            automake \ 
            bash \
            bzip2-dev \
            curl-dev \
            gcc \
            libressl-dev \
            make \ 
            musl-dev \
            ncurses-dev \
            perl \
            tini \ 
            wget \ 
            xz-dev \
            zlib-dev \
	    gcompat \
            libc6-compat \
            openjdk8

# Install samtools and picard

 RUN  mkdir temp; \
	cd temp; \
	wget ${SAMTOOLS_URL}; \
	tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
	cd samtools-${SAMTOOLS_VERSION}; \
	./configure; \
	make; \
	make install; \
	rm -r /usr/gitc/temp
	
Run wget ${PICARD_URL}


# install Drop-seq tools
RUN wget https://github.com/broadinstitute/Drop-seq/releases/download/v2.5.3/Drop-seq_tools-2.5.3.zip; \
	apk add unzip; \
	unzip Drop-seq_tools-2.5.3.zip; \
	rm Drop-seq_tools-2.5.3.zip; \
	echo 'export PATH=$/usr/gitc/Drop-seq_tools-2.5.3:$PATH' >> ~/.bashrc; \
	source ~/.bashrc

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
