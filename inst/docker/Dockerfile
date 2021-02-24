FROM alpine:3.13

RUN apk add --no-cache \
        build-base \
        zlib-dev \
        bzip2-dev \
        xz-dev \
        curl-dev \
        openjdk7-jre \
        bash

## samtools
## ADD rather than RUN curl because macOS is currently broken wrt github.com
ADD https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 /tmp

RUN cd /tmp && \
        tar xjf samtools-1.11.tar.bz2 && \
        cd samtools-1.11 && \
        ./configure --without-curses --prefix=/usr/local && \
        make -j && \
        make install && \
        rm -rf /tmp/samtools*

## bwa
ADD https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 /tmp

RUN cd /tmp && \
        tar xjf bwa-0.7.17.tar.bz2 && \
        cd bwa-0.7.17 && \
        sed -i 's/const uint8_t rle_auxtab/extern const uint8_t rle_auxtab/' rle.h && \
        make -j && \
        ln -s $PWD/bwa /usr/local/bin

## HMMRATAC
ADD https://github.com/LiuLabUB/HMMRATAC/releases/download/1.2.10/HMMRATAC_V1.2.10_exe.jar /
