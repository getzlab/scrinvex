# Dockerfile for scR-Invex
FROM gcr.io/broad-cga-aarong-gtex/rnaseqc:98aebe7e97487669ab9bf9b1368f591a201caeaf
MAINTAINER Aaron Graubert

COPY Makefile /opt/scrinvex/
COPY src /opt/scrinvex/src
RUN cd /opt/scrinvex && ln -s /opt/rnaseqc /opt/scrinvex/rnaseqc && make && \
  ln -s /opt/scrinvex/scrinvex /usr/local/bin/scrinvex && make clean
