# Dockerfile for scR-Invex
FROM gcr.io/broad-cga-aarong-gtex/rnaseqc:5fe48d8c75cd9acc0a8cd186d12afd628926744f
MAINTAINER Aaron Graubert

COPY Makefile /opt/scrinvex/
COPY src /opt/scrinvex/src
RUN cd /opt/scrinvex && ln -s /opt/rnaseqc /opt/scrinvex/rnaseqc && make && \
  ln -s /opt/scrinvex/scrinvex /usr/local/bin/scrinvex && make clean
