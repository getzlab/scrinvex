# Dockerfile for scR-Invex
FROM gcr.io/broad-cga-aarong-gtex/rnaseqc:4aaa2e2846e31cd3714c7b17f215b4556fb5121f
MAINTAINER Aaron Graubert

COPY Makefile /opt/scrinvex/
COPY src /opt/scrinvex/src
RUN cd /opt/scrinvex && ln -s /opt/rnaseqc /opt/scrinvex/rnaseqc && make && \
  ln -s /opt/scrinvex/scrinvex /usr/local/bin/scrinvex && make clean
