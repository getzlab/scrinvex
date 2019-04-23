# Dockerfile for scR-Invex
FROM gcr.io/broad-cga-aarong-gtex/rnaseqc:60c06f2835ffdc8ffd29ab059fd6b2de923e0972
MAINTAINER Aaron Graubert

COPY Makefile /opt/scrinvex/
COPY src /opt/scrinvex/src
RUN cd /opt/scrinvex && ln -s /opt/rnaseqc /opt/scrinvex/rnaseqc && make && \
  ln -s /opt/scrinvex/scrinvex /usr/local/bin/scrinvex && make clean
