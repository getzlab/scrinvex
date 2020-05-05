# Dockerfile for scR-Invex
FROM gcr.io/broad-cga-aarong-gtex/rnaseqc:495d7194f92c71eb8c58d7660e4f6acd1801832a
MAINTAINER Aaron Graubert

COPY Makefile /opt/scrinvex/
COPY src /opt/scrinvex/src
RUN cd /opt/scrinvex && ln -s /opt/rnaseqc /opt/scrinvex/rnaseqc && make && \
  ln -s /opt/scrinvex/scrinvex /usr/local/bin/scrinvex && make clean
