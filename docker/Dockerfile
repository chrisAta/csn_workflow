FROM python:3.7-slim-buster

RUN apt-get update && \
    apt-get install -y \
	     python3.7 \
	     python3-pip && \
    pip3 install \
        networkx \
        tqdm==4.31.1 \
        pandas\
        numpy \
        biopython

RUN mkdir /Work
RUN mkdir /Work/RRCoevNets

COPY ./lib/coev_net_creator.py /usr/bin/coev_net_creator.py
COPY ./lib/create_aln_net.py /usr/bin/create_aln_net.py
COPY ./lib/computeCoevSimilarity.py /usr/bin/computeCoevSimilarity.py
COPY ./lib/csnCreator.py /usr/bin/csnCreator.py

RUN chmod +x /usr/bin/coev_net_creator.py
RUN chmod +x /usr/bin/create_aln_net.py
RUN chmod +x /usr/bin/computeCoevSimilarity.py
RUN chmod +x /usr/bin/csnCreator.py
