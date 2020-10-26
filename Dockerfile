FROM continuumio/miniconda:4.7.12

RUN conda install -y numpy=1.15.2 pandas=0.22
COPY concatenate_chain_results.py /scripts/
COPY collect_migmap_results_BCR.py /scripts/