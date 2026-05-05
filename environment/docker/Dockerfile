FROM anaconda/miniconda:latest

ENV CONDA_NO_PLUGINS=true

WORKDIR /app

COPY environment.yml /tmp/environment.yml

RUN conda config --set solver classic \
 && conda env create -f /tmp/environment.yml \
 && conda clean -afy \
 && rm /tmp/environment.yml

# adjust this to your env name
ENV PATH=/opt/conda/envs/fenics-lkm-env/bin:$PATH

COPY . /app

CMD ["python"]