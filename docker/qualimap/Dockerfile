FROM nfcore/base:1.9
LABEL authors="Uwe Schwartz" \
      description="container for qualimap"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/qualimap/bin:$PATH


# Instruct R processes to use these empty files instead of clashing with a local version
#RUN touch .Rprofile
#RUN touch .Renviron
