FROM biocontainers/danpos:v2.2.2_cv3
LABEL authors="Uwe Schwartz" \
      description="Docker image containing DANPOS software requirements for the nucMACC pipeline"

# Install the wigtobigwig

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/wigToBigWig -P /opt/wig2bigwig \
&& chmod +x /opt/wig2bigwig/wigToBigWig

# Add wigtobigwig dir to PATH
ENV PATH=/opt/wig2bigwig:$PATH
