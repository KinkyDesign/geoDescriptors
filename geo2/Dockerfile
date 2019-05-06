from pantelispanka/biocpu:latest
RUN R -e "install.packages(c('jsonlite', 'RCurl','Matrix', 'vegan'), repos='http://cran.cc.uoc.gr/mirrors/CRAN/')"

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c("org.Hs.eg.db", "GSEABase", "GOstats", "blockcluster", "Category", "GO.db"))'
COPY geoDescriptors.tar.gz /packages/
COPY server.conf /etc/opencpu/

USER root
RUN R CMD INSTALL /packages/geoDescriptors.tar.gz --library=/usr/local/lib/R/site-library
CMD /usr/sbin/apache2ctl -D FOREGROUND
