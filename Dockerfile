from pantelispanka/biocpu
RUN R -e "install.packages(c('jsonlite', 'RCurl','Matrix', 'vegan'), repos='http://cran.cc.uoc.gr/mirrors/CRAN/')"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c('org.Hs.eg.db', 'GSEABase', 'GOstats', 'blockcluster'))'
COPY geoDescriptors.tar.gr /packages/
USER root
RUN R CMD INSTALL /packages/geoDescriptors.tar.gr --library=/usr/local/lib/R/site-library
CMD /usr/sbin/apache2ctl -D FOREGROUND