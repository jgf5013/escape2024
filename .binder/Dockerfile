FROM rocker/binder:4.4.1
LABEL maintainer='Ka Yin Leung'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/rivm-syso/escape2024/raw/main/DESCRIPTION 
 && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2024-12-13/'))" \
 && R -e "devtools::install_deps()" \
 # install escape2024 package
 && R -e "install.packages('pak')" \
 && R -e "pak::pak('rivm-syso/escape2024')" \
 # install R packages needed for deliverable_2024 but not part of escape2024 package dependencies
 && R -e "pak::pak('ggplot2')" \
 && R -e "pak::pak('patchwork')" \
 && R -e "pak::pak('fs')" \
 && R -e "pak::pak('purrr')" \
 && R -e "pak::pak('readr')"

RUN rm DESCRIPTION.1; exit 0
