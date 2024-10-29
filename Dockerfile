
FROM stewartlab/scrnaseq_downstream:v1.1

RUN rm -rf ./src/
COPY src/ ./src/
RUN ls -la ./src/

RUN rm -rf ./environments/
COPY environments/ ./environments/
RUN ls -la ./environments/

RUN rm -rf ./data/
COPY data/ ./data/
RUN ls -la ./data/

RUN rm -rf ./config.json
COPY config.json ./

COPY run_downstream_toolkit.sh ./

CMD ["/bin/bash"]