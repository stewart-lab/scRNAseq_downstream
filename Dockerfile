
FROM stewartlab/scrnaseq_downstream3:v1

RUN rm -rf ./src/
COPY src/ ./src/
RUN ls -la ./src/

RUN rm -rf ./data/
COPY data/ ./data/
RUN ls -la ./data/

RUN rm -rf ./config.json
COPY config.json ./

CMD ["/bin/bash"]