
FROM stewartlab/scrnaseq_downstream3:v1

RUN rm -rf ./src/ ./data/ ./config.json
COPY src/ ./src/
COPY data/ ./data/
COPY environments/ ./environments/
COPY config.json ./

CMD ["/bin/bash"]