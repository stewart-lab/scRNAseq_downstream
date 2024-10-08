
FROM stewartlab/scrnaseq_downstream:v1.1

COPY src/ ./src/
RUN ls -la ./src/

COPY environments/ ./environments/
RUN ls -la ./environments/

COPY data/ ./data/
RUN ls -la ./data/

COPY config.json ./
RUN ls -la

CMD ["/bin/bash"]