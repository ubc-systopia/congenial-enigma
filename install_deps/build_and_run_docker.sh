#!/bin/bash

docker build \
    -t install-deps \
    .
    # --build-arg arg=2.3 \
    # --build-arg pip_reqs=${}
    # .
docker run install-deps