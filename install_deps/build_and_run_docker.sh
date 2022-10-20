#!/bin/bash

docker build -t install-deps .
docker run install-deps