THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))
IMAGE = kalman:1.0
CONTAINER = kalman-container

build: 
	docker build --tag $(IMAGE)  \
	    -f Dockerfile \
        .

run:
	docker run --rm \
		-itd \
		-p 7772:22 \
	    -v $(THIS_FILE_DIR):/opt/project \
	    --name $(CONTAINER) \
	    $(IMAGE)
