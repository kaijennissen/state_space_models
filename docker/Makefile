IMAGE = kalman:1.0
CONTAINER = kalman-container

THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))

build: 
	docker build --tag $(IMAGE)  \
	    -f $(THIS_FILE_DIR)/Dockerfile \
        $(THIS_FILE_DIR)

run:
	docker run --rm \
		-itd \
		-p 7772:22 \
	    -v $(THIS_FILE_DIR)/../:/opt/project \
	    --name $(CONTAINER) \
	    $(IMAGE)
