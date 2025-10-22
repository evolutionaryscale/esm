# OSS-specific variables and commands
DOCKER_TAG ?= dev
DOCKER_IMAGE_OSS=oss_pytests:${DOCKER_TAG}

build-oss-ci:
	docker build \
		--output=type=docker \
		-f oss_pytests/Dockerfile \
		-t $(DOCKER_IMAGE_OSS) \
		oss_pytests

start-docker-oss:
	docker run \
		--rm \
		-e URL=${URL} \
		-e ESM3_FORGE_TOKEN=${ESM3_FORGE_TOKEN} \
		--name=$(USER)-oss_pytests \
		--network=host \
		${DOCKER_IMAGE_OSS}
