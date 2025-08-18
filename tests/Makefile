# OSS-specific variables and commands
DOCKER_TAG ?= dev
DOCKER_IMAGE_OSS=oss_pytests:${DOCKER_TAG}

build-oss-ci:
	docker build -f oss_pytests/Dockerfile oss_pytests -t $(DOCKER_IMAGE_OSS)

start-docker-oss:
	docker run \
		--rm \
		-e URL=${URL} \
		-e ESM3_FORGE_TOKEN=${ESM3_FORGE_TOKEN} \
		--name=$(USER)-oss_pytests \
		--network=host \
		${DOCKER_IMAGE_OSS}
