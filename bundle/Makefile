ROOT = ../include/md/all.hpp
HEADERS = $(shell find ../include/md -name "*.hpp")
TEMPLATE = template.hpp
VERSION ?= unreleased

md.hpp: $(HEADERS) $(TEMPLATE)
	python3 bundle.py --template $(TEMPLATE) --version $(VERSION) $(ROOT) > $@
