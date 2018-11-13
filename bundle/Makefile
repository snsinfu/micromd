ROOT = ../include/md/all.hpp
HEADERS = $(shell find ../include/md -name "*.hpp")
TEMPLATE = template.hpp

md.hpp: $(HEADERS) $(TEMPLATE)
	python3 bundle.py --template $(TEMPLATE) $(ROOT) > $@
