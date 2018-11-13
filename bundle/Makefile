ROOT = ../include/md/all.hpp
HEADERS = $(shell find ../include/md -name "*.hpp")
TEMPLATE = md_template.hpp

md.hpp: $(HEADERS) $(TEMPLATE)
	python3 bundle.py --template $(TEMPLATE) $(ROOT) > $@
