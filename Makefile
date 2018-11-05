ROOT = include/md/all.hpp
HEADERS = $(shell find include/md -name "*.hpp")
TEMPLATE = bundle/md_template.hpp

include/md.hpp: $(HEADERS) $(TEMPLATE)
	python3 bundle/bundle.py --template $(TEMPLATE) $(ROOT) > $@
