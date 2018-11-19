import argparse
import os.path
import re


RE_LOCAL_INCLUDE = re.compile('#include "(.+?)"')
RE_VERSION_COMMENT = re.compile('// Version:\s*(.*)')
RE_BUNDLE_POINT = re.compile('// BUNDLE //')


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=str, default=None, help="Output template file")
    parser.add_argument("--version", type=str, default=None, help="Version string")
    parser.add_argument("root", type=str, help="Root header")
    return vars(parser.parse_args())


def run(root, template, version):
    root_header = Header(root)
    headers = load_headers(root_header.depend_paths)

    for header in headers.values():
        header.resolve_depends(headers)

    bundle = merge(headers.values())

    if template:
        print(render_template(template, bundle, version))
    else:
        print(bundle)


def load_headers(paths, headers=dict()):
    for path in paths:
        if path in headers:
            continue
        header = Header(path)
        headers[path] = header
        load_headers(header.depend_paths, headers)
    return headers


def render_template(template, bundle, version):
    with open(template) as file:
        lines = [line.rstrip() for line in file]

    output = []
    for line in lines:

        if RE_BUNDLE_POINT.match(line):
            output.append(bundle)
            continue

        if RE_VERSION_COMMENT.match(line) and version:
            ma = RE_VERSION_COMMENT.match(line)
            output.append(line[:ma.start(1)] + version)
            continue

        output.append(line)

    return "\n".join(output)


def merge(headers, resolved=set()):
    output = ""

    for header in headers:
        if header in resolved:
            continue

        if header.depends:
            output += merge(header.depends)

        output += header.make_output() + "\n"
        resolved.add(header)

    return output


class Header:
    def __init__(self, path):
        lines = load_code_lines(path)

        self._path = path
        self._code_lines = lines

        rel_depends = extract_local_includes(lines)
        basedir = os.path.dirname(path)
        self.depend_paths = [os.path.normpath(os.path.join(basedir, dep)) for dep in rel_depends]

    def resolve_depends(self, headers):
        self.depends = [headers[path] for path in self.depend_paths]

    def make_output(self):
        header = [
            "",
            "//------------------------------------------------------------------------------",
            "// " + strip_include_prefix(self._path),
            "//------------------------------------------------------------------------------",
            "",
        ]
        return "\n".join(header + remove_local_includes(self._code_lines))


def strip_include_prefix(path):
    prefix = "include/"
    pos = path.find(prefix)
    if pos == -1:
        return path
    pos += len(prefix)
    return path[pos:]


def remove_local_includes(lines):
    output = []
    for line in lines:
        if RE_LOCAL_INCLUDE.match(line):
            output = strip_empty_texts(output)
        else:
            output.append(line)
    return strip_empty_texts(output)


def extract_local_includes(lines):
    depends = []
    for line in lines:
        match = RE_LOCAL_INCLUDE.match(line)
        if match:
            depends.append(match.group(1))
    return depends


def load_code_lines(path):
    with open(path) as file:
        lines = [line.rstrip() for line in file]

    beg = index_if(lines, lambda s: s.startswith("#define")) + 1
    end = rindex_if(lines, lambda s: s.startswith("#endif")) - 1

    return strip_empty_texts(lines[beg:end])


def strip_empty_texts(texts):
    beg = index_if(texts, len) or 0
    end = rindex_if(texts, len) or len(texts)
    return texts[beg:end]


def index_if(texts, pred):
    for i, text in enumerate(texts):
        if pred(text):
            return i
    return None


def rindex_if(texts, pred):
    for i, text in enumerate(texts[::-1]):
        if pred(text):
            return len(texts) - i
    return None


if __name__ == "__main__":
    main()
