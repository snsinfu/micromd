#!/usr/bin/env python3

import os.path
import sys


def main():
    target = None
    depends = dict()

    for line in sys.stdin:
        # We don't interpret the backslash continuation.
        line = line.replace("\\", "").strip()

        colon = line.find(":")
        if colon != -1:
            target = line[:colon]
            line = line[colon + 1:]

        sources = line.split()

        if target not in depends:
            depends[target] = []
        depends[target].extend(sources)

    for target in sorted(depends.keys()):
        sources = depends[target]
        sources = list(sorted(set(os.path.relpath(source) for source in sources)))

        print(f"{target}: \\")
        for i, source in enumerate(sources):
            source = os.path.relpath(source)
            if i == len(sources) - 1:
                end = ""
            else:
                end = " \\"
            print(f"  {source}{end}")


if __name__ == "__main__":
    main()
