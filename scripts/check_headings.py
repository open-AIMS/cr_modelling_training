#!/usr/bin/env python3
"""Reports any heading in a module that Pandoc will not read as a heading.

Pandoc's markdown needs a blank line above an ATX heading. Without one the
`## Title` line renders as literal text inside the paragraph above it, so the
section vanishes from the page and from the sidebar while every chunk inside it
still runs and the render reports success. A scripted replacement of module 5's
Poisson section dropped that blank line on 2026-09-22 and the section was lost
from `docs/` with nothing in the log to say so.

Run after any edit that moves or replaces a section:

    python3 scripts/check_headings.py

Exits 1 where anything is reported. The callout titles in
`0Software-setup.qmd` sit directly under a `::: {.callout-*}` opener, which is
the documented Quarto idiom for a callout title, and are not reported.
"""

import glob
import pathlib
import re
import sys

HEADING = re.compile(r"^#{1,4} \S")
CALLOUT = re.compile(r"^:::+\s*\{\.callout")


def offenders(path):
    lines = pathlib.Path(path).read_text(encoding="utf-8").split("\n")
    in_fence = False
    out = []
    for i, line in enumerate(lines):
        if line.startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence or not HEADING.match(line) or i == 0:
            continue
        above = lines[i - 1].strip()
        if above == "" or CALLOUT.match(above):
            continue
        out.append((i + 1, line))
    return out


def main(paths):
    paths = paths or sorted(glob.glob("vignettes/*.qmd"))
    found = 0
    for path in paths:
        for line_no, text in offenders(path):
            print(f"{path}:{line_no}: no blank line above -- {text}")
            found += 1
    if found:
        print(f"\n{found} heading(s) Pandoc will render as literal text.")
        return 1
    print(f"{len(paths)} file(s) checked, every heading has a blank line above it.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
