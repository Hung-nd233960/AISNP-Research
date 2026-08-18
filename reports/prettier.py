"""Normalize a Markdown report before Pandoc conversion to Word/LaTeX.

Two transforms, applied line-by-line, skipping fenced code blocks:
  1. Drop standalone horizontal-rule lines (a line consisting only of `---`).
  2. Replace em dashes (—) with a plain hyphen (-).

Non-destructive by default: writes to <input>.clean.<ext> unless --in-place
or --output is given, so the hand-authored source keeps its own formatting
and only the pipeline's Pandoc-facing copy is normalized.

Usage:
    python reports/prettier.py reports/REPORT.md
    python reports/prettier.py reports/REPORT.md --in-place
    python reports/prettier.py reports/REPORT.md -o reports/REPORT.pandoc.md
    python reports/prettier.py reports/REPORT.md --check
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

HR_LINE = re.compile(r"^\s{0,3}-{3,}\s*$")
FENCE_LINE = re.compile(r"^\s*```")
EM_DASH = "—"
BLANK_RUN = re.compile(r"\n{3,}")


def clean(text: str) -> str:
    in_fence = False
    out_lines: list[str] = []
    for line in text.split("\n"):
        if FENCE_LINE.match(line):
            in_fence = not in_fence
            out_lines.append(line)
            continue
        if in_fence:
            out_lines.append(line)
            continue
        if HR_LINE.match(line):
            continue
        out_lines.append(line.replace(EM_DASH, "-"))
    return BLANK_RUN.sub("\n\n", "\n".join(out_lines))


def default_output_path(input_path: Path) -> Path:
    return input_path.with_suffix(f".clean{input_path.suffix}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path, help="Markdown file to clean")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Explicit output path")
    parser.add_argument("--in-place", action="store_true", help="Overwrite the input file directly")
    parser.add_argument("--check", action="store_true", help="Exit 1 if the file would change; write nothing")
    args = parser.parse_args()

    original = args.input.read_text(encoding="utf-8")
    cleaned = clean(original)
    changed = cleaned != original

    if args.check:
        if changed:
            print(f"{args.input}: would be modified", file=sys.stderr)
            return 1
        print(f"{args.input}: already clean")
        return 0

    if args.in_place:
        target = args.input
    else:
        target = args.output or default_output_path(args.input)

    target.write_text(cleaned, encoding="utf-8")
    print(f"Wrote {target}" if changed else f"Wrote {target} (no changes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
