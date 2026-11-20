#!/usr/bin/env python3
"""Open a Markdown file as a print-ready page in the browser.

By default the HTML is written only to a temporary file, the system browser
opens it with the print dialog, and that file is deleted. Nothing is left
next to the Markdown source.

Examples
--------
  python code/md2print.py code/PIPELINE.md
  python code/md2print.py CATALOGUE.md --pdf /tmp/catalogue.pdf
  python code/md2print.py README.md -o /tmp/readme-print.html
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PRINT_CSS = """
@page {
  size: A4;
  margin: 18mm 16mm;
}
:root {
  color-scheme: light;
  --fg: #111;
  --muted: #444;
  --border: #ccc;
  --code-bg: #f4f4f4;
  --link: #0645ad;
}
* { box-sizing: border-box; }
html { font-size: 11pt; }
body {
  margin: 0 auto;
  max-width: 48rem;
  padding: 1.25rem 1rem 2rem;
  color: var(--fg);
  font-family: "Source Serif 4", "Iowan Old Style", "Palatino Linotype",
               Palatino, "Book Antiqua", Georgia, serif;
  line-height: 1.45;
}
h1, h2, h3, h4 {
  font-family: "IBM Plex Sans", "Helvetica Neue", Helvetica, Arial, sans-serif;
  line-height: 1.25;
  page-break-after: avoid;
}
h1 { font-size: 1.7rem; margin: 0 0 0.8rem; }
h2 { font-size: 1.3rem; margin: 1.4rem 0 0.55rem; border-bottom: 1px solid var(--border); padding-bottom: 0.2rem; }
h3 { font-size: 1.1rem; margin: 1.1rem 0 0.4rem; }
h4 { font-size: 1rem; margin: 1rem 0 0.35rem; }
p, li { orphans: 3; widows: 3; }
a { color: var(--link); text-decoration: none; }
a:hover { text-decoration: underline; }
code, kbd, samp {
  font-family: "IBM Plex Mono", Menlo, Consolas, monospace;
  font-size: 0.88em;
  background: var(--code-bg);
  padding: 0.08em 0.28em;
  border-radius: 3px;
}
pre {
  background: var(--code-bg);
  border: 1px solid var(--border);
  border-radius: 4px;
  padding: 0.7rem 0.85rem;
  overflow-x: auto;
  font-size: 0.82rem;
  line-height: 1.35;
  page-break-inside: avoid;
}
pre code { background: none; padding: 0; }
blockquote {
  margin: 0.8rem 0;
  padding: 0.2rem 0.9rem;
  border-left: 3px solid var(--border);
  color: var(--muted);
}
table {
  width: 100%;
  border-collapse: collapse;
  margin: 0.8rem 0 1rem;
  font-size: 0.92rem;
  page-break-inside: avoid;
}
th, td {
  border: 1px solid var(--border);
  padding: 0.35rem 0.5rem;
  text-align: left;
  vertical-align: top;
}
th { background: #f0f0f0; font-weight: 600; }
tr:nth-child(even) td { background: #fafafa; }
img { max-width: 100%; height: auto; }
hr { border: 0; border-top: 1px solid var(--border); margin: 1.2rem 0; }
.meta {
  color: var(--muted);
  font-size: 0.85rem;
  margin-bottom: 1rem;
  font-family: "IBM Plex Sans", "Helvetica Neue", Helvetica, Arial, sans-serif;
}
@media print {
  body { max-width: none; padding: 0; }
  a { color: inherit; text-decoration: none; }
  a[href]::after { content: none; }
  pre, table, blockquote, img { break-inside: avoid; }
}
""".strip()


def which(cmd: str) -> str | None:
    return shutil.which(cmd)


def chrome_bin() -> str | None:
    mac = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
    if mac.is_file():
        return str(mac)
    for name in ("google-chrome", "chromium", "chromium-browser", "chrome"):
        found = which(name)
        if found:
            return found
    return None


def md_to_html(md_path: Path, html_path: Path, title: str | None = None) -> None:
    pandoc = which("pandoc")
    title = title or md_path.stem
    css_file = html_path.with_suffix(".print.css")
    css_file.write_text(PRINT_CSS, encoding="utf-8")

    if pandoc:
        cmd = [
            pandoc,
            str(md_path),
            "-f",
            "gfm+tex_math_dollars",
            "-t",
            "html5",
            "--standalone",
            "--metadata",
            f"title={title}",
            "--css",
            css_file.name,
            "-o",
            str(html_path),
        ]
        # Relative CSS href works when HTML and CSS sit side by side.
        subprocess.run(cmd, check=True, cwd=str(html_path.parent))
        # Prefer inlining CSS so the page is a single file for sharing/printing.
        html = html_path.read_text(encoding="utf-8")
        link = f'<link rel="stylesheet" href="{css_file.name}" />'
        style = f"<style>\n{PRINT_CSS}\n</style>"
        if link in html:
            html = html.replace(link, style)
        else:
            html = html.replace("</head>", f"{style}\n</head>")
        html_path.write_text(html, encoding="utf-8")
        css_file.unlink(missing_ok=True)
        return

    # Fallback: Python markdown (no pandoc).
    try:
        import markdown as mdlib
    except ImportError as exc:
        raise SystemExit(
            "Need pandoc or the Python package 'markdown'. "
            "Install one of them and retry."
        ) from exc

    body = mdlib.markdown(
        md_path.read_text(encoding="utf-8"),
        extensions=["extra", "tables", "fenced_code", "toc", "sane_lists", "smarty"],
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{title}</title>
<style>
{PRINT_CSS}
</style>
</head>
<body>
<p class="meta">Source: {md_path.name}</p>
{body}
</body>
</html>
"""
    html_path.write_text(html, encoding="utf-8")
    css_file.unlink(missing_ok=True)


def html_to_pdf(html_path: Path, pdf_path: Path) -> None:
    chrome = chrome_bin()
    if chrome:
        # Chrome needs absolute file:// URLs.
        url = html_path.resolve().as_uri()
        cmd = [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--no-pdf-header-footer",
            f"--print-to-pdf={pdf_path.resolve()}",
            url,
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        return

    pandoc = which("pandoc")
    if pandoc and (which("pdflatex") or which("xelatex") or which("lualatex")):
        engine = "xelatex" if which("xelatex") else ("lualatex" if which("lualatex") else "pdflatex")
        # Re-convert from the sibling .md if present; else from HTML.
        md_guess = html_path.with_suffix(".md")
        src = str(md_guess) if md_guess.is_file() else str(html_path)
        cmd = [pandoc, src, "-o", str(pdf_path), f"--pdf-engine={engine}"]
        subprocess.run(cmd, check=True)
        return

    raise SystemExit(
        "PDF export needs Google Chrome/Chromium or pandoc + a LaTeX engine."
    )


PRINT_ON_LOAD = """
<script>
window.addEventListener("load", function () {
  setTimeout(function () { window.print(); }, 250);
});
</script>
""".strip()


def prepare_browser_print(html_path: Path, md_path: Path) -> None:
    """Point relative links at the Markdown directory and open the print dialog."""
    html = html_path.read_text(encoding="utf-8")
    base = f'<base href="{md_path.parent.resolve().as_uri()}/">'
    if "<base " not in html:
        html = html.replace("<head>", f"<head>\n{base}", 1)
    if "window.print" not in html:
        html = html.replace("</body>", f"{PRINT_ON_LOAD}\n</body>", 1)
    html_path.write_text(html, encoding="utf-8")


def ephemeral_html_path() -> Path:
    handle = tempfile.NamedTemporaryFile(prefix="md2print-", suffix=".html", delete=False)
    path = Path(handle.name)
    handle.close()
    return path


def open_path(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.run(["open", str(path)], check=False)
    elif sys.platform.startswith("linux"):
        subprocess.run(["xdg-open", str(path)], check=False)
    elif sys.platform == "win32":
        os.startfile(str(path))  # type: ignore[attr-defined]


def release_ephemeral(path: Path, wait_s: float = 5.0) -> None:
    """Give the browser time to load the page, then delete the temp HTML."""
    time.sleep(wait_s)
    path.unlink(missing_ok=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Open Markdown in the browser for printing. "
            "The HTML is temporary unless -o is set."
        )
    )
    p.add_argument("markdown", type=Path, help="Input .md file")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Keep a copy of the HTML at this path (default: do not keep a file)",
    )
    p.add_argument(
        "--pdf",
        nargs="?",
        const="",
        default=None,
        metavar="FILE",
        help="Also write a PDF (default path: same stem as the HTML)",
    )
    p.add_argument("--title", help="Document title (default: markdown file stem)")
    p.add_argument("--open", action="store_true", help="Open the HTML (or PDF) when done")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    md_path = args.markdown.expanduser().resolve()
    if not md_path.is_file():
        print(f"error: not a file: {md_path}", file=sys.stderr)
        return 1
    if md_path.suffix.lower() not in {".md", ".markdown", ".txt"}:
        print(f"warning: expected a Markdown file, got {md_path.suffix!r}", file=sys.stderr)

    keep_html = args.output is not None
    html_path = (
        args.output.expanduser().resolve() if keep_html else ephemeral_html_path()
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    delete_html = not keep_html

    try:
        md_to_html(md_path, html_path, title=args.title)
        prepare_browser_print(html_path, md_path)

        pdf_path: Path | None = None
        if args.pdf is not None:
            pdf_path = (
                Path(args.pdf).expanduser().resolve()
                if args.pdf
                else md_path.with_suffix(".pdf")
            )
            html_to_pdf(html_path, pdf_path)
            print(f"PDF: {pdf_path}")

        # A plain run opens the print dialog. -o only keeps a file.
        # --open also opens a kept HTML, or the PDF when --pdf is set.
        open_browser = not keep_html or args.open
        if args.pdf is not None and not args.open:
            open_browser = False
        if open_browser:
            target = pdf_path if (args.open and pdf_path is not None) else html_path
            open_path(target)
            if delete_html and target == html_path:
                print("Opened in the browser for printing.")
                release_ephemeral(html_path)
                delete_html = False
        elif keep_html:
            print(f"HTML: {html_path}")
    finally:
        if delete_html:
            html_path.unlink(missing_ok=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
