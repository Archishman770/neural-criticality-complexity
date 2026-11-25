"""Build the submission .docx from docs/entropy/manuscript.md (python-docx).

Markdown subset: header fields (TITLE/AUTHOR/AFFILIATION/CORRESPONDENCE/ABSTRACT/KEYWORDS),
'# ' / '## ' headings, '- ' and '  - ' bullets, pipe tables, 'FIGURE: path | caption',
'BACKMATTER:', 'REFERENCES' followed by 'key | formatted reference' lines, inline **bold**,
*italic*, and citations [@a; @b]. Citations are numbered in order of first appearance.
The build fails if a cited key has no reference or a reference is never cited.

Usage:  python3 tools/build_manuscript_docx.py [src.md] [out.docx]
"""
from __future__ import annotations

import re
import sys

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

SRC = sys.argv[1] if len(sys.argv) > 1 else "docs/entropy/manuscript.md"
OUT = sys.argv[2] if len(sys.argv) > 2 else "docs/entropy/Vemula_Entropy_manuscript.docx"
CITE = re.compile(r"\[(@[\w]+(?:;\s*@[\w]+)*)\]")
INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*)")


def parse(path: str) -> tuple[dict, list[str], dict]:
    """Split the source into header fields, body lines and the reference map."""
    head, body, refs = {}, [], {}
    in_refs = False
    for ln in open(path, encoding="utf-8").read().split("\n"):
        if ln.strip() == "REFERENCES":
            in_refs = True
            continue
        if in_refs:
            if "|" in ln:
                k, v = ln.split("|", 1)
                refs[k.strip()] = v.strip()
            continue
        m = re.match(r"^(TITLE|AUTHOR|AFFILIATION|CORRESPONDENCE|ABSTRACT|KEYWORDS):\s*(.*)$", ln)
        in_header = not any(b.startswith("# ") for b in body)
        if m and in_header:
            head[m.group(1)] = m.group(2)
        elif in_header and not ln.strip():
            continue
        else:
            body.append(ln)
    return head, body, refs


def number_citations(texts: list[str]) -> list[str]:
    """Citation keys in order of first appearance."""
    order: list[str] = []
    for t in texts:
        for grp in CITE.findall(t):
            for k in re.findall(r"@([\w]+)", grp):
                if k not in order:
                    order.append(k)
    return order


def fmt_nums(nums: list[int]) -> str:
    """[1,2] for pairs, [1–3] for runs of three or more."""
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ",".join(map(str, nums[i:j + 1])))
        i = j + 1
    return "[" + ",".join(out) + "]"


def add_runs(p, text: str, size: float | None = None) -> None:
    """Add text with **bold** and *italic* spans to paragraph p."""
    for part in INLINE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            r = p.add_run(part[2:-2])
            r.bold = True
        elif part.startswith("*") and part.endswith("*"):
            r = p.add_run(part[1:-1])
            r.italic = True
        else:
            r = p.add_run(part)
        if size:
            r.font.size = Pt(size)


def line_numbers(doc: Document) -> None:
    """Continuous line numbering for review."""
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    sect = doc.sections[0]._sectPr
    pg_mar = sect.find(qn("w:pgMar"))      # schema order: ... pgMar, paperSrc, pgBorders, lnNumType ...
    pg_mar.addnext(ln)
    zoom = doc.settings.element.find(qn("w:zoom"))   # python-docx template omits required w:percent
    if zoom is not None and zoom.get(qn("w:percent")) is None:
        zoom.set(qn("w:percent"), "100")


def main() -> None:
    head, body, refs = parse(SRC)
    order = number_citations([head.get("ABSTRACT", "")] + body)
    missing = [k for k in order if k not in refs]
    unused = [k for k in refs if k not in order]
    if missing or unused:
        sys.exit(f"citation check failed: missing={missing} unused={unused}")
    num = {k: i + 1 for i, k in enumerate(order)}

    def cite(t: str) -> str:
        return CITE.sub(lambda m: fmt_nums([num[k] for k in re.findall(r"@([\w]+)", m.group(1))]), t)

    print(f"abstract words: {len(head['ABSTRACT'].split())} (Entropy: about 200 maximum)")

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Inches(1))
    line_numbers(doc)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(head["TITLE"])
    r.bold = True
    r.font.size = Pt(16)
    for key in ("AUTHOR", "AFFILIATION", "CORRESPONDENCE"):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        add_runs(p, head[key])
    add_runs(doc.add_paragraph(), "**Abstract:** " + cite(head["ABSTRACT"]))
    add_runs(doc.add_paragraph(), "**Keywords:** " + head["KEYWORDS"])

    i = 0
    while i < len(body):
        ln = body[i]
        if ln.startswith("# "):
            doc.add_heading(ln[2:], level=1)
        elif ln.startswith("## "):
            doc.add_heading(ln[3:], level=2)
        elif ln.startswith("BACKMATTER:"):
            doc.add_paragraph()
        elif ln.startswith("FIGURE:"):
            path, cap = [s.strip() for s in ln[len("FIGURE:"):].split("|", 1)]
            doc.add_picture(path, width=Inches(6.5))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(doc.add_paragraph(), cite(cap), size=10)
        elif ln.startswith("|"):
            rows = []
            while i < len(body) and body[i].startswith("|"):
                cells = [c.strip() for c in body[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r"-+", c) for c in cells):
                    rows.append(cells)
                i += 1
            t = doc.add_table(rows=len(rows), cols=len(rows[0]))
            t.style = "Table Grid"
            for ri, row in enumerate(rows):
                for ci, c in enumerate(row):
                    add_runs(t.cell(ri, ci).paragraphs[0], ("**" + c + "**") if ri == 0 else cite(c), size=9)
            continue
        elif ln.startswith("  - "):
            add_runs(doc.add_paragraph(style="List Bullet 2"), cite(ln[4:]))
        elif ln.startswith("- "):
            add_runs(doc.add_paragraph(style="List Bullet"), cite(ln[2:]))
        elif ln.strip():
            add_runs(doc.add_paragraph(), cite(ln))
        i += 1

    doc.add_heading("References", level=1)
    for k in order:
        add_runs(doc.add_paragraph(), f"{num[k]}. " + refs[k], size=10)
    doc.save(OUT)
    print(f"wrote {OUT}: {len(order)} references, all cited")


if __name__ == "__main__":
    main()
