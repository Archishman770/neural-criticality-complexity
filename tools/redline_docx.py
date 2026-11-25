"""Tracked-change (redline) helpers for an unpacked .docx word/document.xml.

Not part of the analysis pipeline: used to apply manuscript corrections as Word
tracked changes so the author can accept or reject each one. Works on documents that
already contain tracked changes: matching is done on the accepted view (inserted text
counts, deleted text does not), and only runs overlapping the target text are split.

Usage (from a script):
    from tools.redline_docx import Redline
    r = Redline("unpacked/word/document.xml", author="Claude", date="2026-10-03T00:00:00Z")
    r.replace("old sentence.", "new sentence.")
    r.save()
"""
from __future__ import annotations

import html
import re

P_RE = re.compile(r"<w:p[ >].*?</w:p>", re.S)
# A run: optional rPr, then content up to the closing tag.
RUN_RE = re.compile(r"<w:r>(<w:rPr>(?:(?!</w:rPr>).)*</w:rPr>)?((?:(?!</w:r>).)*)</w:r>", re.S)


def _esc(s: str) -> str:
    return html.escape(s, quote=False)


def _strip_deleted(p: str) -> str:
    return re.sub(r"<w:del [^>]*>.*?</w:del>", "", p, flags=re.S)


def accepted_text(p: str) -> str:
    """Paragraph text with tracked deletions removed and insertions kept."""
    q = _strip_deleted(p)
    return html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>", q, re.S)))


class Redline:
    """Apply tracked replacements to one document.xml."""

    def __init__(self, path: str, author: str, date: str, first_id: int = 20000) -> None:
        self.path = path
        self.auth = f'w:author="{author}" w:date="{date}"'
        self.x = open(path, encoding="utf-8").read()
        used = [int(i) for i in re.findall(r'w:id="(\d+)"', self.x)]
        self.next_id = max([first_id] + [u + 1 for u in used])

    def _nid(self) -> int:
        self.next_id += 1
        return self.next_id

    def _run(self, t: str, rpr: str) -> str:
        return f'<w:r>{rpr}<w:t xml:space="preserve">{_esc(t)}</w:t></w:r>' if t else ""

    def _ins(self, t: str, rpr: str) -> str:
        return f'<w:ins w:id="{self._nid()}" {self.auth}>{self._run(t, rpr)}</w:ins>' if t else ""

    def _del(self, t: str, rpr: str) -> str:
        if not t:
            return ""
        return (f'<w:del w:id="{self._nid()}" {self.auth}><w:r>{rpr}'
                f'<w:delText xml:space="preserve">{_esc(t)}</w:delText></w:r></w:del>')

    def replace(self, old: str, new: str, expect: int = 1) -> None:
        """Tracked replacement of `old` (accepted-view text) by `new` in each matching paragraph."""
        hits = [m for m in P_RE.finditer(self.x) if old in accepted_text(m.group(0))]
        if len(hits) != expect:
            raise SystemExit(f"{len(hits)} hits (expected {expect}) for: {old[:80]}")
        for m in reversed(hits):
            p = m.group(0)
            # Candidate runs: runs not inside a w:del element.
            spans = []
            for r in RUN_RE.finditer(p):
                before = p[: r.start()]
                inside_del = before.rfind("<w:del ") > before.rfind("</w:del>")
                if inside_del or "<w:delText" in r.group(0):
                    continue
                t = html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>",
                                                     r.group(0), re.S)))
                spans.append((r, t))
            full = "".join(t for _, t in spans)
            s0 = full.index(old)
            s1 = s0 + len(old)
            offsets, acc = [], 0
            for _, t in spans:
                offsets.append(acc)
                acc += len(t)
            for (r, t), a0 in reversed(list(zip(spans, offsets))):
                a1 = a0 + len(t)
                if a1 <= s0 or a0 >= s1:
                    continue
                rpr = r.group(1) or ""
                lo, hi = max(s0, a0) - a0, min(s1, a1) - a0
                tail = self._ins(new, rpr) if a0 < s1 <= a1 else ""
                repl = (self._run(t[:lo], rpr) + self._del(t[lo:hi], rpr) + tail
                        + self._run(t[hi:], rpr))
                p = p[: r.start()] + repl + p[r.end():]
            self.x = self.x[: m.start()] + p + self.x[m.end():]

    def save(self) -> None:
        open(self.path, "w", encoding="utf-8").write(self.x)
