"""Small helper for writing TOPAS parameter files.

TOPAS parameters have the form   <type>:<Name> = <value> [unit]
  s  string        d  dimensioned double (needs a unit)   u  unitless double
  i  integer       b  boolean
  sv/dv/uv/iv/bv   vectors: first token is the vector length
A value may also be an expression that references another parameter
(e.g. ``Ge/SAD - 10 cm``) -- pass such expressions as a :class:`Expr`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence


@dataclass(frozen=True)
class Expr:
    """A raw TOPAS expression that is written verbatim (no quoting, no unit added)."""

    text: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.text


def _num(x: float) -> str:
    """Format a number compactly but without losing precision."""
    if isinstance(x, Expr):
        return str(x)
    if isinstance(x, bool):
        raise TypeError("booleans must be written with ParamWriter.b()")
    if isinstance(x, int):
        return str(x)
    s = f"{float(x):.10g}"
    if "e" not in s and "." not in s:
        s += "."
    return s


def _quote(s: str) -> str:
    if '"' in s:
        raise ValueError(f"TOPAS string values cannot contain double quotes: {s!r}")
    return f'"{s}"'


@dataclass
class ParamWriter:
    """Accumulates TOPAS parameter lines and renders them as text."""

    lines: list[str] = field(default_factory=list)

    # -- structure -------------------------------------------------------
    def comment(self, text: str = "") -> "ParamWriter":
        for line in str(text).splitlines() or [""]:
            self.lines.append(f"# {line}".rstrip())
        return self

    def blank(self) -> "ParamWriter":
        self.lines.append("")
        return self

    def section(self, title: str) -> "ParamWriter":
        self.blank()
        self.lines.append("#" + "-" * 78)
        self.lines.append(f"# {title}")
        self.lines.append("#" + "-" * 78)
        return self

    def include(self, *paths: str) -> "ParamWriter":
        self.lines.append("includeFile = " + " ".join(paths))
        return self

    def raw(self, line: str) -> "ParamWriter":
        self.lines.append(line)
        return self

    # -- scalars ---------------------------------------------------------
    def s(self, name: str, value: str) -> "ParamWriter":
        self.lines.append(f"s:{name} = {_quote(value)}")
        return self

    def d(self, name: str, value, unit: str, changeable: bool = False) -> "ParamWriter":
        t = "dc" if changeable else "d"
        v = value if isinstance(value, Expr) else float(value)
        self.lines.append(f"{t}:{name} = {_num(v)} {unit}")
        return self

    def u(self, name: str, value) -> "ParamWriter":
        v = value if isinstance(value, Expr) else float(value)
        self.lines.append(f"u:{name} = {_num(v)}")
        return self

    def i(self, name: str, value) -> "ParamWriter":
        if not isinstance(value, Expr):
            if int(value) != value:
                raise ValueError(f"integer parameter {name} got non-integer {value}")
            value = int(value)
        self.lines.append(f"i:{name} = {value}")
        return self

    def b(self, name: str, value: bool) -> "ParamWriter":
        self.lines.append(f'b:{name} = "{"True" if value else "False"}"')
        return self

    # -- vectors ---------------------------------------------------------
    def sv(self, name: str, values: Sequence[str]) -> "ParamWriter":
        self.lines.append(f"sv:{name} = {len(values)} " + " ".join(_quote(v) for v in values))
        return self

    def dv(self, name: str, values: Iterable, unit: str) -> "ParamWriter":
        vals = [v if isinstance(v, Expr) else float(v) for v in values]
        self.lines.append(f"dv:{name} = {len(vals)} " + " ".join(_num(v) for v in vals) + f" {unit}")
        return self

    def uv(self, name: str, values: Iterable) -> "ParamWriter":
        vals = [v if isinstance(v, Expr) else float(v) for v in values]
        self.lines.append(f"uv:{name} = {len(vals)} " + " ".join(_num(v) for v in vals))
        return self

    def iv(self, name: str, values: Iterable) -> "ParamWriter":
        vals = [int(v) for v in values]
        self.lines.append(f"iv:{name} = {len(vals)} " + " ".join(str(v) for v in vals))
        return self

    # -- output ----------------------------------------------------------
    def extend(self, other: "ParamWriter") -> "ParamWriter":
        self.lines.extend(other.lines)
        return self

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"

    def write(self, path) -> None:
        from pathlib import Path

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(self.text())
