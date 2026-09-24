# -*- coding: utf-8 -*-
r"""
make_chapters.py -- turn the 17 standalone papers of papers/ into chapters of AG_monograph.tex.

For every paper it
  * keeps the body (between \begin{document} and \end{document}), drops \maketitle,
  * turns \title into \chapter[short]{full}\label{chap:chNN} and the abstract into a chapter summary,
  * prefixes every \label and every \ref / \eqref / \pageref with "chNN:" (19 labels collide),
  * turns citations of other papers of the book (keys BioN, CC26) into chapter references, and the
    words "Bio~N" in the text into "Chapter~\ref{chap:chMM}",
  * moves every other bibliography entry into AG_monograph.bib, merging duplicates,
  * points \includegraphics at ../papers/<folder>/,
and collects every paper's macros into AG_macros.tex, resolving conflicts.  The papers are only read.

    python make_chapters.py        (from monograph/)

Standard library only.  The chapter files it writes are generated; the manual harmonisation of the
chapters (step 4 of the roadmap) starts from them.
"""
import os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
PAPERS = os.path.join(HERE, "..", "papers")
OUT = os.path.join(HERE, "chapters")

# chapter, folder, tex file, short title (headers and contents), paper number (None: joint paper)
CH = [
    (1, "Ch01_Bio_01_Critical_Groups", "dbg_sandpile.tex", "Critical groups of de Bruijn graphs", 1),
    (2, "Ch02_Bio_02_Rotor_Routing", "rotor_assembly.tex", "The sandpile group as a coordinate system", 2),
    (3, "Ch03_Bio_03_Repeat_Splitting", "repeat_splitting.tex", "The sandpile group splits along the repeats", 3),
    (4, "Ch04_Bio_04_Multicopy_Repeats", "multicopy.tex", "Multi-copy repeats", 4),
    (5, "Ch05_Bio_06_RC_Double_Cover", "rc_double_cover.tex", "The reverse-complement double cover", 6),
    (6, "Ch06_Bio_07_Signed_Quotient", "rc_quotient_signed.tex", "The quotient is a signed graph", 7),
    (7, "Ch07_Bio_08_DS_Lattice_Sum", "anccft21.tex", "Double-stranded assembly as a lattice sum", 8),
    (8, "Ch08_Bio_12_Frontier_Elimination", "anccft25.tex", "Frontier elimination", 12),
    (9, "Ch09_Bio_13_Fermionic_Partition", "anccft26.tex", "A fermionic partition function", 13),
    (10, "Ch10_Bio_19_SharpP_Hardness", "anccft32.tex", "Double-stranded reconstruction is \\#P-hard", 19),
    (11, "Ch11_Bio_17_Inversions_Fibre", "anccft30.tex", "Inversions at inverted repeats", 17),
    (12, "Ch12_Spectral_Fibres", "spectral_fibres.tex", "What a local encoder cannot see", None),
    (13, "Ch13_Bio_09_Polyploid_Phasing", "anccft22.tex", "Polyploid phasing", 9),
    (14, "Ch14_Bio_10_Pangenome_Sheaf", "anccft23.tex", "The recombination group of a pan-genome", 10),
    (15, "Ch15_Bio_05_Ecoli_Decomposition", "ecoli_decomposition.tex", "The bundle-chain decomposition of E.~coli", 5),
    (16, "Ch16_Bio_14_Snarl_Schur", "anccft27.tex", "The BEST determinant along the snarl tree", 14),
    (17, "Ch17_Bio_11_Isoform_Polytope", "anccft24.tex", "The isoform decomposition polytope", 11),
]
BIO2CH = {bio: ch for ch, _, _, _, bio in CH if bio is not None}
JOINT_KEYS = {"CC26"}                        # the joint paper, chapter 12
RENAME = {(14, "Rec"): "RecGrp"}             # macro conflicts resolved by renaming in one chapter
report = []


# ------------------------------------------------------------------ small TeX parsing helpers

def braced(s, i):
    """s[i] == '{': return (content, index after the matching '}')"""
    assert s[i] == "{", s[i:i + 20]
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError("unbalanced braces")


def macros_of(preamble):
    """[(name, nargs, definition)] of \newcommand / \renewcommand in a preamble"""
    out = []
    for m in re.finditer(r"\\(?:re)?newcommand\*?\s*\{?\\([A-Za-z]+)\}?\s*(\[(\d)\])?\s*(\[[^\]]*\])?\s*(?=\{)",
                         preamble):
        body, _ = braced(preamble, m.end())
        out.append((m.group(1), m.group(3) or "", body.strip()))
    return out


def canon(defn):
    d = re.sub(r"\\mathrm\{([A-Za-z]+)\}", r"\\operatorname{\1}", defn)
    return re.sub(r"\s+", "", d)


def title_of(preamble):
    m = re.search(r"\\title\s*\{", preamble)
    t, _ = braced(preamble, m.end() - 1)
    t = re.sub(r"\\bfseries|\\large|\\Large", "", t)
    t = re.sub(r"\\\\(\[[^\]]*\])?", " ", t)
    return " ".join(t.split())


def bib_entries(body):
    m = re.search(r"\\begin\{thebibliography\}\{[^}]*\}(.*?)\\end\{thebibliography\}", body, re.S)
    if not m:
        return {}, body
    items = {}
    for it in re.finditer(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}(.*?)(?=\\bibitem|\Z)", m.group(1), re.S):
        items[it.group(1).strip()] = " ".join(it.group(2).split())
    return items, body[:m.start()] + body[m.end():]


def internal_chapter(key, text):
    """the chapter a bibliography entry refers to, if it is a paper of this book"""
    if key in JOINT_KEYS or "What a local encoder cannot see" in text:
        return 12
    m = re.fullmatch(r"Bio(\d+)", key)
    if m and int(m.group(1)) in BIO2CH:
        return BIO2CH[int(m.group(1))]
    return None


# ------------------------------------------------------------------ bibliography across chapters

GLOBAL = {}            # normalised text -> global key
BIBTEXT = {}           # global key -> text


def global_key(key, text, ch):
    norm = re.sub(r"[^a-z0-9]", "", re.sub(r"\\[a-zA-Z]+", "", text.lower()))[:90]
    if norm in GLOBAL:
        return GLOBAL[norm]
    g = key if key not in BIBTEXT else f"{key}ch{ch:02d}"
    GLOBAL[norm] = g
    BIBTEXT[g] = text
    return g


# ------------------------------------------------------------------ conversion

def chapter_list_ref(nums):
    chs = [BIO2CH[n] for n in nums]
    refs = [f"\\ref{{chap:ch{c:02d}}}" for c in chs]
    if len(refs) == 1:
        return f"Chapter~{refs[0]}"
    return "Chapters~" + ", ".join(refs[:-1]) + " and~" + refs[-1]


def convert(ch, folder, tex, short, bio, macro_table):
    src = open(os.path.join(PAPERS, folder, tex), encoding="utf-8").read()
    pre, rest = src.split("\\begin{document}", 1)
    body = rest.rsplit("\\end{document}", 1)[0]
    title = title_of(pre)
    for name, nargs, defn in macros_of(pre):
        macro_table[name].append((ch, nargs, defn))
    body = body.replace("\\maketitle", "")
    items, body = bib_entries(body)
    local = {}
    for k, t in items.items():
        c = internal_chapter(k, t)
        local[k] = ("chap", c) if c else ("bib", global_key(k, t, ch))

    def cite(m):
        opt, keys = m.group(1), [k.strip() for k in m.group(2).split(",")]
        chs, bibs = [], []
        for k in keys:
            if k not in local:
                report.append(f"ch{ch:02d}: \\cite key '{k}' has no bibliography entry")
                bibs.append(k)
            elif local[k][0] == "chap":
                chs.append(local[k][1])
            else:
                bibs.append(local[k][1])
        parts = []
        if chs:
            refs = [f"\\ref{{chap:ch{c:02d}}}" for c in chs]
            parts.append(("Chapter~" if len(refs) == 1 else "Chapters~") + ", ".join(refs))
        if bibs:
            parts.append(f"\\cite{{{','.join(bibs)}}}")
        s = " and ".join(parts)
        if opt:
            s += f", {opt[1:-1]}" if chs and not bibs else ""
            if bibs and not chs:
                s = f"\\cite{opt}{{{','.join(bibs)}}}"
        return s

    body = re.sub(r"\\cite(\[[^\]]*\])?\{([^}]*)\}", cite, body)

    # the words "Bio~8", "Bio 8", "Bio~12 and~13", "Bio~8, 12 and~13" in the text
    def words(m):
        nums = [int(x) for x in re.findall(r"\d+", m.group(0))]
        if not all(n in BIO2CH for n in nums):
            return m.group(0)
        return chapter_list_ref(nums)
    body = re.sub(r"Bio[~ ]\d+(?:(?:,\s*~?|\s*,?\s*and~?\s*)\d+)*(?![\d/])", words, body)

    # labels and references, prefixed by chapter
    def lab(m):
        cmd, star, key = m.group(1), m.group(2) or "", m.group(3)
        if key.startswith("chap:ch"):
            return m.group(0)
        return f"\\{cmd}{star}{{ch{ch:02d}:{key}}}"
    body = re.sub(r"\\(label|ref|eqref|pageref|autoref|nameref)(\*?)\{([^}]+)\}", lab, body)

    # figures live beside the paper
    body = re.sub(r"\\includegraphics(\[[^\]]*\])?\{([^}]+)\}",
                  lambda m: f"\\includegraphics{m.group(1) or ''}{{../papers/{folder}/{m.group(2)}}}", body)

    # macro renames for resolved conflicts
    for (c, old), new in RENAME.items():
        if c == ch:
            body = re.sub(r"\\" + old + r"(?![A-Za-z])", lambda _m: "\\" + new, body)

    # abstract -> chapter summary
    body = re.sub(r"\\begin\{abstract\}\s*(\\noindent)?", r"\\begin{chaptersummary}", body)
    body = body.replace("\\end{abstract}", "\\end{chaptersummary}")

    head = (f"% Chapter {ch}, generated by make_chapters.py from papers/{folder}/{tex}.\n"
            f"% Source paper: " + (f"Bio {bio}" if bio else "Z. Chen and R. Chen (joint)") + "\n")
    chap = f"\\chapter[{short}]{{{title}}}\\label{{chap:ch{ch:02d}}}\n"
    note = ("\\noindent{\\small\\emph{Joint work with Zhengyi Chen, published separately; included with the "
            "co-author's agreement.}}\\par\\medskip\n" if bio is None else
            f"\\noindent{{\\small\\emph{{Paper Bio~{bio} of the series. Verification script and transcript: "
            f"\\texttt{{papers/{folder.replace('_', chr(92) + '_')}/}}.}}}}\\par\\medskip\n")
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, f"ch{ch:02d}.tex"), "w", encoding="utf-8").write(head + chap + note + body.strip() + "\n")
    return title, len(items), sum(1 for v in local.values() if v[0] == "chap")


def main():
    macro_table = collections.defaultdict(list)
    print("chapter  bib entries (internal)  title")
    for ch, folder, tex, short, bio in CH:
        title, nb, ni = convert(ch, folder, tex, short, bio, macro_table)
        print(f"  ch{ch:02d}   {nb:3d} ({ni:2d})   {title[:80]}")
    # macros
    lines = ["% AG_macros.tex -- every chapter's macros, generated by make_chapters.py.",
             "% Conflicts are resolved below; see the comments.", "\\usepackage{tikz}",
             "\\usetikzlibrary{arrows.meta,positioning}"]
    for name in sorted(macro_table):
        defs = macro_table[name]
        forms = collections.OrderedDict()
        for ch, nargs, d in defs:
            forms.setdefault((nargs, canon(d)), []).append((ch, d))
        (nargs, _), users = next(iter(forms.items()))
        d = users[0][1]
        if "\\operatorname" in "".join(u[1] for f in forms.values() for u in f) and "\\mathrm{" in d \
                and len({k[1] for k in forms}) == 1:
            d = re.sub(r"\\mathrm\{([A-Za-z]+)\}", r"\\operatorname{\1}", d)
        lines.append(f"\\newcommand{{\\{name}}}" + (f"[{nargs}]" if nargs else "") + f"{{{d}}}")
        if len(forms) > 1:
            for (na, _), us in list(forms.items())[1:]:
                chs = [u[0] for u in us]
                new = next((RENAME[(c, name)] for c in chs if (c, name) in RENAME), None)
                if new:
                    lines.append(f"\\newcommand{{\\{new}}}" + (f"[{na}]" if na else "") +
                                 f"{{{us[0][1]}}}   % \\{name} of chapter(s) {chs}, renamed")
                    report.append(f"macro \\{name}: chapter(s) {chs} use \\{new}")
                else:
                    report.append(f"macro \\{name}: UNRESOLVED conflict in chapter(s) {chs}: {us[0][1]}")
    open(os.path.join(HERE, "AG_macros.tex"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    # bibliography
    with open(os.path.join(HERE, "AG_monograph.bib"), "w", encoding="utf-8") as f:
        for k, t in BIBTEXT.items():
            f.write(f"@misc{{{k},\n  key  = {{{k}}},\n  note = {{{t}}}\n}}\n\n")
    print(f"\nAG_macros.tex: {len(macro_table)} macros;  AG_monograph.bib: {len(BIBTEXT)} entries")
    for r in report:
        print("  note:", r)


if __name__ == "__main__":
    main()
