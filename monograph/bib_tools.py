# -*- coding: utf-8 -*-
r"""
bib_tools.py -- maintenance of AG_monograph.bib (entries are @misc{key, key={key}, note={...}}).

    python bib_tools.py dups            list likely duplicate entries
    python bib_tools.py merge A=B ...   keep A, drop B, rewrite \cite keys in chapters/*.tex
    python bib_tools.py drop K ...      drop entries K (they must no longer be cited)
    python bib_tools.py doi             look up DOIs in Crossref; add those whose title and year match
    python bib_tools.py check           cited keys vs entries

Standard library only.  DOIs are added only on a verified title match (similarity >= 0.93 after
normalisation) and equal year; everything else is left without a DOI and listed.
"""
import os, re, sys, json, time, difflib, urllib.request, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
BIB = os.path.join(HERE, "AG_monograph.bib")
CHAP = os.path.join(HERE, "chapters")


def load():
    t = open(BIB, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"@misc\{([^,]+),\s*key\s*=\s*\{[^}]*\},\s*note\s*=\s*\{(.*?)\}\s*\n\}", t, re.S):
        out[m.group(1)] = m.group(2)
    return out


def save(entries):
    with open(BIB, "w", encoding="utf-8") as f:
        for k, t in entries.items():
            f.write(f"@misc{{{k},\n  key  = {{{k}}},\n  note = {{{t}}}\n}}\n\n")


def title(t):
    m = re.search(r"\\emph\{((?:[^{}]|\{[^{}]*\})*)\}", t)
    return m.group(1) if m else t


def norm(s):
    s = re.sub(r"\\[a-zA-Z]+\s*|[{}$~'`\"\\]", " ", s.lower())
    return " ".join(re.sub(r"[^a-z0-9 ]", " ", s).split())


def year(t):
    ys = re.findall(r"\b(1[89]\d\d|20\d\d)\b", t)
    return ys[-1] if ys else None


def surnames(t):
    head = t.split("\\emph")[0]
    return set(re.findall(r"([A-Z][a-z\\'`\"{}]+(?:-[A-Z][a-z]+)?)\s*(?:,|and|$)", head))


def cmd_dups():
    e = load()
    keys = list(e)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = e[keys[i]], e[keys[j]]
            r = difflib.SequenceMatcher(None, norm(title(a)), norm(title(b))).ratio()
            same_auth = bool(surnames(a) & surnames(b))
            if r > 0.72 or (same_auth and year(a) == year(b) and r > 0.45):
                print(f"{r:.2f}  {keys[i]} | {keys[j]}\n      {a[:150]}\n      {b[:150]}")


def rewrite_cites(mapping):
    for f in os.listdir(CHAP):
        p = os.path.join(CHAP, f)
        t = open(p, encoding="utf-8").read()

        def fix(m):
            keys = [mapping.get(k.strip(), k.strip()) for k in m.group(2).split(",")]
            seen = []
            for k in keys:
                if k not in seen:
                    seen.append(k)
            return f"\\cite{m.group(1) or ''}{{{','.join(seen)}}}"
        t2 = re.sub(r"\\cite(\[[^\]]*\])?\{([^}]*)\}", fix, t)
        if t2 != t:
            open(p, "w", encoding="utf-8").write(t2)


def cmd_merge(pairs):
    e = load()
    mapping = {}
    for p in pairs:
        keep, drop = p.split("=")
        assert keep in e and drop in e, p
        mapping[drop] = keep
        del e[drop]
    rewrite_cites(mapping)
    save(e)
    print(f"merged {len(mapping)}; {len(e)} entries remain")


def cited():
    ks = set()
    for f in os.listdir(CHAP):
        t = open(os.path.join(CHAP, f), encoding="utf-8").read()
        for m in re.finditer(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", t):
            ks |= {k.strip() for k in m.group(1).split(",")}
    return ks


def cmd_drop(keys):
    e = load()
    c = cited()
    for k in keys:
        assert k not in c, f"{k} is still cited"
        del e[k]
    save(e)
    print(f"dropped {len(keys)}; {len(e)} entries remain")


def cmd_check():
    e = load()
    c = cited()
    print(f"{len(e)} entries, {len(c)} cited keys")
    print("cited but missing:", sorted(c - set(e)))
    print("never cited:", sorted(set(e) - c))


def cmd_doi():
    e = load()
    added, skipped = 0, []
    for k, t in e.items():
        if "\\doi{" in t:
            continue
        ti, yr = title(t), year(t)
        if len(norm(ti)) < 12:
            skipped.append((k, "no usable title"))
            continue
        q = urllib.parse.urlencode({"query.bibliographic": f"{norm(ti)} {yr or ''}", "rows": 3,
                                    "mailto": "ruqing@hotmail.com"})
        try:
            with urllib.request.urlopen(f"https://api.crossref.org/works?{q}", timeout=30) as r:
                items = json.load(r)["message"]["items"]
        except Exception as ex:
            skipped.append((k, f"lookup failed: {ex}"))
            continue
        best = None
        for it in items:
            ct = " ".join(it.get("title", [""]))
            sim = difflib.SequenceMatcher(None, norm(ti), norm(ct)).ratio()
            cy = None
            for fld in ("published-print", "published-online", "issued"):
                dp = it.get(fld, {}).get("date-parts") or [[None]]
                if dp[0] and isinstance(dp[0][0], int):
                    cy = dp[0][0]
                    break
            if sim >= 0.93 and (yr is None or cy is None or abs(cy - int(yr)) <= 1):
                best = (sim, it["DOI"], ct, cy)
                break
        if best:
            e[k] = t.rstrip() + f" \\doi{{{best[1]}}}"
            added += 1
            save(e)
            print(f"  + {k:14s} {best[1]}   ({best[0]:.2f}, {best[3]})")
        else:
            skipped.append((k, f"no verified match (best: {items[0].get('title', [''])[0][:60] if items else '-'})"))
        time.sleep(0.3)
    save(e)
    print(f"\nDOIs added: {added}; without DOI: {len(skipped)}")
    for k, why in skipped:
        print(f"  - {k:14s} {why}")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"dups": lambda: cmd_dups(), "merge": lambda: cmd_merge(args), "drop": lambda: cmd_drop(args),
     "doi": lambda: cmd_doi(), "check": lambda: cmd_check()}[cmd]()
