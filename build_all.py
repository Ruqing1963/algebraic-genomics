# -*- coding: utf-8 -*-
r"""
build_all.py -- one-click build and verification of the Algebraic Genomics repository
(adapted from AG_build.py to the repository layout: papers/<chapter>/, monograph/, lib/).

For every chapter: compile the paper with pdflatex (twice), run its verification script(s), and
compare the result with the saved reference transcript.  Writes AG_build_report.txt and the new
transcripts under _build/logs/.  Reference transcripts are never modified: every .txt file of a
chapter folder is snapshotted before its scripts run and restored afterwards, and checkpoints that
would let a run skip its work are moved aside for the duration.

    python -u build_all.py --quick        papers + fast checks            (about 30-40 minutes)
    python -u build_all.py --full         papers + every computation      (about 4 hours)
    options:  --only "Bio_13,Bio_17"     restrict to some chapters
              --no-tex / --no-scripts    skip one half
              --inventory                list what would be run, run nothing

A script PASSES when it exits with code 0 and its last line of the form "TOTAL a/b passed"
(or "a/b passed", or "ALL CHECKS PASSED") has a == b.  A paper PASSES when pdflatex reports no
error and no undefined reference; overfull boxes are reported as warnings.  In --full mode the
transcript is also compared with the reference after removing timings; differences are reported,
not failed, because some references were edited by hand (e.g. omitted progress lines).
Standard library only.
"""
import os, sys, re, time, shutil, subprocess, difflib, argparse

ROOT = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(ROOT, "_build")
LOGS = os.path.join(BUILD, "logs")

# (chapter, folder, tex or None, [scripts]); a script is a dict:
#   file, ref (reference transcript), quick (args, or None = skipped in quick), full (args),
#   timeout (s, full mode), hide (files moved aside while running), out_flag (write to a file)
H = 3600
CHAPTERS = [
    ("1", "papers/Bio_01_Critical_Groups",  "dbg_sandpile.tex",       [dict(file="dbg_sandpile.py", ref="dbg_sandpile_output.txt", quick=[], full=[], timeout=H)]),
    ("2", "papers/Bio_02_Rotor_Routing",  "rotor_assembly.tex",     [dict(file="rotor_assembly.py", ref="rotor_assembly_output.txt", quick=[], full=[], timeout=H)]),
    ("3", "papers/Bio_03_Repeat_Splitting",  "repeat_splitting.tex",   [dict(file="repeat_splitting.py", ref="repeat_splitting_output.txt", quick=[], full=[], timeout=H)]),
    ("4", "papers/Bio_04_Multicopy_Repeats",  "multicopy.tex",          [dict(file="multicopy.py", ref="multicopy_output.txt", quick=[], full=[], timeout=H)]),
    ("15", "papers/Bio_05_Ecoli_Decomposition",  "ecoli_decomposition.tex", [dict(file="ecoli_sandpile.py", ref="ecoli_sandpile_output.txt", quick=[], full=[], timeout=H)]),
    ("5", "papers/Bio_06_RC_Double_Cover",  "rc_double_cover.tex",    [dict(file="rc_double_cover.py", ref="rc_double_cover_output.txt", quick=[], full=[], timeout=H)]),
    ("6", "papers/Bio_07_Signed_Quotient",  "rc_quotient_signed.tex", [dict(file="bidirected_sandpile.py", ref="bidirected_sandpile_output.txt", quick=None, full=[], timeout=2 * H)]),
    ("7", "papers/Bio_08_DS_Lattice_Sum",  "anccft21.tex",           [dict(file="ds_lattice.py", ref="ds_lattice_output.txt", quick=[], full=[], timeout=H),
                                                dict(file="ds_flips.py", ref="ds_flips_output.txt", quick=[], full=[], timeout=H)]),
    ("13", "papers/Bio_09_Polyploid_Phasing",  "anccft22.tex",           [dict(file="polyploid_phasing.py", ref="polyploid_phasing_output.txt", quick=[], full=[], timeout=H)]),
    ("14", "papers/Bio_10_Pangenome_Sheaf", "anccft23.tex",           [dict(file="pangenome_sheaf.py", ref="pangenome_sheaf_output.txt", quick=[], full=[], timeout=H),
                                                dict(file="real_pangenome.py", ref="real_pangenome_output.txt", quick=None, full=[], timeout=H, out_flag="--out")]),
    ("17", "papers/Bio_11_Isoform_Polytope", "anccft24.tex",           [dict(file="isoform_polytope.py", ref="isoform_polytope_output.txt", quick=[], full=[], timeout=H)]),
    ("8", "papers/Bio_12_Frontier_Elimination", "anccft25.tex",           [dict(file="ds_transfer_matrix.py", ref="ds_transfer_matrix_output.txt", quick=["--quick"], full=[], timeout=4 * H)]),
    ("9", "papers/Bio_13_Fermionic_Partition", "anccft26.tex",           [dict(file="ds_fermion.py", ref="ds_fermion_output.txt", quick=["--quick"], full=[], timeout=H),
                                                dict(file="ds_count.py", ref="ds_count_output.txt", quick=["--quick"], full=[], timeout=3 * H,
                                                     hide=["ds_count_k10_checkpoint.json", "ds_k10_checkpoint.json"])]),
    ("16", "papers/Bio_14_Snarl_Schur", "anccft27.tex",           [dict(file="bio14_verify.py", ref="bio14_output.txt", quick=[], full=[], timeout=H)]),
    ("11", "papers/Bio_17_Inversions_Fibre", "anccft30.tex",           [dict(file="bio17_verify.py", ref="bio17_output.txt", quick=["--quick"], full=[], timeout=H)]),
    ("10", "papers/Bio_19_SharpP_Hardness", "anccft32.tex",           [dict(file="bio19_verify.py", ref="bio19_output.txt", quick=[], full=[], timeout=H)]),
    ("12", "papers/Spectral_Fibres", "spectral_fibres.tex",
                                               [dict(file="spectral_fibres.py", ref="spectral_fibres_output.txt", quick=None, full=[], timeout=2 * H),
                                                dict(file="make_figures.py", ref=None, quick=None, full=[], timeout=H)]),
    ("-", "monograph",      "AG_contents.tex",        []),
    ("-", "monograph",      "AG_titlepage.tex",       []),
    ("F", "code",           None,                     [dict(file="generate_figures.py", ref=None, quick=[], full=[], timeout=H)]),
]

PASS_RE = [re.compile(r"TOTAL\s+(\d+)\s*/\s*(\d+)\s+passed"), re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s+passed", re.M)]
TIME_RE = [re.compile(r"\[[^\]]*\d+(\.\d+)?\s*s[^\]]*\]"), re.compile(r"elapsed\s+[\d.]+\s*s"),
           re.compile(r"eta ~?\d+ ?min"), re.compile(r"\d+(\.\d+)?s\b"), re.compile(r"python \d+\.\d+\.\d+"),
           re.compile(r"\d+\.\d+ GB")]


def read_text(path):
    b = open(path, "rb").read()
    if b[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return b.decode("utf-16", "replace")
    if b[:3] == b"\xef\xbb\xbf":
        b = b[3:]
    t = b.decode("utf-8", "replace")
    if t.count("\x00") > len(t) // 4:           # UTF-16 without BOM
        t = b.decode("utf-16-le", "replace")
    return t


def pass_counts(text):
    for rx in PASS_RE:
        m = rx.findall(text)
        if m:
            a, b = m[-1][:2]
            return int(a), int(b)
    if "ALL CHECKS PASSED" in text:
        return 1, 1
    return None


def normalise(text):
    out = []
    for line in text.splitlines():
        line = line.replace("\t", " ")
        if line.lstrip().startswith(("(base)", "Loading personal", "[Transcript", " the same command", "[intermediate")):
            continue
        for rx in TIME_RE:
            line = rx.sub("", line)
        line = " ".join(line.split())
        if line:
            out.append(line)
    return out


def snapshot(folder):
    snap = {}
    for f in os.listdir(folder):
        p = os.path.join(folder, f)
        if os.path.isfile(p) and f.lower().endswith(".txt"):
            snap[p] = open(p, "rb").read()
    return snap


def restore(snap):
    changed = []
    for p, b in snap.items():
        if not os.path.exists(p) or open(p, "rb").read() != b:
            open(p, "wb").write(b)
            changed.append(os.path.basename(p))
    return changed


def compile_tex(folder, tex, log):
    if shutil.which("pdflatex") is None:
        return "SKIP", "pdflatex not found", 0.0
    t0 = time.time()
    base = tex[:-4]
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode", tex], cwd=folder,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
    lp = os.path.join(folder, base + ".log")
    L = read_text(lp) if os.path.exists(lp) else ""
    errors = len(re.findall(r"^! ", L, re.M))
    undef = len(re.findall(r"(Reference|Citation) `[^']*' .*undefined", L))
    over = len(re.findall(r"^Overfull \\hbox", L, re.M))
    m = re.search(r"Output written on .*?\((\d+) pages?", L)
    pages = m.group(1) if m else "?"
    shutil.copy(lp, os.path.join(log, base + ".log")) if os.path.exists(lp) else None
    status = "PASS" if (errors == 0 and undef == 0 and m) else "FAIL"
    detail = f"{pages} pages, {errors} errors, {undef} undefined refs, {over} overfull"
    return status, detail, time.time() - t0


def run_script(folder, s, mode, log):
    args = s["quick"] if mode == "quick" else s["full"]
    if args is None:
        return "SKIP", "full mode only", 0.0, None
    name = s["file"]
    logfile = os.path.join(log, name[:-3] + f"_{mode}.txt")
    args = list(args)
    if s.get("out_flag"):
        args += [s["out_flag"], logfile]
    hidden = []
    for h in s.get("hide", []):
        p = os.path.join(folder, h)
        if os.path.exists(p):
            os.replace(p, p + ".build_hidden")
            hidden.append(p)
    snap = snapshot(folder)
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", MPLBACKEND="Agg")
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable, "-u", name] + args, cwd=folder, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=s["timeout"] if mode == "full" else 1800)
        out, code = r.stdout.decode("utf-8", "replace"), r.returncode
    except subprocess.TimeoutExpired as e:
        out, code = (e.stdout or b"").decode("utf-8", "replace") + "\n[TIMEOUT]", -9
    dt = time.time() - t0
    restored = restore(snap)
    for p in hidden:
        os.replace(p + ".build_hidden", p)
    if s.get("out_flag") and os.path.exists(logfile):
        out = read_text(logfile) + "\n" + out
    open(logfile, "w", encoding="utf-8").write(out)
    pc = pass_counts(out)
    if code != 0:
        status, detail = "FAIL", f"exit code {code}; last line: {out.strip().splitlines()[-1][:90] if out.strip() else ''}"
    elif s["ref"] is None:
        status, detail = "PASS", "ran (no pass line expected)"
    elif pc is None:
        status, detail = "FAIL", "no pass line found"
    else:
        a, b = pc
        status = "PASS" if a == b else "FAIL"
        detail = f"{a}/{b} passed"
    cmp = None
    if s["ref"]:
        rp = os.path.join(folder, s["ref"])
        if os.path.exists(rp):
            ref = read_text(rp)
            rc_ = pass_counts(ref)
            if rc_:
                detail += f" (reference {rc_[0]}/{rc_[1]})"
            if mode == "full":
                A, B = normalise(ref), normalise(out)
                ndiff = sum(1 for d in difflib.ndiff(A, B) if d[:1] in "+-")
                cmp = "identical after removing timings" if ndiff == 0 else f"{ndiff} lines differ from the reference"
        else:
            detail += " (reference missing)"
    if restored:
        detail += f"; restored {', '.join(restored)}"
    return status, detail, dt, cmp


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--quick", action="store_true")
    g.add_argument("--full", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--no-tex", action="store_true")
    ap.add_argument("--no-scripts", action="store_true")
    ap.add_argument("--inventory", action="store_true")
    a = ap.parse_args()
    mode = "full" if a.full else "quick"
    only = {x.strip() for x in a.only.split(",") if x.strip()}
    os.makedirs(LOGS, exist_ok=True)
    report = []

    def say(s=""):
        print(s, flush=True)
        report.append(s)

    t_all = time.time()
    say(f"=== Algebraic Genomics build ({mode} mode), python {sys.version.split()[0]} ===")
    say(f"    root {ROOT}")
    say("")
    rows = []
    for ch, folder, tex, scripts in CHAPTERS:
        if only and not any(o in folder for o in only):
            continue
        fp = os.path.join(ROOT, folder)
        label = f"ch.{ch:>2} {os.path.basename(folder)[:15]}" if ch not in ("-", "F") else \
            ("figures" if ch == "F" else "monograph")
        if not os.path.isdir(fp):
            rows.append((label, "folder", "FAIL", "missing", 0.0))
            say(f"{label}: folder missing")
            continue
        items = []
        if tex:
            items.append(("tex", tex, os.path.exists(os.path.join(fp, tex))))
        for s in scripts:
            items.append(("py", s["file"], os.path.exists(os.path.join(fp, s["file"]))))
            if s["ref"]:
                items.append(("ref", s["ref"], os.path.exists(os.path.join(fp, s["ref"]))))
        if a.inventory:
            say(f"{label}: " + "; ".join(f"{k}:{n}{'' if ok else ' MISSING'}" for k, n, ok in items))
            continue
        if tex and not a.no_tex:
            if not os.path.exists(os.path.join(fp, tex)):
                rows.append((label, tex, "FAIL", "missing", 0.0))
            else:
                st, de, dt = compile_tex(fp, tex, LOGS)
                rows.append((label, tex, st, de, dt))
                say(f"{label:22s} {tex:28s} {st:4s}  {de}  [{dt:.0f}s]")
        if not a.no_scripts:
            for s in scripts:
                if not os.path.exists(os.path.join(fp, s["file"])):
                    rows.append((label, s["file"], "FAIL", "missing", 0.0))
                    continue
                say(f"{label:22s} {s['file']:28s} running ...")
                st, de, dt, cmp = run_script(fp, s, mode, LOGS)
                rows.append((label, s["file"], st, de + (f"; {cmp}" if cmp else ""), dt))
                say(f"{label:22s} {s['file']:28s} {st:4s}  {de}  [{dt:.0f}s]" + (f"\n{'':52s}{cmp}" if cmp else ""))
    if a.inventory:
        return
    say("")
    say("=== summary ===")
    say(f"{'item':52s} {'status':6s} {'time':>7s}  detail")
    for label, item, st, de, dt in rows:
        say(f"{label + ' / ' + item:52s} {st:6s} {dt:6.0f}s  {de}")
    npass = sum(1 for r in rows if r[2] == "PASS")
    nfail = sum(1 for r in rows if r[2] == "FAIL")
    nskip = sum(1 for r in rows if r[2] == "SKIP")
    say("")
    say(f"TOTAL {npass} passed, {nfail} failed, {nskip} skipped ({mode} mode)   "
        f"elapsed {time.time()-t_all:.0f}s")
    say("ALL ITEMS PASSED" if nfail == 0 else "FAILURES: " + ", ".join(f"{r[0]} / {r[1]}" for r in rows if r[2] == "FAIL"))
    open(os.path.join(ROOT, "build_report.txt"), "w", encoding="utf-8").write("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
