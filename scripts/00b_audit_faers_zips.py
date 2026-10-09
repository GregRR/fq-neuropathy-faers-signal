"""
00b_audit_faers_zips.py

Checks every downloaded FAERS ZIP before parsing, so file-name and layout surprises
show up in one report instead of halfway through script 01.

For each ZIP in data/raw/ it reports:
  * the DEMO, DRUG, REAC and Deleted-list files that script 01 will use
    (flagging zero or several candidates)
  * the header of each table, and any required column that is missing
  * whether the first data row has the same number of '$'-separated fields as the header
    (allowing the one trailing '$' the files carry)
  * the Deleted-list files and their first rows
Then it groups ZIPs by identical header, so a layout change shows up as a new group.

Usage (from the repository root):
    python3 scripts/00b_audit_faers_zips.py            # fast: reads headers only
    python3 scripts/00b_audit_faers_zips.py --crc      # also verifies every ZIP end to end
Writes data/interim/00b_zip_audit.txt as well as printing.
"""
import importlib.util, sys, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("s01", HERE / "01_parse_and_link_core_files.py")
s01 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s01)

REQUIRED = {
    ("L", "demo"): {"isr", "case", "fda_dt", "age", "age_cod", "gndr_cod", "reporter_country"},
    ("F", "demo"): {"primaryid", "caseid", "caseversion", "fda_dt", "init_fda_dt", "age", "age_cod"},
    ("L", "drug"): {"isr", "drug_seq", "role_cod", "drugname"},
    ("F", "drug"): {"primaryid", "drug_seq", "role_cod", "drugname"},
    ("L", "reac"): {"isr", "pt"},
    ("F", "reac"): {"primaryid", "pt"},
}
lines = []


def out(msg=""):
    print(msg, flush=True)
    lines.append(msg)


def first_lines(zf, name, n=2):
    with zf.open(name) as fh:
        raw = [fh.readline() for _ in range(n)]
    return [r.decode("utf-8-sig", errors="replace").rstrip("\r\n") for r in raw]


def main(crc):
    zips = sorted(s01.RAW.glob("*_ascii_*.zip"))
    if not zips:
        sys.exit(f"No ZIPs in {s01.RAW}")
    groups, problems = {}, []
    for zp in zips:
        label = zp.stem.split("_")[-1]
        schema = "L" if zp.name.startswith("aers") else "F"
        out(f"== {zp.name}")
        try:
            zf = zipfile.ZipFile(zp)
            if crc:
                bad = zf.testzip()
                if bad:
                    problems.append(f"{zp.name}: corrupt member {bad}")
        except zipfile.BadZipFile as e:
            problems.append(f"{zp.name}: not a valid ZIP ({e})")
            out(f"  BAD ZIP: {e}")
            continue
        found = {}
        for n in zf.namelist():
            m = s01.MEMBER.match(Path(n).name)
            if m:
                found.setdefault(m.group(1).lower(), []).append(n)
        mem = {}
        for key in ("demo", "drug", "reac"):
            names = found.get(key, [])
            if len(names) != 1:
                problems.append(f"{zp.name}: {len(names)} candidates for {key.upper()}: {names}")
                out(f"  {key.upper()}: {len(names)} candidates {names}")
                continue
            mem[key] = names[0]
            head = first_lines(zf, names[0])
            cols = [c.strip().lower() for c in head[0].split("$")]
            if cols and cols[-1] == "":
                cols = cols[:-1]
            missing = REQUIRED[(schema, key)] - set(cols)
            if key == "demo" and schema == "F" and not ({"sex", "gndr_cod"} & set(cols)):
                missing |= {"sex|gndr_cod"}
            nfields = len(head[1].split("$")) if len(head) > 1 and head[1] else None
            ok_width = nfields in (len(cols), len(cols) + 1)
            out(f"  {key.upper():5s} {Path(names[0]).name}  cols={len(cols)}"
                f"{'  MISSING ' + str(sorted(missing)) if missing else ''}"
                f"{'' if ok_width else f'  FIRST-ROW-FIELDS={nfields}'}")
            if missing:
                problems.append(f"{zp.name}: {key.upper()} missing {sorted(missing)}")
            if not ok_width:
                problems.append(f"{zp.name}: {key.upper()} first row has {nfields} fields, header {len(cols)}")
            groups.setdefault((key, schema, tuple(cols)), []).append(label)
        dels = sorted(n for n in zf.namelist()
                      if len(n.split("/")) > 1 and n.split("/")[-2].lower() == "deleted"
                      and n.lower().endswith(".txt"))
        for d in dels:
            h = first_lines(zf, d, 3)
            out(f"  DELETED {d}  first rows: {h}")
    out("\n==== Header groups (a new group means a layout change) ====")
    for (key, schema, cols), qs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[1][0])):
        out(f"{key.upper()} [{'legacy' if schema == 'L' else 'FAERS'}] {len(qs)} ZIPs: {qs[0]} .. {qs[-1]}")
        out(f"    {'$'.join(cols)}")
    out("\n==== Problems ====")
    for p in problems:
        out("  " + p)
    if not problems:
        out("  none")
    dest = s01.ROOT / "data" / "interim" / "00b_zip_audit.txt"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines))
    print(f"\nWritten to {dest}")


if __name__ == "__main__":
    main("--crc" in sys.argv)
