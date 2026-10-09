"""
00_download_faers_quarters.py

Downloads FDA AERS/FAERS quarterly ASCII ZIPs into data/raw/.
Run from the Flox project folder:  python3 scripts/00_download_faers_quarters.py

Source page: https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
URL pattern:
  2004q1-2012q3 : https://fis.fda.gov/content/Exports/aers_ascii_YYYYqN.zip   (legacy AERS, keyed on ISR)
  2012q4-       : https://fis.fda.gov/content/Exports/faers_ascii_YYYYqN.zip  (FAERS, keyed on PRIMARYID)
The quarter letter is lowercase 'q' for most files and uppercase 'Q' for some
(2019-2021 all quarters, 2022Q3-Q4, 2023Q3-Q4, 2024Q4, 2025Q4), per the FDA page.
This script tries the listed case first and the other case as a fallback.

Windows:
  validation : 2007q1 .. 2024q3  (the 2026 Shamim et al. window)
  extension  : 2024q4 .. 2026q2  (latest posted quarter at time of writing)
"""
import sys, time, urllib.request, urllib.error
from pathlib import Path

BASE = "https://fis.fda.gov/content/Exports/"
RAW = Path(__file__).resolve().parent.parent / "data" / "raw"

UPPER = {(2019, q) for q in range(1, 5)} | {(2020, q) for q in range(1, 5)} | \
        {(2021, q) for q in range(1, 5)} | {(2022, 3), (2022, 4), (2023, 3), (2023, 4),
                                            (2024, 4), (2025, 4)}

def quarters(start, end):
    y, q = start
    while (y, q) <= end:
        yield y, q
        q += 1
        if q == 5:
            y, q = y + 1, 1

def names(y, q):
    prefix = "aers" if (y, q) <= (2012, 3) else "faers"
    first = "Q" if (y, q) in UPPER else "q"
    second = "q" if first == "Q" else "Q"
    return [f"{prefix}_ascii_{y}{c}{q}.zip" for c in (first, second)]

def fetch(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)

def main():
    RAW.mkdir(parents=True, exist_ok=True)
    wanted = list(quarters((2007, 1), (2026, 2)))
    failed = []
    for y, q in wanted:
        cands = names(y, q)
        dest = RAW / cands[0].lower()          # normalize local name to lowercase
        if dest.exists() and dest.stat().st_size > 1_000_000:
            print(f"skip   {dest.name}")
            continue
        for n in cands:
            try:
                print(f"get    {n}", flush=True)
                fetch(BASE + n, dest)
                break
            except urllib.error.HTTPError as e:
                if dest.exists(): dest.unlink()
                if e.code != 404:
                    print(f"  HTTP {e.code} on {n}")
            except Exception as e:
                if dest.exists(): dest.unlink()
                print(f"  error on {n}: {e}")
        else:
            failed.append(cands[0])
        time.sleep(1)
    print(f"\nDone. {len(wanted) - len(failed)}/{len(wanted)} files in {RAW}")
    if failed:
        print("Failed (download manually from the FDA page):")
        for n in failed: print("  ", n)
        sys.exit(1)

if __name__ == "__main__":
    main()
