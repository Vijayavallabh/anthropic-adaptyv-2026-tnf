#!/usr/bin/env python3
"""Fetch an unpaired MSA (a3m) from the ColabFold MMseqs2 API for ONE sequence.

AlphaFold3 degrades badly single-sequence and the local genetic DBs are absent, so we source the
target's MSA from the public ColabFold server (https-only, de-novo binders stay single-sequence).
Usage: fetch_colabfold_msa.py <seq|fastafile> <name> <out.a3m> [--mode env] [--timeout 1800]
"""
import sys, time, json, tarfile, io, os, urllib.request, urllib.parse

HOST = "https://api.colabfold.com"
UA = {"User-Agent": "adaptyv-tnf-af3-research (noncommercial academic)"}

def _get(url, data=None):
    req = urllib.request.Request(url, data=data, headers=UA)
    return urllib.request.urlopen(req, timeout=300)

def submit(seq, mode):
    data = urllib.parse.urlencode({"q": f">query\n{seq}\n", "mode": mode}).encode()
    r = json.load(_get(HOST + "/ticket/msa", data))
    tid = r.get("id")
    if not tid:   # RATELIMIT / MAINTENANCE / error responses omit "id" -> fail with the real cause, not a KeyError
        raise RuntimeError(f"ColabFold MSA submit returned no ticket id (status={r.get('status', r)})")
    return tid

def poll(tid, timeout):
    t0 = time.time()
    while time.time() - t0 < timeout:
        st = json.load(_get(HOST + f"/ticket/{tid}")).get("status")
        if st == "COMPLETE":
            return
        if st in ("ERROR", "MAINTENANCE", "UNKNOWN"):
            raise RuntimeError(f"ColabFold ticket {st}")
        time.sleep(10)
    raise TimeoutError("MSA poll timeout")

def download(tid):
    blob = _get(HOST + f"/result/download/{tid}").read()
    tf = tarfile.open(fileobj=io.BytesIO(blob))
    return {m.name: tf.extractfile(m).read().decode() for m in tf.getmembers() if m.name.endswith(".a3m")}

def main():
    arg, name, out = sys.argv[1], sys.argv[2], sys.argv[3]
    def _flagval(flag, default):   # safe: error cleanly if flag is the last token with no value
        if flag not in sys.argv: return default
        i = sys.argv.index(flag)
        if i + 1 >= len(sys.argv): sys.exit(f"fetch_colabfold_msa: {flag} requires a value")
        return sys.argv[i + 1]
    mode = _flagval("--mode", "env")
    try:
        timeout = float(_flagval("--timeout", 1800))
    except ValueError:
        sys.exit("fetch_colabfold_msa: --timeout must be a number")
    raw = open(arg).read() if os.path.exists(arg) else arg
    seq = "".join(l.strip() for l in raw.splitlines() if l.strip() and not l.startswith(">"))
    print(f"[msa] {name}: submit len={len(seq)} mode={mode}", flush=True)
    tid = submit(seq, mode); print(f"[msa] {name}: ticket {tid}", flush=True)
    poll(tid, timeout); print(f"[msa] {name}: COMPLETE", flush=True)
    a3ms = download(tid)
    if not a3ms:
        raise RuntimeError("no a3m in download")
    key = next((k for k in a3ms if "uniref" in k.lower()), None) or max(a3ms, key=lambda k: len(a3ms[k]))
    a3m = a3ms[key]
    # ensure the query (exact input seq, ungapped uppercase) is the first record (AF3 requirement)
    recs = a3m.replace("\r", "").split("\n>")
    if recs and seq not in recs[0]:
        a3m = f">query\n{seq}\n" + a3m
    open(out, "w").write(a3m)
    print(f"[msa] {name}: wrote {out} from {key} — {a3m.count('>')} seqs, {len(a3m)} chars", flush=True)

if __name__ == "__main__":
    main()
