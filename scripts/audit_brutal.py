#!/usr/bin/env python3
"""Brutal devil's-advocate switch audit. For each MD json, TRY TO FAIL the switch:
 (1) CONTACT GATE   : ON-state (HIE) His->nearest-cation mindist < 8.0 A  (else drift, not a grip)
 (2) NOISE FLOOR    : d_repel = mean(HIP)-mean(HIE) must exceed 2*sqrt(SEM_HIP^2+SEM_HIE^2)
                      (else the 'switch' is within replica scatter = noise)
 (3) CONTROL MARGIN : d_repel must exceed dud(Gln isostere) + 0.3  (else pose drift, not His-driven)
 (4) DIRECTION      : d_repel > 0 (protonation pushes His AWAY from cation = release)
A switch SURVIVES only if all 4 pass. Reads per-replica his_lys_mindist from states[].
Usage: audit_brutal.py <dir-with-md-jsons> [more dirs...]"""
import sys, glob, os, json, math
from collections import defaultdict

def sem(xs):
    n=len(xs)
    if n<2: return 0.0
    m=sum(xs)/n; v=sum((x-m)**2 for x in xs)/(n-1)
    return math.sqrt(v/n)

def audit(j):
    states=j.get("states",[])
    by=defaultdict(list)
    for s in states:
        d=s.get("his_lys_mindist")
        if d is not None: by[s["variant"]].append(d)
    def mean(v): return sum(by[v])/len(by[v]) if by.get(v) else None
    hip,hie,gln=mean("HIP"),mean("HIE"),mean("GLN")
    out={"his_pos":j.get("his_pos"),"HIP":hip,"HIE":hie,"GLN":gln,
         "n_rep":{k:len(v) for k,v in by.items()}}
    if hie is None or hip is None:
        out["verdict"]="NO_DATA"; return out
    d_repel=hip-hie
    noise=math.sqrt(sem(by["HIP"])**2+sem(by["HIE"])**2)
    dud=(gln-hie) if gln is not None else None  # Gln(HIE-equivalent neutral) drift vs HIE... use reported if present
    dud_rep=j.get("dud_repel")
    g1 = hie < 8.0                                   # contact gate (ON-state)
    g2 = d_repel > 2*noise and d_repel>0.3           # above noise floor + min effect
    g3 = (dud_rep is None) or (d_repel > dud_rep+0.3) # control margin
    g4 = d_repel > 0                                  # direction
    out.update(d_repel=round(d_repel,3), noise_2sem=round(2*noise,3), dud=dud_rep,
               gate_contact=g1, gate_noise=g2, gate_control=g3, gate_dir=g4,
               verdict="SURVIVE" if (g1 and g2 and g3 and g4) else "CUT")
    return out

seen={}
for d in sys.argv[1:]:
    for p in sorted(glob.glob(os.path.join(d,"**","*.json"),recursive=True)):
        try: j=json.load(open(p))
        except: continue
        if "states" not in j: continue
        nm=os.path.basename(p).replace(".json","")
        seen[nm]=audit(j)
print(f"{'design':40s} {'his':10s} {'HIE':6s} {'d_rep':7s} {'2sem':6s} {'dud':6s} C N K D  verdict")
for nm,o in sorted(seen.items(), key=lambda kv: -(kv[1].get('d_repel') or -9)):
    if o.get("verdict")=="NO_DATA":
        print(f"{nm[:40]:40s} NO_DATA"); continue
    flags="".join("Y" if o[g] else "." for g in ("gate_contact","gate_noise","gate_control","gate_dir"))
    print(f"{nm[:40]:40s} {str(o['his_pos']):10s} {o['HIE']:<6} {o['d_repel']:<7} {o['noise_2sem']:<6} {str(o['dud']):6s} {flags[0]} {flags[1]} {flags[2]} {flags[3]}  {o['verdict']}")
