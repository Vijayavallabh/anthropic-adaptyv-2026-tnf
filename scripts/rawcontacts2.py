import glob,sys,numpy as np
from collections import defaultdict
HOT={('B',114),('B',115),('B',141),('B',65)}
def load(fn):
    at=defaultdict(list); hot=[]
    for ln in open(fn):
        if ln[:4]!="ATOM": continue
        ch=ln[21]; rn=int(ln[22:26]); xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
        at[ch].append(xyz)
        if (ch,rn) in HOT: hot.append(xyz)
    return at,(np.array(hot) if hot else np.zeros((0,3)))
pat=sys.argv[1]
rows=[]
for fn in sorted(glob.glob(pat)):
    at,hot=load(fn)
    if 'C' not in at: continue
    C=np.array(at['C']); T=np.concatenate([np.array(at[c]) for c in at if c!='C'])
    d=np.sqrt(((C[:,None]-T[None])**2).sum(-1)); ncont=int((d<5.0).sum())
    heng=0; hmin=np.nan
    if len(hot):
        dh=np.sqrt(((C[:,None]-hot[None])**2).sum(-1)); hmin=float(dh.min())
        heng=int((dh.min(0)<6.0).sum())  # hotspot atoms within 6A of binder
    rows.append((ncont,hmin,heng,fn))
if not rows: print("no designs"); sys.exit()
import statistics as st
ncs=[r[0] for r in rows]
print(f"n={len(rows)}  contacts<5A: min={min(ncs)} med={int(st.median(ncs))} max={max(ncs)} mean={st.mean(ncs):.0f}  | >20contacts:{sum(1 for n in ncs if n>20)}/{len(rows)}  >40:{sum(1 for n in ncs if n>40)}")
print(f"  hotspot min-dist med={st.median([r[1] for r in rows if not np.isnan(r[1])]):.1f}A  top3 by contacts: {[r[0] for r in sorted(rows,reverse=True)[:3]]}")
