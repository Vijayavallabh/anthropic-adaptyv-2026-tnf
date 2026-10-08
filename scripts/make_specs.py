#!/usr/bin/env python3
"""Emit a diverse set of BoltzGen binder specs for the EGFR Domain III campaign.

Three hotspot-focus variants x two length bands = 6 specs, all sharing the
cross-reactivity steer (not_binding = human/mouse-divergent cetuximab contacts).
Ordinals are read from data/targets/6ARU_domainIII_map.json so this stays in sync
with the crop.

  H  His433 sensor core            (P00533 408,432,433,436)
  A  His433 + Glu496 acidic pair   (P00533 432,433,436,489,496)
  C  full conserved patch          (P00533 408,432,433,436,464,489,496)
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MAP = json.loads((REPO / "data" / "targets" / "6ARU_domainIII_map.json").read_text())
PDB = str(REPO / "data" / "targets" / "6ARU_domainIII.pdb")
SPECS = REPO / "specs"
SPECS.mkdir(exist_ok=True)

cons = {v["p00533"]: v["ordinal"] for v in MAP["conserved_binding"].values()}
div = {v["p00533"]: v["ordinal"] for v in MAP["divergent_not_binding"].values()}
NOTBIND = ",".join(str(div[p]) for p in sorted(div))

VARIANTS = {
    "H": [408, 432, 433, 436],
    "A": [432, 433, 436, 489, 496],
    "C": [408, 432, 433, 436, 464, 489, 496],
}
DESCR = {
    "H": "His433 sensor core",
    "A": "His433 + Glu496 acidic pair",
    "C": "full conserved cross-reactive patch",
}
BANDS = {"s": "55..80", "l": "80..110"}


def main() -> None:
    written = []
    for vk, ps in VARIANTS.items():
        missing = [p for p in ps if p not in cons]
        if missing:
            raise SystemExit(f"variant {vk}: P00533 {missing} not in conserved map")
        bind = ",".join(str(cons[p]) for p in ps)
        for bk, rng in BANDS.items():
            name = f"egfr_dom3_{vk}_{bk}"
            text = (
                f"# EGFR Domain III de novo binder — {DESCR[vk]}; binder length {rng}.\n"
                f"# binding = conserved anchors P00533 {ps} (His433=ord100, Glu496=ord163);\n"
                f"# not_binding = divergent cetuximab contacts (human/mouse) for cross-reactivity.\n"
                "entities:\n"
                "  - protein:\n"
                "      id: B\n"
                f"      sequence: {rng}\n"
                "  - file:\n"
                f"      path: {PDB}\n"
                "      include:\n"
                "        - chain:\n"
                "            id: A\n"
                "      binding_types:\n"
                "        - chain:\n"
                "            id: A\n"
                f'            binding: "{bind}"\n'
                f'            not_binding: "{NOTBIND}"\n'
                '      structure_groups: "all"\n'
            )
            (SPECS / f"{name}.yaml").write_text(text)
            written.append((name, bind))
            print(f"wrote {name}: binding={bind}")
    print(f"not_binding (all): {NOTBIND}")
    print(f"{len(written)} specs -> {SPECS}")


if __name__ == "__main__":
    main()
