"""Compare frozen predictions against separate, existing same-assay panels."""
from itertools import combinations, product
import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT.parent
PREDICTIONS = ANALYSIS / "sequence_complementarity_20261108/paired_samples.csv"
MASTOPARAN = ANALYSIS / "cys0_membrane_evidence_20261107/mastoparan_panel_sequence_overlap.csv"
PONERICIN = ANALYSIS / "membrane_expansion_20261025/same_paper_matches.json"
MODELS = ["p_1d", "esm2_mean", "egnn_existing", "esmif_mean"]


def pair_results(data, stronger_weaker):
    results = []
    for strong, weak in stronger_weaker:
        a, b = data.loc[strong], data.loc[weak]
        results.append({"stronger": strong, "weaker": weak, "same_outer_fold": bool(a.fold == b.fold),
                        "ranking": {m: "right" if a[m] > b[m] else "reversed" if a[m] < b[m] else "tie"
                                    for m in MODELS}})
    return results


def interval_pairs(bounds):
    for (a, (alo, ahi)), (b, (blo, bhi)) in combinations(bounds.items(), 2):
        if ahi is not None and ahi < blo:
            yield a, b
        elif bhi is not None and bhi < alo:
            yield b, a


def main():
    pred = pd.read_csv(PREDICTIONS).set_index("sample_id")
    assert len(pred) == 6387 and pred.index.is_unique
    targets = [r["sample_id"] for r in json.loads((OUT / "audit.json").read_text(encoding="utf-8"))["records"]
               if r["all_models_call_nontoxic"] and r["category"] == "damage"]
    assert set(targets) == {f"peptide_{i}" for i in (44, 59, 97, 102, 251, 347, 424, 463)}
    mast = pd.read_csv(MASTOPARAN)
    assert len(mast) == 55
    mast = mast.dropna(subset=["sample_id"]).set_index("sample_id")
    assert len(mast) == 23 and mast.index.is_unique
    pd.testing.assert_series_equal(mast.sequence, pred.loc[mast.index, "sequence"])
    assert (mast.p_1d - pred.loc[mast.index, "p_1d"]).abs().max() < 1e-10
    pct = mast.human_rbc_hemolysis_320_ug_ml.str.rstrip("%").astype(float)
    ec = pd.to_numeric(mast.human_rbc_ec50_uM.str.extract(r"^\s*(\d+(?:\.\d+)?)\s*±")[0], errors="coerce").dropna()
    assert len(ec) == 12
    mast_ranks = {
        "matched_n": len(mast), "group_n": int(pred.loc[mast.index].group_id.nunique()),
        "n_with_ec50": len(ec),
        "spearman_percent": {m: float(spearmanr(pred.loc[pct.index, m], pct).statistic) for m in MODELS},
        "spearman_negative_ec50": {m: float(spearmanr(pred.loc[ec.index, m], -ec).statistic) for m in MODELS},
        "cases": {sid: {"name": mast.loc[sid, "peptide"], "form": mast.loc[sid, "form_in_table"],
                        "human_hemolysis_percent": float(pct.loc[sid]),
                        "human_ec50_uM": float(ec.loc[sid]) if sid in ec.index else None}
                  for sid in mast.index},
    }

    nc = {row["sample"]["sample_id"]: row["sequence"] for row in json.loads(PONERICIN.read_text(encoding="utf-8-sig"))
          if row["present_in_training"]}
    assert set(nc) == {f"peptide_{i}" for i in (249, 3912, 463, 444)}
    assert all(pred.loc[sid, "sequence"] == seq for sid, seq in nc.items())
    cc50 = {"peptide_249": (4.5, 5.9), "peptide_3912": (12.6, 13.5),
            "peptide_463": (48.2, 56.9), "peptide_444": (300, None)}
    hc50 = {"peptide_249": (39.9, 61.5), "peptide_3912": (28.6, 48.2),
            "peptide_463": (250, None), "peptide_444": (300, None)}
    # Only rank pairs whose reported intervals or censoring bounds cannot overlap.
    assert len(list(interval_pairs(cc50))) == 6 and len(list(interval_pairs(hc50))) == 4
    fish_detected = ["peptide_58", "peptide_59"]
    fish_not_detected = ["peptide_1864", "peptide_3463", "peptide_5786"]
    profile = pd.read_csv(ANALYSIS / "six_seed_1d_weakness_20261107/sample_profiles.csv").set_index("sample_id")
    for sid, name in {"peptide_58": "Gs 1", "peptide_59": "Gs 2", "peptide_1864": "Gs A",
                      "peptide_3463": "Gs C", "peptide_5786": "Gs B", "peptide_347": "Pp 3"}.items():
        assert name in profile.loc[sid, "name"] and profile.loc[sid, "sequence"] == pred.loc[sid, "sequence"]
    result = {
        "mastoparan": mast_ranks,
        "ponericin": {
            "HEK293_CC50": pair_results(pred, interval_pairs(cc50)),
            "human_RBC_HC50": pair_results(pred, interval_pairs(hc50)),
            "CC50_reported_ranges_uM": cc50, "HC50_reported_ranges_uM": hc50,
        },
        "grammistin": {
            "detected_vs_not_detected_2005": pair_results(pred, product(fish_detected, fish_not_detected)),
            "Gs2_vs_Gs1_2000": pair_results(pred, [("peptide_59", "peptide_58")]),
            "Pp3_terminal_pair_2006": {"sample_id": "peptide_347",
                                      "sequence_only_same_score_by_construction": True},
        },
        "H1_H2": {"ids": ["peptide_251", "peptide_424"], "quantitative_pair_included": False,
                  "reason": "Primary quantitative table and terminal-state conflict unresolved."},
        "eight_target_ids": sorted(targets),
        "test1_test2_scores_used": False,
    }
    scored = set(mast.index) | set(nc) | set(fish_detected + fish_not_detected) | set(targets)
    result["frozen_scores"] = pred.loc[sorted(scored), ["fold", "sequence", "label", "group_id", *MODELS]].to_dict("index")
    (OUT / "same_assay_comparisons.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    for section, key in [("ponericin", "HEK293_CC50"), ("ponericin", "human_RBC_HC50"),
                         ("grammistin", "detected_vs_not_detected_2005"), ("grammistin", "Gs2_vs_Gs1_2000")]:
        pairs = result[section][key]
        print(key, len(pairs), {m: sum(row["ranking"][m] == "right" for row in pairs) for m in MODELS})
    print({key: value for key, value in mast_ranks.items() if key != "cases"})
    for sid in sorted(scored - set(mast.index)):
        print(sid, pred.loc[sid, ["fold", *MODELS]].to_dict())


if __name__ == "__main__":
    main()
