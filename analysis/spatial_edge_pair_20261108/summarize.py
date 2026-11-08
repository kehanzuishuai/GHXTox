"""Score the frozen test benchmarks and the matched-edge fold-0 pilot."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT.parent
REVIEW = ANALYSIS / "reviewer_response_20261102"
IDS = ["sample_id", "sequence", "label"]


def metrics(y, p, pred):
    return {"mcc": matthews_corrcoef(y, pred), "auroc": roc_auc_score(y, p),
            "auprc": average_precision_score(y, p), "brier": brier_score_loss(y, p),
            "fn": int(((y == 1) & (pred == 0)).sum()), "fp": int(((y == 0) & (pred == 1)).sum())}


def differences(y, baseline, candidate):
    return {"corrected": int(((baseline != y) & (candidate == y)).sum()),
            "added_errors": int(((baseline == y) & (candidate != y)).sum())}


def summarize_test_sets():
    results = {}
    for split, n in [("test1", 1126), ("test2", 582)]:
        old = pd.read_csv(ANALYSIS / f"test_comparison_current_20261106/{split}_paired.csv")
        six = pd.read_csv(REVIEW / f"three_seed_full_received_20261106/{split}_predictions.csv")
        joined = old.merge(six, on=IDS, validate="one_to_one")
        assert len(joined) == len(old) == len(six) == n and joined.sample_id.is_unique
        y = joined.label.to_numpy()
        sequence = joined.six_sequence.to_numpy()
        fusion = joined.six_plus_3d.to_numpy()
        threshold = 0.47022256531284556
        assert np.allclose(fusion, .72 * sequence + .28 * joined.p3d)
        baseline = sequence >= threshold
        models = {"six_1d": (sequence, baseline), "six_plus_3d": (fusion, fusion >= threshold)}
        for name in ["GHXTox", "ToxMSRC", "ToxPLTC", "StrucToxNet", "PeptiTox"]:
            models[name] = (joined[f"{name}_probability"].to_numpy(), joined[f"{name}_prediction"].to_numpy())
        results[split] = {name: metrics(y, p, pred) | differences(y, baseline, pred)
                          for name, (p, pred) in models.items()}
    (OUT / "historical_test_comparison.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


def summarize_fold():
    run = OUT / "fold0"
    assert json.loads((run / "metadata.json").read_text())["status"] == "complete"
    package = REVIEW / "six_seed_nested_autodl_20261106"
    spec = importlib.util.spec_from_file_location("calibrator", package / "compare.py")
    calibrator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(calibrator)
    manifest = pd.read_csv(package / "data/nested_folds.csv")
    held = manifest[(manifest.outer_fold == 0) & manifest.role.eq("test")][IDS]
    cal = manifest[(manifest.outer_fold == 0) & manifest.role.eq("calibration")][IDS]
    profile = pd.read_csv(ANALYSIS / "six_seed_1d_weakness_20261107/sample_profiles.csv")
    data = held.merge(profile[profile.fold.eq(0)], on=IDS, validate="one_to_one")
    assert len(data) == len(held) == 1271
    for mode in ("matched_sequence", "matched_spatial"):
        cal_score, held_score = [pd.read_csv(run / mode / f"{role}_predictions.csv") for role in ("calibration", "test")]
        assert len(cal_score) == len(cal) and len(held_score) == len(held)
        calibrated = cal.merge(cal_score[IDS + ["toxicity_probability"]], on=IDS, validate="one_to_one")
        assert len(calibrated) == len(cal)
        parameters = calibrator.fit(calibrated.toxicity_probability.to_numpy(), calibrated.label.to_numpy())
        data = data.merge(held_score[IDS + ["toxicity_probability"]].rename(columns={"toxicity_probability": mode}),
                          on=IDS, validate="one_to_one")
        assert len(data) == len(held)
        data[mode] = calibrator.apply(data[mode].to_numpy(), parameters)
        data[f"six_plus_{mode}"] = .72 * data.p_1d + .28 * data[mode]
    y = data.label.to_numpy()
    base = data.six_plus_matched_sequence.to_numpy() >= .5
    model_names = ["p_1d", "matched_sequence", "matched_spatial", "six_plus_matched_sequence", "six_plus_matched_spatial"]
    result = {name: metrics(y, data[name].to_numpy(), data[name].to_numpy() >= .5) |
                    differences(y, base, data[name].to_numpy() >= .5) for name in model_names}
    spatial = data.six_plus_matched_spatial.to_numpy() >= .5
    groups = [rows.index.to_numpy() for _, rows in data.groupby("group_id")]
    rng = np.random.default_rng(20261108)
    deltas = []
    for _ in range(2000):
        sample = np.concatenate([groups[i] for i in rng.integers(len(groups), size=len(groups))])
        deltas.append(matthews_corrcoef(y[sample], spatial[sample]) - matthews_corrcoef(y[sample], base[sample]))
    result["matched_fusion_comparison"] = {
        "delta_mcc": matthews_corrcoef(y, spatial) - matthews_corrcoef(y, base),
        "group_bootstrap_95ci": np.quantile(deltas, [.025, .975]).tolist(),
        "groups": len(groups), "draws": len(deltas)}
    cohorts = {"weaker_reverse": data.band.eq("weaker_reverse"),
               "cys0_toxic_strong_reverse": data.cys.eq(0) & data.label.eq(1) & data.band.eq("strong_reverse"),
               "cys0_nontoxic": data.cys.eq(0) & data.label.eq(0)}
    for name, mask in cohorts.items():
        cohort = data[mask]
        result[name] = {"n": len(cohort)}
        for mode in model_names:
            p = np.where(cohort.label.eq(1), cohort[mode], 1 - cohort[mode])
            result[name][mode] = {"correct": int(((cohort[mode] >= .5) == cohort.label).sum()),
                                  "bce": float(-np.log(np.maximum(1e-12, p)).mean())}
    data[IDS + ["group_id", "p_1d", "band", "cys", *model_names[1:]]].to_csv(run / "paired_scores.csv", index=False)
    (run / "summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scope", choices=["test", "fold"])
    (summarize_test_sets if parser.parse_args().scope == "test" else summarize_fold)()
