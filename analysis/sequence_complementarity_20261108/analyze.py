"""Compare cached sequence and structure representations on six-seed 1D errors."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
ASSETS = ROOT.parent / "GHXTox-ziliao"
REVIEW = ROOT / "analysis/reviewer_response_20261102"
PACKAGE = REVIEW / "six_seed_nested_autodl_20261106"
OLD = REVIEW / "nested_1d_received_20261106/nested_1d_results/runs"
NEW = REVIEW / "six_seed_received_20261106/six_seed_new_results/runs"
SHOOTOUT = ROOT / "runs/representation_shootout_v1"
STRUCTURE = OUT.parent / "structure_expert_20261107"
PROFILE = OUT.parent / "six_seed_1d_weakness_20261107/sample_profiles.csv"
SPLIT = PACKAGE / "data/nested_folds.csv"
BRIDGE = ROOT / "data/folds/prottox_tsk_train_prott5_official_mean_folds_v1.csv"
ESM = ASSETS / "train_chemical_sites_final_esm2_multisource_1_30_full2795.pt"
T5 = ASSETS / "train_prott5_official_mean.pt"
IF1 = STRUCTURE / "esm_if_train_mean.npz"
KEYS = ["sample_id", "sequence", "label", "fold", "group_id"]
SEEDS = (42, 123, 2025, 314, 2718, 2026)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tensor_cache(path):
    return torch.load(path, map_location="cpu", weights_only=False, mmap=True)


def cache_index(ids, sequences, data, labels=None):
    assert len(ids) == len(set(ids)) == len(data)
    rows = {sid: i for i, sid in enumerate(ids)}
    indices = np.array([rows[sid] for sid in data.sample_id])
    np.testing.assert_array_equal(np.asarray(sequences)[indices], data.sequence)
    if labels is not None:
        np.testing.assert_array_equal(np.asarray(labels).reshape(-1)[indices], data.label)
    return indices


def sequence_features(data, bridge):
    esm = tensor_cache(ESM)
    assert esm["plm_model"] == "esm2_t33_650M_UR50D"
    rows = esm["records"]
    index = cache_index([r["sample_id"] for r in rows], [r["sequence"] for r in rows],
                        data, [r["label"] for r in rows])
    pooled = []
    for i in index:
        row = rows[i]
        tokens = row["plm_features"]
        assert tokens.shape == (len(row["sequence"]), 1280) and row["plm_model"] == esm["plm_model"]
        pooled.append(tokens.float().mean(0).numpy())
    features = {"esm2_mean": np.stack(pooled)}
    del rows, esm

    t5 = tensor_cache(T5)
    official = data.assign(sample_id=data.sample_id.map(bridge.set_index("historical_sample_id").sample_id))
    index = cache_index(t5["sample_ids"], t5["sequences"], official, t5["labels"])
    features["prott5_mean"] = np.asarray(t5["features"])[index]
    with np.load(IF1, allow_pickle=False) as z:
        index = cache_index(z["sample_ids"], z["sequences"], data)
        features["esmif_mean"] = z["features"][index]
    history, source_paths = {}, []
    bridge_index = bridge.set_index("sample_id")
    for size in ("large", "small"):
        path = SHOOTOUT / f"pepbert_{size}_uniparc_features.pt"
        saved = tensor_cache(path)["records"]
        index = cache_index([r["sample_id"] for r in saved], [r["sequence"] for r in saved],
                            official, [r["label"] for r in saved])
        tri = np.stack([np.asarray(saved[i]["tristat"]) for i in index])
        features[f"pepbert_{size}_mean"] = tri[:, :tri.shape[1] // 3]
        scores = []
        for seed in (20260905, 20260906):
            file = SHOOTOUT / f"pepbert_{size}_uniparc_seed{seed}_v1/oof_predictions.csv"
            source_paths.append(file)
            run = pd.read_csv(file).set_index("sample_id").loc[official.sample_id]
            assert len(run) == len(data) and run.index.is_unique
            np.testing.assert_array_equal(run.label, data.label)
            np.testing.assert_array_equal(run.fold, data.fold)
            np.testing.assert_array_equal(run.sequence_sha256, data.sequence.map(
                lambda s: hashlib.sha256(s.encode()).hexdigest()))
            np.testing.assert_array_equal(bridge_index.loc[run.index].fold, run.fold)
            scores.append(run.probability.to_numpy())
        history[f"pepbert_{size}_historical"] = np.mean(scores, axis=0)
        source_paths.append(path)
    features["prott5_esm2_mean"] = np.concatenate([features["prott5_mean"], features["esm2_mean"]], axis=1)
    features["esm2_esmif_mean"] = np.concatenate([features["esm2_mean"], features["esmif_mean"]], axis=1)
    assert all(len(x) == len(data) and np.isfinite(x).all() for x in features.values())
    return features, history, source_paths


def equal_fp_threshold(probability, allowance):
    ranked = np.sort(probability)[::-1]
    assert 0 <= allowance < len(ranked)
    threshold = float(np.nextafter(ranked[allowance], np.inf))
    assert int((probability >= threshold).sum()) <= allowance
    return threshold


def fit_probes(data, manifest, features, prior):
    meta = json.loads((PROFILE.parent / "provenance.json").read_text())
    base_calibration = {entry["fold"]: entry for entry in meta["calibration"]}
    data["egnn_existing_matched_decision"] = np.nan
    for col in features:
        data[col] = np.nan
        data[f"{col}_matched_decision"] = np.nan
    audits = []
    for fold, roles in manifest.groupby("outer_fold"):
        roles = {name: rows.set_index("sample_id") for name, rows in roles.groupby("role")}
        groups = [set(rows.group_id) for rows in roles.values()]
        assert sum(map(len, groups)) == len(set.union(*groups))
        indices = {name: data.index.get_indexer(rows.index) for name, rows in roles.items()}
        for name, rows in roles.items():
            assert (indices[name] >= 0).all() and rows.index.is_unique
            for key in ("sequence", "label", "group_id"):
                np.testing.assert_array_equal(data.loc[rows.index, key], rows[key])
        train, cal, test = [indices[k] for k in ("train", "calibration", "test")]
        assert (data.iloc[test].fold == fold).all() and len(test) == (data.fold == fold).sum()
        y_train = data.iloc[train].label.to_numpy()
        y_cal = data.iloc[cal].label.to_numpy()
        raw = []
        for seed in SEEDS:
            run = OLD if seed in SEEDS[:3] else NEW
            file = run / f"seed{seed}_fold{fold}/calibration_predictions.csv"
            assert sha(file) == meta["seed_files_sha256"][str(file.relative_to(ROOT)).replace("/", "\\")]
            pred = pd.read_csv(file).set_index("sample_id").loc[roles["calibration"].index]
            np.testing.assert_array_equal(pred.label, y_cal)
            raw.append(pred.probability.to_numpy())
        baseline_raw = np.mean(raw, axis=0)
        saved_cal = base_calibration[fold]
        assert (saved_cal["calibration_n"], saved_cal["held_n"]) == (len(cal), len(test))
        base_params = (saved_cal["slope"], saved_cal["intercept"])
        base_cal = prior.apply(baseline_raw, base_params)
        np.testing.assert_allclose(
            prior.apply(data.iloc[test].p_1d_raw_mean.to_numpy(), base_params),
            data.iloc[test].p_1d, atol=1e-8, rtol=0)
        allowance = int(((base_cal >= .5) & (y_cal == 0)).sum())
        structure_cal = REVIEW / f"six_seed_received_20261106/structure/fold{fold}/calibration.csv"
        assert sha(structure_cal) == prior.STRUCT_SHA[(fold, "calibration")]
        structure_test = structure_cal.with_name("test.csv")
        assert sha(structure_test) == prior.STRUCT_SHA[(fold, "test")]
        raw_3d = pd.read_csv(structure_cal).set_index("sample_id").loc[roles["calibration"].index]
        raw_3d_test = pd.read_csv(structure_test).set_index("sample_id").loc[roles["test"].index]
        np.testing.assert_array_equal(raw_3d.label, y_cal)
        np.testing.assert_array_equal(raw_3d.sequence, data.iloc[cal].sequence)
        np.testing.assert_array_equal(raw_3d_test.sequence, data.iloc[test].sequence)
        np.testing.assert_array_equal(raw_3d_test.label, data.iloc[test].label)
        neg_3d = raw_3d.toxicity_probability.to_numpy()[y_cal == 0]
        threshold_3d = equal_fp_threshold(neg_3d, allowance)
        data.loc[data.index[test], "egnn_existing_matched_decision"] = (
            raw_3d_test.toxicity_probability.to_numpy() >= threshold_3d).astype(int)
        audits.append(dict(model="egnn_existing", fold=int(fold), train=None, cal=len(cal), test=len(test),
                           base_cal_fp=allowance, candidate_cal_fp=int((neg_3d >= threshold_3d).sum()),
                           matched_raw_threshold=threshold_3d))
        for name, x in features.items():
            model = make_pipeline(StandardScaler(), LogisticRegression(C=.1, max_iter=1000))
            model.fit(x[train], y_train)
            assert model[-1].n_iter_[0] < 1000
            cal_raw = model.predict_proba(x[cal])[:, 1]
            test_raw = model.predict_proba(x[test])[:, 1]
            params = prior.fit(cal_raw, y_cal)
            probability = prior.apply(test_raw, params)
            threshold = equal_fp_threshold(cal_raw[y_cal == 0], allowance)
            data.loc[data.index[test], name] = probability
            data.loc[data.index[test], f"{name}_matched_decision"] = (test_raw >= threshold).astype(int)
            audits.append(dict(model=name, fold=int(fold), train=len(train), cal=len(cal), test=len(test),
                               base_cal_fp=allowance, candidate_cal_fp=int((cal_raw[y_cal == 0] >= threshold).sum()),
                               matched_raw_threshold=threshold, platt_slope=params[0], platt_intercept=params[1]))
    return audits


def matched_controls(data, cases, label, low_score):
    pool = data[data.label.eq(label) & (data.p_1d.le(.1) if low_score else data.p_1d.ge(.9))]
    matched = set()
    for case in cases.itertuples():
        eligible = pool[pool.fold.eq(case.fold) & pool.group_id.ne(case.group_id)
                        & pool.cys.eq(case.cys) & pool.length.sub(case.length).abs().le(3)]
        matched.update(eligible.sample_id)
    return data.index.isin(matched)


def summarize(data, models):
    label = data.label
    cys0 = data.cys.eq(0)
    cys7 = data.cys.ge(7)
    strong = data.band.eq("strong_reverse")
    low = label.eq(1) & cys0 & strong
    high = label.eq(0) & cys7 & strong
    masks = {
        "all": np.ones(len(data), dtype=bool), "cys0_toxic_strong": low,
        "cys7_nontoxic_strong": high, "all_cys0_toxic": label.eq(1) & cys0,
        "all_cys7_nontoxic": label.eq(0) & cys7, "all_cys0_nontoxic": label.eq(0) & cys0,
        "all_cys7_toxic": label.eq(1) & cys7,
        "weaker_reverse": data.band.eq("weaker_reverse"), "near_neutral": data.band.eq("near_neutral"),
        "high_seed_disagreement": data.high_seed_disagreement,
        "matched_cys0_nontoxic": matched_controls(data, data[low], 0, True),
        "matched_cys7_toxic": matched_controls(data, data[high], 1, False),
    }
    assert [int(masks[k].sum()) for k in list(masks)[:3]] == [6387, 67, 6]
    assert masks["matched_cys0_nontoxic"].sum() == 1889
    for cohort, mask in masks.items():
        data[f"cohort_{cohort}"] = np.asarray(mask)
    baseline = data.p_1d.ge(.5).eq(label)
    summary = []
    family_rows = []
    for model, probability, decision, mode in models:
        assert data[probability].between(0, 1).all()
        d = data[probability].ge(.5) if decision is None else data[decision]
        assert d.isin([0, 1]).all()
        correct = d.eq(label)
        outcome = np.select([~baseline & correct, baseline & ~correct, baseline],
                            ["rescued", "added_error", "both_correct"], default="both_wrong")
        data[f"{model}_{mode}_outcome"] = outcome
        for cohort, mask in masks.items():
            subset = data[np.asarray(mask)]
            for fold in ("all", 0, 1, 2, 3, 4):
                part = subset if fold == "all" else subset[subset.fold.eq(fold)]
                idx = part.index
                y, pred = label.loc[idx], d.loc[idx]
                rescued = int((~baseline.loc[idx] & correct.loc[idx]).sum())
                added = int((baseline.loc[idx] & ~correct.loc[idx]).sum())
                err = int((pred != y).sum())
                assert err == int((~baseline.loc[idx]).sum()) - rescued + added
                p_true = np.where(y.eq(1), part[probability], 1 - part[probability])
                summary.append(dict(model=model, mode=mode, cohort=cohort, fold=fold, n=len(idx),
                                    groups=part.group_id.nunique(), baseline_errors=int((~baseline.loc[idx]).sum()),
                                    errors=err, rescued=rescued, added_error=added,
                                    false_negative=int(((pred == 0) & y.eq(1)).sum()),
                                    false_positive=int(((pred == 1) & y.eq(0)).sum()),
                                    mean_bce=float(-np.log(np.clip(p_true, 1e-12, 1)).mean()) if len(idx) else None,
                                    mcc=matthews_corrcoef(y, pred) if y.nunique() == 2 else None,
                                    auprc=average_precision_score(y, part[probability]) if y.nunique() == 2 else None))
            if cohort in ("cys0_toxic_strong", "cys7_nontoxic_strong", "weaker_reverse"):
                for family, part in subset.groupby(subset.family.fillna("unannotated")):
                    outcomes = data.loc[part.index, f"{model}_{mode}_outcome"]
                    family_rows.append(dict(model=model, mode=mode, cohort=cohort, family=family,
                                            n=len(part), rescued=int(outcomes.eq("rescued").sum()),
                                            added_error=int(outcomes.eq("added_error").sum()),
                                            folds_with_rescue=part.loc[outcomes.eq("rescued"), "fold"].nunique()))
    return pd.DataFrame(summary), pd.DataFrame(family_rows)


def main():
    spec = importlib.util.spec_from_file_location("nested_compare", PACKAGE / "compare.py")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    cols = KEYS + ["cys", "length", "p_1d", "p_1d_raw_mean", "band", "high_seed_disagreement",
                   "family", "fragment", "conflict", "uniform_flag", "evidence_tier"]
    data = pd.read_csv(PROFILE, usecols=cols).set_index("sample_id", drop=False)
    assert len(data) == 6387 and data.index.is_unique and data.groupby("group_id").fold.nunique().eq(1).all()
    data.index.name = None
    manifest = pd.read_csv(SPLIT)
    bridge = pd.read_csv(BRIDGE)
    cross = bridge.set_index("historical_sample_id").loc[data.index]
    for key in ("sequence", "label", "fold", "group_id"):
        np.testing.assert_array_equal(cross[key], data[key])
    structure_file = STRUCTURE / "esm_if_probe_predictions.csv"
    struct = pd.read_csv(structure_file).set_index("sample_id").loc[data.index]
    assert struct.index.is_unique
    for key in ("sequence", "label", "fold", "group_id", "cys"):
        np.testing.assert_array_equal(struct[key], data[key])
    data["egnn_existing"] = struct.p_old3d
    features, historical, paths = sequence_features(data, bridge)
    for key, value in historical.items():
        data[key] = value
    audits = fit_probes(data, manifest, features, prior)
    models = [(key, key, None, "p50") for key in ["p_1d", "egnn_existing", *historical, *features]]
    models += [(key, key, f"{key}_matched_decision", "cal_matched_fp") for key in ["egnn_existing", *features]]
    summary, family = summarize(data, models)
    data.to_csv(OUT / "paired_samples.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "cohort_fold_summary.csv", index=False, encoding="utf-8-sig")
    family.to_csv(OUT / "family_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(audits).to_csv(OUT / "calibration_audit.csv", index=False, encoding="utf-8-sig")
    inputs = [PROFILE, SPLIT, BRIDGE, ESM, T5, IF1, structure_file, *paths]
    hashes = {str(p): sha(p) for p in inputs}
    known_esm = json.loads((STRUCTURE / "geometry_provenance.json").read_text())["sha256"]
    known_if1 = json.loads((STRUCTURE / "esm_if_features_manifest.json").read_text())["features_sha256"]
    assert hashes[str(ESM)] == known_esm and hashes[str(IF1)] == known_if1
    (OUT / "provenance.json").write_text(json.dumps({"samples": len(data), "inputs": hashes,
        "head": "StandardScaler + logistic regression C=0.1", "test1_test2_used": False}, indent=2))
    print(summary[summary.fold.eq("all") & summary.cohort.isin(
        ["all", "cys0_toxic_strong", "cys7_nontoxic_strong", "matched_cys0_nontoxic"]
    )].to_string(index=False))


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
