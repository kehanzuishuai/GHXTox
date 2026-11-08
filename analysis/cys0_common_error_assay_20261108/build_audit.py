"""Join frozen GHXTox scores to manually adjudicated experimental evidence."""
from collections import Counter
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT.parent
AUDIT = ANALYSIS / "error_evidence_linkage_20261025/all_training_evidence_index.csv"
ASSAYS = ANALYSIS / "experimental_dataset_inventory_20261031/curated/experimental_records_all.csv"
MODEL_TABLE = ANALYSIS / "sequence_complementarity_20261108/paired_samples.csv"
PROFILE = ANALYSIS / "six_seed_1d_weakness_20261107/sample_profiles.csv"
MODELS = {
    "egnn_existing": "EGNN", "esm2_mean": "ESM-2",
    "pepbert_large_historical": "PepBERT-L", "pepbert_small_historical": "PepBERT-S",
    "esmif_mean": "ESM-IF1", "esm2_esmif_mean": "ESM-2+ESM-IF1",
}
LABELS = {
    "damage": "直接宿主损伤实验", "specific": "特定功能或条件效应",
    "fragment": "部分序列／对象不符", "precursor": "前体／成熟体不符",
    "analogue": "天然体／实测类似物不符", "identity": "仅身份或直接实验未核实", "limited_negative": "条件受限的弱或阴性结果",
    "unresolved": "未能确认身份",
}


def reference_urls(raw):
    return [
        f"https://pubmed.ncbi.nlm.nih.gov/{item}/" if item.isdigit() else item
        for item in raw.split(";") if item
    ]


def read_fasta(path):
    for record in path.read_text(encoding="utf-8-sig").split(">")[1:]:
        header, *sequence = record.splitlines()
        yield header.split()[0], "".join(sequence).strip()


def upstream_membership(target):
    sources = OUT / "sources"
    files = {
        "csm_positive": "csm_uniprot_positive_curated.fasta",
        "csm_blind_positive": "csm_uniprot_blind_positive_curated.fasta",
        "csm_positive_cdhit70": "csm_uniprot_positive_curated_cdhit70.fasta",
        "csm_blind_positive_cdhit70": "csm_uniprot_blind_positive_curated_cdhit70.fasta",
        "toxinpred2_positive": "toxinpred2_positive.fasta",
    }
    fasta = {key: list(read_fasta(sources / file)) for key, file in files.items()}
    atse = {split: pd.read_csv(sources / f"atse_{split}.csv") for split in ("train", "test")}
    assert len(atse["train"]) == 3284 and len(atse["test"]) == 580
    assert len(fasta["toxinpred2_positive"]) == 8233
    all_matches = {}
    for row in target.itertuples():
        entry = {key: [identifier for identifier, seq in rows if seq == row.sequence]
                 for key, rows in fasta.items()}
        for split, frame in atse.items():
            entry[f"atse_{split}_labels"] = frame.loc[frame.sequence.eq(row.sequence), "label"].astype(int).tolist()
        all_matches[row.sample_id] = entry
    return all_matches


def main():
    reviews = {}
    for line in (OUT / "manual_evidence.psv").read_text(encoding="utf-8-sig").splitlines():
        number, category, citations, note = line.split("|", 3)
        sid = f"peptide_{number}"
        assert sid not in reviews and category in LABELS
        reviews[sid] = (category, reference_urls(citations), note)
    data = pd.read_csv(MODEL_TABLE)
    target = data[data.cohort_cys0_toxic_strong].copy()
    assert len(data) == 6387 and data.sample_id.is_unique
    assert len(target) == len(reviews) == 67 and set(target.sample_id) == set(reviews)
    assert target.p_1d.le(.10).all() and target.label.eq(1).all() and target.cys.eq(0).all()
    seq = pd.read_csv(PROFILE).set_index("sample_id")
    assert seq.loc[target.sample_id, "seeds_support_label_side"].eq(0).all()
    audit = pd.read_csv(AUDIT).set_index("sample_id")
    archive = pd.read_csv(ASSAYS, low_memory=False)
    measurements = archive[archive.sequence.isin(target.sequence)].groupby("sequence")
    measurement_count = measurements.size().to_dict()
    uniprot = json.loads((OUT / "sources/uniprot_identity_snapshot.json").read_text(encoding="utf-8"))
    membership = upstream_membership(target)
    records = []
    for r in target.sort_values("sample_id", key=lambda s: s.str.split("_").str[-1].astype(int)).itertuples():
        category, urls, note = reviews[r.sample_id]
        a = audit.loc[r.sample_id]
        assert a.sequence == r.sequence and a.original_label == r.label and a.fold == r.fold
        accessions = [] if pd.isna(a.automatic_annotation_accessions) else a.automatic_annotation_accessions.split(";")
        for accession in accessions:
            assert r.sequence in uniprot[accession]["sequence"]["value"]
        calls = {model: float(getattr(r, model)) for model in MODELS}
        exceptions = [MODELS[model] for model, probability in calls.items() if probability >= .5]
        records.append(dict(
            sample_id=r.sample_id, sequence=r.sequence, fold=int(r.fold), group_id=r.group_id,
            p_prott5_six=float(r.p_1d), alternative_probabilities=calls, alternatives_calling_toxic=exceptions,
            all_models_call_nontoxic=not exceptions, fragment_flag=int(r.fragment),
            existing_annotation_accessions=accessions, prior_name=None if pd.isna(a["name"]) else a["name"],
            exact_sequence_archived_assays=measurement_count.get(r.sequence, 0),
            category=category, original_label_trace="unresolved", note=note, references=urls,
            upstream_exact_membership=membership[r.sample_id],
        ))
    assert len({r["sample_id"] for r in records}) == 67
    assert sum(r["all_models_call_nontoxic"] for r in records) == 62
    assert sum(r["fragment_flag"] for r in records) == 22
    assert {r["sample_id"] for r in records if r["fragment_flag"]} == {
        r["sample_id"] for r in records if r["category"] == "fragment"
    }
    groups = Counter(r["category"] for r in records)
    shared = Counter(r["category"] for r in records if r["all_models_call_nontoxic"])
    matched_by_source = {
        name: sum(bool(membership[sid][name]) for sid in membership)
        for name in next(iter(membership.values()))
    }
    matched_any = sum(any(r["upstream_exact_membership"].values()) for r in records)
    result = {
        "rows": len(records),
        "groups": dict(groups),
        "common_errors": sum(shared.values()),
        "common_error_groups": dict(shared),
        "exact_sequence_archived_assay_samples": sum(r["exact_sequence_archived_assays"] > 0 for r in records),
        "exact_sequence_archived_assay_rows": sum(r["exact_sequence_archived_assays"] for r in records),
        "predictions_use_test1_test2": False,
        "upstream_exact_membership_counts": matched_by_source,
        "upstream_any_exact_membership_count": matched_any,
        "records": records,
    }
    (OUT / "audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# 67 条 Cys=0 强反向肽的逐样本实验溯源", "",
        "替代模型在 0.5 判有毒的列为空时，该样本在所列六个替代模型与六种子 ProtT5 基线下共同漏判。"
        "这里的漏判相对于数据集原标签。每条只分配一个首要证据类型，化学修饰和实验终点仍可在备注内同时记录。",
        "“旧实验记录”仅计精确裸序列匹配，不能核定实测分子、结果或原标签来源。上游匹配只说明该序列出现过，具体来源见 audit.json。", "",
        "| 样本 | 序列 | 折 | 六种子有毒概率 | 判有毒的替代模型 | 旧实验记录 | 证据类型 | 对象、终点和证据 |",
        "|---|---|---:|---:|---|---:|---|---|",
    ]
    for r in records:
        refs = " ".join(f"[来源 {i}]({url})" for i, url in enumerate(r["references"], 1))
        if not refs and r["existing_annotation_accessions"]:
            refs = " ".join(f"[{a}](https://www.uniprot.org/uniprotkb/{a}/entry)"
                            for a in r["existing_annotation_accessions"])
        alt = "、".join(r["alternatives_calling_toxic"]) or "—"
        lines.append(
            f"| {r['sample_id']} | {r['sequence']} | {r['fold']} | {r['p_prott5_six']:.4f} | {alt} | "
            f"{r['exact_sequence_archived_assays']} | {LABELS[r['category']]} | {r['note']} {refs} |"
        )
    (OUT / "逐样本核查.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print({key: value for key, value in result.items() if key != "records"})


if __name__ == "__main__":
    main()
