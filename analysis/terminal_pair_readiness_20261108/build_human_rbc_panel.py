"""Assemble measured human RBC terminal comparisons as paired assay records."""

import csv
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPLIT = HERE.parent / "reviewer_response_20261102/six_seed_nested_autodl_20261106/data/nested_folds.csv"


def read_json(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def measured_pairs(legacy, new):
    for peptide in legacy["pairs"]:
        if peptide["tier"] != "quantitative" or peptide["rbc_species"] != ["human"]:
            continue
        buffers = {
            "17980935": "Tris saline, pH 7.2",
            "10879468": None,
            "28230084": "HEPES saline, pH 7.4",
        }
        pmid = peptide["pmid"]
        for point in peptide["matched_percent_hemolysis"]:
            yield {
                "study_id": pmid, "source": peptide["source"], "locator": peptide["source_locator"],
                "series": "1", "name": peptide["name"], "sequence": peptide["sequence"],
                "assay": {"species": "human", "temperature_c": 37, "minutes": 60,
                          "buffer": buffers[pmid], "concentration_unit": "uM"},
                "concentration": point["dose_uM"],
                "forms": [dict(form, hemolysis_percent=value) for form, value in zip(
                    peptide["forms"], point["percent_by_form"], strict=True
                )],
                "origin": "previously_audited",
            }

    for study in new["studies"]:
        for peptide in study["pairs"]:
            for dose, acid, amide, *errors in peptide["points"]:
                forms = [
                    {"n_terminal": "free", "c_terminal": terminal, "hemolysis_percent": value}
                    for terminal, value in (("acid", acid), ("amide", amide))
                ]
                if errors:
                    for form, error in zip(forms, errors, strict=True):
                        form["reported_error"] = error
                yield {
                    "study_id": study["id"], "source": study["source"], "locator": study["locator"],
                    "series": peptide["series"], "name": peptide["name"], "sequence": peptide["sequence"],
                    "assay": study["assay"], "concentration": dose, "forms": forms,
                    "origin": "new_primary_source",
                }


def main():
    new = read_json("human_rbc_primary_sources.json")
    pairs = list(measured_pairs(read_json("pair_ledger.json"), new))
    identities = set()
    counts = Counter()
    sequences = {row["sequence"] for row in pairs}
    with SPLIT.open(encoding="utf-8-sig", newline="") as handle:
        overlaps = {
            (row["sequence"], row["sample_id"])
            for row in csv.DictReader(handle)
            if row["sequence"] in sequences
        }
    for row in pairs:
        a, b = row["forms"]
        changed, = [field for field in ("n_terminal", "c_terminal") if a[field] != b[field]]
        key = (row["study_id"], row["series"], row["sequence"], changed,
               row["concentration"], row["assay"]["concentration_unit"])
        assert key not in identities
        identities.add(key)
        assert row["assay"]["species"] == "human" and row["concentration"] > 0
        assert all(0 <= form["hemolysis_percent"] <= 100 for form in (a, b))
        assert set(row["sequence"]) <= set("ACDEFGHIKLMNPQRSTVWY")
        diff = Decimal(str(b["hemolysis_percent"])) - Decimal(str(a["hemolysis_percent"]))
        row["changed_terminal"] = changed
        row["delta_percentage_points"] = float(diff)
        counts[changed] += 1
    summary = {
        "studies": len({row["study_id"] for row in pairs}),
        "unique_sequences": len(sequences),
        "dose_pairs": len(pairs),
        "dose_pairs_by_terminal": dict(counts),
        "c_terminal_sequences": len({row["sequence"] for row in pairs if row["changed_terminal"] == "c_terminal"}),
        "new_quantitative_sequences": len({row["sequence"] for row in pairs if row["origin"] == "new_primary_source"}),
        "new_quantitative_dose_pairs": sum(row["origin"] == "new_primary_source" for row in pairs),
        "six_seed_training_manifest_exact_overlaps": [
            {"sequence": sequence, "sample_id": sample_id} for sequence, sample_id in sorted(overlaps)
        ],
    }
    result = {"summary": summary, "pairs": pairs, "other_evidence": new["other_evidence"]}
    (HERE / "human_rbc_matched_pairs.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
