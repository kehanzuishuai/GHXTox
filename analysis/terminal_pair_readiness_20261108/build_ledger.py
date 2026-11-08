"""Build the terminal chemistry diagnostic ledger from audited records."""

import json
from collections import Counter
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
PREVIOUS = ANALYSIS / "label_ambiguity_map_20261025" / "public_terminal_pair_blocks.json"
GRAMMISTIN = "https://link.springer.com/article/10.1111/j.1444-2906.2006.01264.x"
GRAMMISTIN_SEQUENCES = "https://www.jstage.jst.go.jp/article/fishsci1994/67/1/67_1_163/_pdf"

SOURCES = {
    "17980935": ("BP76/BP77", "BP76", "https://www.academia.edu/18839892/A_library_of_linear_undecapeptides_with_bactericidal_activity_against_phytopathogenic_bacteria", "Table 1 and hemolysis methods"),
    "10879468": ("Mast 21", "Mast 21", "https://www.jstage.jst.go.jp/article/bbb/64/5/64_5_985/_pdf", "Tables 1–2 and hemolysis methods"),
    "28230084": ("MreB1–9 series", "MreB1–9 series", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5322399/", "Tables 1 and 5 and hemolysis methods"),
    "29108968": ("Pteroicidin alpha", "Pteroicidin alpha", "https://sygal.normandie-univ.fr/informations/fichiers/telecharger/303", "Primary article reproduced in thesis, section 3.4 and Figure 4"),
    "34361810": ("Kn2-5R", "Kn2-5R", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8348142/", "Table 4 and hemolysis methods"),
    "34819923": ("RR4", "RR4", "https://doi.org/10.3389/fmicb.2021.753874", "Supplementary Tables 1 and 4 and hemolysis methods"),
}

GRAMMISTIN_PAIRS = (
    ("Gs1", "LFGFLIKLIPSLFGALSNIGRNRNQ", [], "no_observed_difference"),
    ("Gs2", "LFGFLIPLLPHIIGAIPQVIGAIR", ["Pp4a"], "amidated_higher"),
    ("Pp1", "FIGGIISFFKRLF", ["GsD"], "amidated_higher"),
    ("Pp2b", "FIGGIISFIKKLF", ["GsE"], "amidated_higher"),
    ("Pp3", "NWRKILGQIASVGAGLLGSLLAGYE", [], "amidated_higher"),
)


def add_previous_pairs(pairs):
    blocks = json.loads(PREVIOUS.read_text(encoding="utf-8"))
    for block in blocks:
        sequence, pmid = block["sequence"], block["pmid"]
        name, family, source, locator = SOURCES[pmid]
        records = sorted(
            block["records"], key=lambda row: (row["n_terminal"] != "free", row["c_terminal"] != "acid")
        )
        assert len(records) == 2
        states = [(row["n_terminal"], row["c_terminal"]) for row in records]
        is_figure = pmid == "29108968"
        pair = pairs.setdefault(sequence, {
            "name": records[0]["name"] if pmid == "28230084" else name,
            "sequence": sequence,
            "group": family,
            "source": source,
            "source_locator": locator,
            "pmid": pmid,
            "tier": "figure" if is_figure else "quantitative",
            "forms": [{"n_terminal": n, "c_terminal": c} for n, c in states],
            "rbc_species": ["human", "sea bass", "lesser-spotted dogfish"] if is_figure else [block["species"]],
            "matched_percent_hemolysis": [],
            "legacy_record_ids": [],
        })
        assert pair["forms"] == [{"n_terminal": n, "c_terminal": c} for n, c in states]
        record_ids = [row["id"] for row in records]
        pair["legacy_record_ids"].append(record_ids)
        if is_figure:
            pair["evidence_note"] = "Do not reuse legacy fish 50% vs 90% as an exact species-matched pair. Reported 90% is an upper value across conditions."
            continue
        dose = block["dose_uM"]
        measurement = {
            "dose_uM": dose,
            "percent_by_form": [row["hemolysis_percent"] for row in records],
            "legacy_record_ids": record_ids,
        }
        if pmid == "17980935" and record_ids == ["1843", "1847"]:
            measurement.update(dose_uM=150.0, legacy_dose_uM=dose)
        pair["matched_percent_hemolysis"].append(measurement)
    assert len(blocks) == 19 and len(pairs) == 7
    rr4 = pairs["WLRRIKAWLRRIKA"]
    rr4["matched_percent_hemolysis"].insert(1, {
        "dose_uM": 250.0, "percent_by_form": [9.8, 85.0], "source_locator": "Supplementary Table 4",
    })


def main():
    pairs = {}
    add_previous_pairs(pairs)
    for name, sequence, aliases, result in GRAMMISTIN_PAIRS:
        assert sequence not in pairs
        pairs[sequence] = {
            "name": name, "aliases": aliases, "sequence": sequence, "group": "Grammistin 2006 study",
            "source": GRAMMISTIN, "sequence_source": GRAMMISTIN_SEQUENCES,
            "source_locator": "Abstract", "tier": "abstract",
            "forms": [{"n_terminal": None, "c_terminal": c} for c in ("acid", "amide")],
            "rbc_species": [], "reported_activity_direction": result, "matched_percent_hemolysis": [],
        }
    for sequence, pair in pairs.items():
        first, second = pair["forms"]
        assert set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
        assert sum(first[key] != second[key] for key in first) == 1
        assert len({(sequence, form["n_terminal"], form["c_terminal"]) for form in pair["forms"]}) == 2
        measurements = pair["matched_percent_hemolysis"]
        assert len({row["dose_uM"] for row in measurements}) == len(measurements)
    numeric = [pair for pair in pairs.values() if pair["tier"] == "quantitative"]
    human = [pair for pair in numeric if pair["rbc_species"] == ["human"]]
    summary = {
        "distinct_bare_sequences": len(pairs),
        "backbones_by_tier": dict(Counter(pair["tier"] for pair in pairs.values())),
        "quantitative_studies": len({pair["pmid"] for pair in numeric}),
        "quantitative_dose_blocks": sum(len(pair["matched_percent_hemolysis"]) for pair in numeric),
        "quantitative_human_backbones": len(human),
        "quantitative_human_studies": len({pair["pmid"] for pair in human}),
        "quantitative_human_c_terminal_backbones": sum(
            pair["forms"][0]["c_terminal"] != pair["forms"][1]["c_terminal"] for pair in human
        ),
    }
    (HERE / "pair_ledger.json").write_text(
        json.dumps({"summary": summary, "pairs": list(pairs.values())}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
