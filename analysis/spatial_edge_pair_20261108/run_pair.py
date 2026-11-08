"""Fixed fold-0 pilot comparing equal-budget spatial and sequential extra edges."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import torch

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
REVIEW = ROOT / "analysis/reviewer_response_20261102"
PREP = REVIEW / "paper_ablation_prepared"
CACHE = Path("D:/dachuang2026/GHXTox-ziliao/train_chemical_sites_final_esm2_multisource_1_30_full2795.pt")
MANIFEST = PREP / "nested_folds.csv"
MODES = ("matched_sequence", "matched_spatial")


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else streaming_sha(stream)


def streaming_sha(stream):
    digest = hashlib.sha256()
    for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()


def verify_graphs():
    from ghxtox.models.layers import PLDDTAwareEGNNLayer

    generator = torch.Generator().manual_seed(42)
    coords = torch.randn(5, 20, 3, generator=generator)
    valid = torch.arange(20)[None, :] < torch.tensor([1, 2, 5, 11, 20])[:, None]
    distance = torch.cdist(coords, coords)
    separation = (torch.arange(20)[:, None] - torch.arange(20)[None, :]).abs()
    pairs = valid[:, :, None] & valid[:, None, :]
    adjacent = pairs & (separation == 1)
    expected = adjacent.sum(dim=-1) + (pairs & (separation > 1)).sum(dim=-1).clamp(max=4)
    graphs = []
    for mode in MODES:
        layer = PLDDTAwareEGNNLayer(hidden_dim=16, rbf_bins=16, dropout=0, graph_mode=mode, spatial_top_k=4)
        graph = layer._graph_mask(distance, valid)
        assert torch.equal(graph.sum(dim=-1), expected)
        assert (graph | adjacent).equal(graph)
        assert not (graph & (~pairs | (separation == 0))).any()
        graphs.append(graph)
    assert not torch.equal(*graphs)
    return {"row_edge_counts_equal": True, "changed_directed_edges": int((graphs[0] ^ graphs[1]).sum())}


def main():
    assert sha(CACHE) == "f8a6ddf7dd81198142ad8c30a0f7926de1967618b01720bd90cb483a45dfdb5b"
    assert sha(MANIFEST) == sha(REVIEW / "six_seed_nested_autodl_20261106/data/nested_folds.csv")
    base = json.loads((PREP / "A_backbone.json").read_text())
    base["model"]["spatial_top_k"] = 4
    base["train"].update(epochs=3, selection_mode="fixed_epoch", coord_noise_std=0)
    record = {"fold": 0, "epochs": 3, "seed": 42, "graph_preflight": verify_graphs(),
              "cache_sha256": sha(CACHE), "split_sha256": sha(MANIFEST), "status": "running"}
    run = OUT / "fold0"
    run.mkdir(exist_ok=False)
    metadata = run / "metadata.json"
    metadata.write_text(json.dumps(record, indent=2))
    for mode in MODES:
        base["model"]["graph_mode"] = mode
        config = run / f"{mode}.json"
        config.write_text(json.dumps(base, indent=2))
        output = run / mode
        common = ["--nested-manifest", str(MANIFEST), "--outer-fold", "0", "--device", "cuda"]
        subprocess.run([sys.executable, "-u", "-m", "ghxtox.train", "--train", str(CACHE),
                        "--config", str(config), "--output-dir", str(output), *common], cwd=ROOT, check=True)
        for role in ("calibration", "test"):
            subprocess.run([sys.executable, "-u", "-m", "ghxtox.evaluate",
                            "--checkpoint", str(output / "best_model.pt"), "--processed", str(CACHE),
                            "--nested-role", role, "--output", str(output / f"{role}_metrics.json"),
                            "--predictions", str(output / f"{role}_predictions.csv"), *common], cwd=ROOT, check=True)
    record["status"] = "complete"
    metadata.write_text(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
