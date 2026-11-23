import json
import os
from pathlib import Path
import numpy as np
from SdfSampler.sdf_sampler import SDFSampler, SDFBase

root = Path(os.environ["NTE_ROOT"])
data = root / "data" / "smoke"
exp = root / "experiments" / "smoke"
exp.mkdir(parents=True, exist_ok=True)

class SphereSDF(SDFBase):
    def __init__(self, radius):
        self.radius = radius
    def __call__(self, queries):
        return (np.linalg.norm(queries, axis=1) - self.radius).reshape(-1, 1)

np.random.seed(42)
sampler = SDFSampler(
    outdir=str(data / "SdfSamples"),
    splitdir=str(exp),
)
info = {"dataset_name": "Smoke", "class_name": "sphere"}
ids = sampler.sample_sdfs(
    [SphereSDF(0.4), SphereSDF(0.6)],
    info, n_samples=20000,
)
sampler.write_json("train.json", info, ids)

for p in data.rglob("*.npz"):
    with np.load(p) as arrays:
        for key in ("pos", "neg"):
            a = arrays[key]
            assert a.ndim == 2 and a.shape[1] == 4
            assert len(a) > 0 and np.isfinite(a).all()

spec = json.loads(
    (root / "DeepSDF/examples/sofas/specs.json").read_text()
)
spec.update({
    "Description": "Installation smoke test on two spheres",
    "DataSource": str(data),
    "TrainSplit": str(exp / "train.json"),
    "CodeLength": 16,
    "NetworkSpecs": {"dims": [64, 64, 64]},
    "NumEpochs": 2,
    "SnapshotFrequency": 1,
    "AdditionalSnapshots": [],
    "LogFrequency": 1,
    "SamplesPerScene": 1024,
    "ScenesPerBatch": 2,
    "DataLoaderThreads": 0,
})
spec.pop("TestSplit", None)
(exp / "specs.json").write_text(json.dumps(spec, indent=2))
print("Smoke experiment:", exp)