#!/usr/bin/env python3
"""Render the zero surfaces of DeepSDF training instances as standalone HTML."""

import argparse
import importlib
import json
from pathlib import Path
import sys

import numpy as np
from skimage.measure import marching_cubes
import torch
import trimesh
from trimesh.viewer.notebook import scene_to_html


PROJECT = Path(__file__).resolve().parent
DEEPSDF = PROJECT / "DeepSDF"


def load_experiment(experiment, checkpoint, device):
    specs = json.loads((experiment / "specs.json").read_text())
    split_path = DEEPSDF / specs["TrainSplit"]
    split = json.loads(split_path.read_text())
    instances = [
        Path(dataset) / class_name / name
        for dataset, classes in split.items()
        for class_name, names in classes.items()
        for name in names
    ]

    sys.path.insert(0, str(DEEPSDF))
    arch = importlib.import_module("networks." + specs["NetworkArch"])
    decoder = arch.Decoder(specs["CodeLength"], **specs["NetworkSpecs"])
    filename = f"{checkpoint}.pth"
    model = torch.load(
        experiment / "ModelParameters" / filename,
        map_location="cpu",
        weights_only=True,
    )
    latent = torch.load(
        experiment / "LatentCodes" / filename,
        map_location="cpu",
        weights_only=True,
    )
    if model["epoch"] != latent["epoch"]:
        raise ValueError("Model and latent codes are from different epochs.")

    decoder.load_state_dict({
        key.removeprefix("module."): value
        for key, value in model["model_state_dict"].items()
    })
    codes = latent["latent_codes"]
    if isinstance(codes, dict):
        codes = codes["weight"]
    if codes.ndim == 3 and codes.shape[1] == 1:
        codes = codes.squeeze(1)
    if codes.shape != (len(instances), specs["CodeLength"]):
        raise ValueError("Latent codes do not match the training split or code length.")

    return decoder.to(device).eval(), codes.to(device), instances, model["epoch"]


@torch.inference_mode()
def sample_field(decoder, code, resolution, batch_size=32768):
    """Sample an (x, y, z) volume over the normalized cube [-1, 1]^3."""
    axis = torch.linspace(-1, 1, resolution, device=code.device, dtype=code.dtype)
    count = resolution ** 3
    plane_size = resolution ** 2
    values = np.empty(count, dtype=np.float32)
    for start in range(0, count, batch_size):
        end = min(start + batch_size, count)
        index = torch.arange(start, end, device=code.device)
        xyz = torch.stack((
            axis[index // plane_size],
            axis[(index // resolution) % resolution],
            axis[index % resolution],
        ), dim=1)
        inputs = torch.cat((code.expand(len(xyz), -1), xyz), dim=1)
        values[start:end] = decoder(inputs).flatten().cpu().numpy()
    return values.reshape((resolution,) * 3)


def extract_zero_surface(field):
    if not np.isfinite(field).all():
        raise ValueError("Predicted SDF contains non-finite values.")
    if not field.min() < 0 < field.max():
        return None

    step = 2.0 / (field.shape[0] - 1)
    vertices, faces, _, _ = marching_cubes(
        field,
        level=0.0,
        spacing=(step, step, step),
        gradient_direction="ascent",
        allow_degenerate=False,
    )
    mesh = trimesh.Trimesh(vertices=vertices - 1.0, faces=faces, process=False)
    mesh.fix_normals()
    mesh.visual.face_colors = [71, 155, 207, 255]
    return mesh


def write_html(mesh, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(scene_to_html(mesh.scene()), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-e", "--experiment", type=Path, default=PROJECT / "experiments" / "smoke",
        help="Experiment containing specs.json, ModelParameters and LatentCodes.",
    )
    parser.add_argument(
        "-c", "--checkpoint", default="latest", help="Epoch number or latest."
    )
    parser.add_argument(
        "--resolution", type=int, default=128, help="Grid samples per axis (default: 128)."
    )
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    if args.checkpoint != "latest" and not args.checkpoint.isdecimal():
        parser.error("--checkpoint must be an epoch number or latest")
    if args.resolution < 2:
        parser.error("--resolution must be at least 2")
    if args.device == "cuda" and not torch.cuda.is_available():
        parser.error("--device cuda was requested, but CUDA is unavailable")

    decoder, codes, instances, epoch = load_experiment(
        args.experiment, args.checkpoint, args.device
    )
    output_dir = args.experiment / "Visualizations" / str(epoch)
    written = 0
    for instance, code in zip(instances, codes):
        field = sample_field(decoder, code, args.resolution)
        print(f"{instance}: SDF range [{field.min():.6g}, {field.max():.6g}]")
        mesh = extract_zero_surface(field)
        if mesh is None:
            print("  No sign change on this grid. No zero surface was exported.")
            continue

        output = (output_dir / instance).with_suffix(".html")
        write_html(mesh, output)
        print(f"  {len(mesh.vertices)} vertices, {len(mesh.faces)} faces, "
              f"watertight={mesh.is_watertight}")
        print(f"  {output}")
        written += 1

    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main())
