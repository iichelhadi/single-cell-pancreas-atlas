#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import time

import pandas as pd


EXPR_FILE = "expression_matrix.csv"
TF_FILE = "1.1_inputTFs.txt"
OUT_DIR = "grn_runs"
TIMING_FILE = os.path.join(OUT_DIR, "timing.tsv")

REQUIRED_COLUMNS = ["TF", "target", "importance"]


def validate_inputs() -> str:
    executable = shutil.which("arboreto_with_multiprocessing.py")

    if executable is None:
        raise FileNotFoundError(
            "arboreto_with_multiprocessing.py was not found in the active environment."
        )

    for path in [EXPR_FILE, TF_FILE]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing required input: {path}")

    return executable


def validate_output(path: str) -> None:
    output = pd.read_csv(path, sep="\t")

    if not all(column in output.columns for column in REQUIRED_COLUMNS):
        raise ValueError(
            f"Unexpected columns in {path}: {list(output.columns)}. "
            f"Expected {REQUIRED_COLUMNS}."
        )

    if output.empty:
        raise ValueError(f"{path} contains no GRN edges.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-runs", type=int, default=50)
    parser.add_argument("--seed-start", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--calibrate", action="store_true")
    args = parser.parse_args()

    executable = validate_inputs()
    os.makedirs(OUT_DIR, exist_ok=True)

    number_of_runs = 1 if args.calibrate else args.n_runs
    timing_rows = []

    for index in range(number_of_runs):
        seed = args.seed_start + index
        output_file = os.path.join(
            OUT_DIR,
            f"grn_output_run{seed}.tsv"
        )

        if os.path.exists(output_file):
            print(f"Skipping existing output: {output_file}", flush=True)
            continue

        command = [
            executable,
            EXPR_FILE,
            TF_FILE,
            "--method", "grnboost2",
            "--output", output_file,
            "--num_workers", str(args.workers),
            "--seed", str(seed),
        ]

        print(
            f"Starting run {index + 1}/{number_of_runs}, seed={seed}",
            flush=True
        )
        print("Command:", " ".join(command), flush=True)

        start = time.time()
        subprocess.run(command, check=True)
        elapsed = time.time() - start

        validate_output(output_file)

        timing_rows.append({
            "seed": seed,
            "seconds": round(elapsed, 1),
            "workers": args.workers,
            "output": output_file,
        })

        print(
            f"Completed seed={seed} in {elapsed / 3600:.2f} hours",
            flush=True
        )

    if timing_rows:
        new_timings = pd.DataFrame(timing_rows)

        if os.path.exists(TIMING_FILE):
            old_timings = pd.read_csv(TIMING_FILE, sep="\t")
            new_timings = pd.concat(
                [old_timings, new_timings],
                ignore_index=True
            ).drop_duplicates(subset=["seed"], keep="last")

        new_timings.to_csv(
            TIMING_FILE,
            sep="\t",
            index=False
        )

    if args.calibrate and timing_rows:
        seconds = timing_rows[0]["seconds"]

        print(
            "\nCALIBRATION RESULT\n"
            f"One run using {args.workers} workers: "
            f"{seconds / 3600:.2f} hours\n"
            f"Estimated 20 sequential runs: "
            f"{seconds * 20 / 3600:.1f} hours\n"
            f"Estimated 50 sequential runs: "
            f"{seconds * 50 / 3600:.1f} hours\n"
            f"Estimated 100 sequential runs: "
            f"{seconds * 100 / 3600:.1f} hours",
            flush=True
        )


if __name__ == "__main__":
    main()

