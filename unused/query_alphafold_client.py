#!/usr/bin/env python3

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Optional

import requests

# AlphaFold Server API endpoint
DEFAULT_ALPHAFOLD_API_URL = "https://alphafoldserver.com/api/v1"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit a protein sequence to AlphaFold Server and retrieve predictions (no templates)."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the input FASTA file containing the protein sequence.",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Path to the output ZIP archive containing the prediction results.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=120,
        help="Maximum number of retry attempts (default: 120, ~30 min with 15s interval).",
    )
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=15.0,
        help="Polling interval in seconds (default: 15).",
    )
    parser.add_argument(
        "--api-url",
        default=os.environ.get("ALPHAFOLD_API_URL", DEFAULT_ALPHAFOLD_API_URL),
        help="Base URL for the AlphaFold Server API (defaults to ALPHAFOLD_API_URL or the public endpoint).",
    )
    return parser.parse_args()


def read_fasta(fasta_path: Path) -> str:
    """Read the first sequence from a FASTA file."""
    sequence = ""
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    break
            else:
                sequence += line
    if not sequence:
        raise ValueError(f"No sequence found in {fasta_path}")
    return sequence


def submit_prediction(sequence: str, api_url: str) -> str:
    """Submit a sequence to AlphaFold Server and return the job ID."""
    payload = {
        "sequences": [{"sequence": sequence}],
        "template_mode": "none",  # No templates
        "multimer": False,
    }

    headers = {"Content-Type": "application/json"}

    try:
        response = requests.post(
            f"{api_url}/predict",
            json=payload,
            headers=headers,
            timeout=30,
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(
            f"Failed to submit prediction to {api_url}/predict: {exc}. "
            "The endpoint may be unavailable, require authentication, or use a different API path."
        )

    result = response.json()
    if "job_id" not in result:
        raise RuntimeError(f"Unexpected response from AlphaFold Server: {result}")

    job_id = result["job_id"]
    print(f"Submitted prediction job: {job_id}")
    return job_id


def check_job_status(job_id: str, api_url: str) -> dict:
    """Check the status of a prediction job."""
    try:
        response = requests.get(
            f"{api_url}/predict/{job_id}",
            timeout=30,
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(
            f"Failed to check job status for {job_id} at {api_url}/predict/{job_id}: {exc}"
        )

    return response.json()


def download_results(job_id: str, output_zip: Path, api_url: str) -> None:
    """Download the prediction results as a ZIP archive."""
    try:
        response = requests.get(
            f"{api_url}/predict/{job_id}/download",
            timeout=120,
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(
            f"Failed to download results for {job_id} from {api_url}/predict/{job_id}/download: {exc}"
        )

    output_zip.parent.mkdir(parents=True, exist_ok=True)
    with open(output_zip, "wb") as f:
        f.write(response.content)

    print(f"Downloaded results to: {output_zip}")


def wait_for_completion(
    job_id: str,
    max_retries: int = 120,
    poll_interval: float = 15.0,
    api_url: str = DEFAULT_ALPHAFOLD_API_URL,
) -> None:
    """Poll the job status until completion."""
    for attempt in range(max_retries):
        status_data = check_job_status(job_id, api_url)
        status = status_data.get("status", "unknown").lower()

        print(f"[Attempt {attempt + 1}/{max_retries}] Job status: {status}")

        if status == "completed":
            print("Prediction completed successfully!")
            return
        elif status in ("failed", "error"):
            raise RuntimeError(f"Prediction job failed: {status_data}")
        elif status == "running":
            print(f"Job still running. Waiting {poll_interval}s before next check...")
            time.sleep(poll_interval)
        else:
            print(f"Unknown status: {status}. Waiting...")
            time.sleep(poll_interval)

    raise TimeoutError(
        f"Job {job_id} did not complete within {max_retries * poll_interval}s"
    )


def main():
    args = parse_args()
    input_fasta = Path(args.input)
    output_zip = Path(args.output)
    api_url = args.api_url.rstrip("/")

    if not input_fasta.exists():
        print(f"Error: Input FASTA file not found: {input_fasta}", file=sys.stderr)
        sys.exit(1)

    try:
        # Read the sequence
        sequence = read_fasta(input_fasta)
        print(f"Loaded sequence of length {len(sequence)} from {input_fasta}")

        print(f"Using AlphaFold API base URL: {api_url}")

        # Submit the prediction
        job_id = submit_prediction(sequence, api_url)

        # Wait for completion
        wait_for_completion(job_id, args.max_retries, args.poll_interval, api_url)

        # Download the results
        download_results(job_id, output_zip, api_url)
        print(f"AlphaFold prediction completed. Results saved to: {output_zip}")

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
