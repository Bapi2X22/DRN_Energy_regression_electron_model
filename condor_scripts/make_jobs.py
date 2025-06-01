# make_jobs.py
import os
import math
import sys

input_file = sys.argv[1]   # input_files.txt
files_per_job = int(sys.argv[2])  # e.g. 5

with open(input_file) as f:
    all_files = [line.strip() for line in f if line.strip()]

n_jobs = math.ceil(len(all_files) / files_per_job)
os.makedirs("batch_jobs_electron", exist_ok=True)

for i in range(n_jobs):
    batch = all_files[i * files_per_job : (i + 1) * files_per_job]
    job_file = f"batch_jobs_electron/job_{i}.txt"
    with open(job_file, "w") as out:
        out.write("\n".join(batch))
