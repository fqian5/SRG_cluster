#!/bin/bash

# Submit all SRG jobs to the cluster

echo "Submitting all SRG jobs..."
echo "=========================="

# Submit H2 job
echo "Submitting H2 job..."
sbatch submit_h2.sh

# Submit all 10 random Hamiltonian jobs
for i in {0..9}; do
    echo "Submitting random Hamiltonian job $i..."
    sbatch submit_random_${i}.sh
done

echo "=========================="
echo "All jobs submitted!"
echo ""
echo "To check job status, use: squeue -u \$USER"
echo "To cancel all jobs, use: scancel -u \$USER"
