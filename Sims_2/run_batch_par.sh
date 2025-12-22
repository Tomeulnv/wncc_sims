#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

nohup R CMD BATCH InterMatching_Par.R &
sim_pid=$!
wait "$sim_pid"

R CMD BATCH display_results.R