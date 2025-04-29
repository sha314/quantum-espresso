#!/bin/bash


# .
# ├── input/                   # input files
# │   ├── scf.in
# │   ├── nscf.in
# │   ├── bands.in
# │   ├── dos.in
# ├── run_qe.sh                # << THIS SCRIPT
# └── tmp/                     # temporary output dir



# === CONFIGURATION ===
BIN_DIR="/opt/qe741/bin"                      # Path to QE binaries (pw.x, bands.x, etc.)
PSEUDO_DIR="$HOME/pseudo"                     # Path to pseudopotentials
INPUT_DIR="./input"
TMP_DIR="./tmp"
OUTPUT_DIR="./output"
NPROC=1                                      # Number of MPI processes. Keep it =1 for GPU computation

# === ENV SETUP ===
export ESPRESSO_PSEUDO=$PSEUDO_DIR
mkdir -p "$TMP_DIR" "$OUTPUT_DIR"


# === WORKFLOW ===
echo "SCF : Current time: $(date +"%T")"
mpirun -np $NPROC "$BIN_DIR/pw.x" -in "${INPUT_DIR}/scf.in" > "${OUTPUT_DIR}/scf.out"

echo "NSCF : Current time: $(date +"%T")"
mpirun -np $NPROC "$BIN_DIR/pw.x" -in "${INPUT_DIR}/nscf.in" > "${OUTPUT_DIR}/nscf.out"

echo "bands : Current time: $(date +"%T")"
mpirun -np $NPROC "$BIN_DIR/pw.x" -in "${INPUT_DIR}/bands.in" > "${OUTPUT_DIR}/bands.out"

# === POST-PROCESSING ===
# you use bands.pp.in with bands.x to post-process the data and extract eigenvalues, suitable for plotting.
echo "bands post processing : Current time: $(date +"%T")"
mpirun -np $NPROC "$BIN_DIR/bands.x" -in "${INPUT_DIR}/bands.pp.in" > "${OUTPUT_DIR}/bands_post.out"

echo "dos post processing : Current time: $(date +"%T")"
mpirun -np $NPROC "$BIN_DIR/dos.x" -in "${INPUT_DIR}/dos.in" > "${OUTPUT_DIR}/dos_post.out"

echo "Current time: $(date +"%T")"
echo ">>> All calculations complete!"