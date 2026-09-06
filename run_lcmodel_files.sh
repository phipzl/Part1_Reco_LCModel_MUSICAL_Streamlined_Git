#!/bin/bash
# Run LCModel on the input files a finished Part1 run left behind.
#
#   run_lcmodel_files.sh <output directory> [cores]
#
# Part1 writes one .RAW and one .control per voxel into <output directory>/spectra
# and keeps them when run with -d. This fits them the way step 7 does, so a run
# made with -l or -Q can be completed with LCModel afterwards and Part2 has its
# tables. The directory may have been copied or moved since: the control files
# name the directory they were written for, and are rewritten when it differs.
set -euo pipefail

out_path=$(cd "${1:?usage: run_lcmodel_files.sh <output directory> [cores]}" && pwd)
cores=${2:-$(nproc)}
spectra=$out_path/spectra

ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# It pulls in the MINC toolkit config, which reads variables that may be unset.
if [[ -f $ScriptDir/InstallProgramPaths.sh ]]; then
    set +u; source "$ScriptDir/InstallProgramPaths.sh" >/dev/null; set -u
fi
lcmodel=${LCM_Path:-$(command -v lcmodel || true)}
if [[ ! -x $lcmodel ]]; then
    echo "LCModel was not found. Set LCM_Path to the lcmodel binary." >&2
    exit 1
fi

n=$(find "$spectra" -maxdepth 1 -name '*.control' 2>/dev/null | wc -l)
if [[ $n -eq 0 ]]; then
    echo "No .control files in $spectra." >&2
    echo "Part1 deletes them at the end unless it was run with -d." >&2
    exit 1
fi

# The control files carry the absolute paths of the directory they were written
# for. If the directory has moved since, point them at where it is now.
first=$(find "$spectra" -maxdepth 1 -name '*.control' | head -1)
written_for=$(grep -m1 "^ FILRAW=" "$first" | sed "s/^ FILRAW='//; s|/spectra/[^/]*'\$||")
if [[ -n $written_for && $written_for != "$out_path" ]]; then
    echo "The control files were written for $written_for, rewriting them for $out_path."
    find "$spectra" -maxdepth 1 -name '*.control' -print0 \
        | xargs -0 sed -i "s|'$written_for/spectra/|'$out_path/spectra/|g"
fi
mkdir -p "$spectra/CoordFiles"

echo "Fitting $n voxels with LCModel on $cores cores."
if command -v parallel >/dev/null; then
    find "$spectra" -maxdepth 1 -name '*.control' -print0 \
        | parallel -0 -j "$cores" "$lcmodel < {} 2>/dev/null"
else
    find "$spectra" -maxdepth 1 -name '*.control' -print0 \
        | xargs -0 -P "$cores" -I{} sh -c "$lcmodel < '{}' 2>/dev/null"
fi

tables=$(find "$spectra" -maxdepth 1 -name '*.table' | wc -l)
echo "Done: $tables tables for $n voxels."
[[ $tables -gt 0 ]]
