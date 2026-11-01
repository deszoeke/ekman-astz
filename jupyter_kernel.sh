#!/bin/bash
# Install the "Julia" Jupyter kernel (name julia) that runs kernel.sh, on the HPC:
#   ./jupyter_kernel.sh
# Then ./precompile.sh builds the caches for the project and plot/, and in Jupyter
# select the "Julia" kernel. Re-run only if the repository moves.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
dir=${JUPYTER_DATA_DIR:-$HOME/.local/share/jupyter}/kernels/julia
mkdir -p "$dir"
cat > "$dir/kernel.json" <<EOF
{
  "argv": ["$here/kernel.sh", "{connection_file}"],
  "display_name": "Julia",
  "language": "julia",
  "interrupt_mode": "signal",
  "env": {}
}
EOF
chmod +x "$here/kernel.sh"
echo "installed $dir/kernel.json"
