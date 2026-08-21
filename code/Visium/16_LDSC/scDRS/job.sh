#!/bin/bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --job-name=operon-80f5c6ba
#SBATCH --output='/users/mtotty/claude_scratch/.claude-science/jobs/80f5c6ba-5327-41e2-ba75-633ebba10c32'/slurm-%j.out
set -euo pipefail
cd '/users/mtotty/claude_scratch/.claude-science/jobs/80f5c6ba-5327-41e2-ba75-633ebba10c32'
# We are the session leader (setsid bash -lc 'exec bash job.sh' — exec replaces
# the -lc shell so this script IS the leader). Record our own pgid+starttime
# so status()/cancel() target the real process group — `$!` of the launching
# shell would be setsid's transient PID, not us.
echo "$$ $(awk '{print $22}' /proc/$$/stat 2>/dev/null || echo 0)" > _pgid
status() { printf '{"state":"%s","exitCode":%s,"ts":%s}\n' "$1" "${2:-null}" "$(date +%s)" > _status.json; }
trap 'rc=$?; [ $rc -eq 0 ] && status done $rc || status failed $rc' EXIT
# bash only processes signals during `wait`, not while a foreground command
# runs — so cancel()'s SIGTERM is invisible until cmd.sh exits unless we
# background+wait. The TERM trap forces a non-zero rc into the EXIT trap.
trap 'exit 143' TERM
status running
# Conda/Lmod activate scripts routinely reference unset vars and may have
# benign non-zero exits (e.g. from `test`); -u/-e would abort the job before
# cmd.sh runs. Relax for activation only, restore after.
set +eu

set -eu
trap 'kill -TERM ${_cmd:-} 2>/dev/null || true; exit 143' TERM
if command -v timeout >/dev/null 2>&1; then
  timeout 2700 bash -eo pipefail ./cmd.sh &
else
  bash -eo pipefail ./cmd.sh &
fi
_cmd=$!
wait $_cmd
