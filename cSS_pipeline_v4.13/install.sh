#!/bin/zsh
# usage (from the package root):  zsh install.sh
B=$HOME/css_project; SRC=$(cd "$(dirname "$0")" && pwd)/scripts
mkdir -p $B/scripts $B/data $B/synthseg $B/work $B/review $B/results
if ls $B/scripts/*.py >/dev/null 2>&1; then
  BK=$B/scripts/backup_$(date +%Y%m%d_%H%M); mkdir -p $BK; cp $B/scripts/*.py $B/scripts/*.sh $BK/ 2>/dev/null
  echo "old scripts backed up to $BK"
fi
# v4.18.3: copy to a temporary name, then rename - a batch that is running keeps reading its old file
# (cp over a running zsh script changes it under the shell's feet)
for f in $SRC/*.py $SRC/*.sh $SRC/BUILD; do   # BUILD: commit id, filled in by git archive (v5.1)
  t=$B/scripts/.$(basename $f).new; cp $f $t && chmod +x $t && mv -f $t $B/scripts/$(basename $f) || { echo "install FAILED at $f"; exit 1; }
done
grep -q 'css_project/scripts' ~/.zshrc || echo 'export PATH="$HOME/css_project/scripts:$PATH"' >> ~/.zshrc
echo "installed to $B/scripts  - open a new Terminal window"
