#!/bin/zsh
# usage (from the package root):  zsh install.sh
B=$HOME/css_project; SRC=$(cd "$(dirname "$0")" && pwd)/scripts
mkdir -p $B/scripts $B/data $B/synthseg $B/work $B/review $B/results
if ls $B/scripts/*.py >/dev/null 2>&1; then
  BK=$B/scripts/backup_$(date +%Y%m%d_%H%M); mkdir -p $BK; cp $B/scripts/*.py $B/scripts/*.sh $BK/ 2>/dev/null
  echo "old scripts backed up to $BK"
fi
cp $SRC/*.py $SRC/*.sh $B/scripts/ && chmod +x $B/scripts/*.sh $B/scripts/*.py
grep -q 'css_project/scripts' ~/.zshrc || echo 'export PATH="$HOME/css_project/scripts:$PATH"' >> ~/.zshrc
echo "installed to $B/scripts  - open a new Terminal window"
