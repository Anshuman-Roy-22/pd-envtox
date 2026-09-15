# Transfer this branch to the existing Windows clone

Download `pd-envtox-natb-mechanism.bundle` into Downloads. It contains the complete branch update since `4b2c005`, including the intervening analyses and the NatB extension. The required starting commit is already on the verified remote branch.

In PowerShell, run each command in order and stop if any command reports an error. Review local changes before switching or merging; do not discard your research files.

```powershell
cd "$env:USERPROFILE\pd-envtox"
git status --short
git switch aan-v2-overhaul

git fetch "$env:USERPROFILE\Downloads\pd-envtox-natb-mechanism.bundle" "refs/heads/aan-v2-overhaul:refs/remotes/bundle/aan-v2-overhaul"
git merge --ff-only refs/remotes/bundle/aan-v2-overhaul
git push origin aan-v2-overhaul

git rev-parse HEAD
git ls-remote --heads origin main aan-v2-overhaul
```

The remote `aan-v2-overhaul` tip should match your local HEAD and the tip reported by `git bundle list-heads` on the downloaded bundle. Remote `main` should remain `04df4ed7075c1056b5893777e23e16cdc5162e62`. These commands do not push to `main`.

The standalone analysis command is documented in `results/v2/natb_selectivity/README.md`. This bundle does not contain a completed whole-project R clean-clone verification or an AAN submission release tag.
