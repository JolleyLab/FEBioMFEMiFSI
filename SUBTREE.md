# Working with Git Subtree in FEBioMFEMiFSI

This repository uses **Git subtree** to include the iFSI code. Subtree keeps the code integrated in the main repo while allowing you to push and pull changes to the original iFSI repository.

---

## Cloning the Repository
No special steps are needed for subtree. Just clone as usual:
```bash
git clone git@github.com:JolleyLab/FEBioMFEMiFSI.git
```
All files, including iFSI, will be present.

---

## Pulling Updates from iFSI
To update the iFSI subtree from its remote repository:
```bash
git subtree pull --prefix=iFSI git@github.com:ryanblack129/iFSI.git main --squash
```
**Explanation:**
- `--prefix=iFSI` → folder where subtree lives
- `main` → branch to pull from
- `--squash` → optional, reduces history to one commit

---

## Pushing Changes Back to iFSI
If you make changes in `iFSI/` and want to push them back to the iFSI repo:
```bash
git subtree push --prefix=iFSI git@github.com:ryanblack129/iFSI.git main
```

---

## Workflow Summary
- Edit files in `iFSI/` directly.
- Commit changes in the main repo as usual.
- Use `git subtree push` to sync changes back to iFSI.
- Use `git subtree pull` to bring in updates from iFSI.

---

## Best Practices
- Use `--squash` when pulling to keep history clean.
- Document any subtree updates in commit messages.