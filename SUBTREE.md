# Working with Git Subtree in FEBioMFEMiFSI

This repository uses **Git subtree** to include the MFEMFSI code. Subtree keeps the code integrated in the main repo while allowing you to push and pull changes to the original MFEMFSI repository.

---

## Cloning the Repository
No special steps are needed for subtree. Just clone as usual:
```bash
git clone git@github.com:JolleyLab/FEBioMFEMiFSI.git
```
All files, including MFEMFSI, will be present.

---

## Pulling Updates from MFEMFSI
To update the MFEMFSI subtree from its remote repository:
```bash
git subtree pull --prefix=MFEMFSI git@github.com:ryanblack129/MFEMFSI.git main --squash
```
**Explanation:**
- `--prefix=MFEMFSI` → folder where subtree lives
- `main` → branch to pull from
- `--squash` → optional, reduces history to one commit

---

## Pushing Changes Back to MFEMFSI
If you make changes in `MFEMFSI/` and want to push them back to the MFEMFSI repo:
```bash
git subtree push --prefix=MFEMFSI git@github.com:ryanblack129/MFEMFSI.git main
```

---

## Workflow Summary
- Edit files in `MFEMFSI/` directly.
- Commit changes in the main repo as usual.
- Use `git subtree push` to sync changes back to MFEMFSI.
- Use `git subtree pull` to bring in updates from MFEMFSI.

---

## Best Practices
- Use `--squash` when pulling to keep history clean.
- Document any subtree updates in commit messages.