# Git Subtree Guide for MFEMFSI folder (iFSI Upstream)

This project includes the external **iFSI** codebase using git subtree, stored locally in the directory **MFEMFSI/**

Upstream repository details:

- **Repository:** https://github.com/ryanblack129/iFSI  
- **Local remote name:** `iFSI`  
- **Default branch:** `main`

This is a developer-only workflow for adding, updating, and pushing the subtree.

---

## 1. Initial Setup (run once)

If you do not yet have the `iFSI` remote configured:

```bash
git remote add iFSI https://github.com/ryanblack129/iFSI.git
git fetch iFSI
```

You can verify with:

```bash
git remote -v
```

---

## 2. Adding the Subtree

If the directory `MFEMFSI/` does **not** yet exist:

```bash
git subtree add --prefix=MFEMFSI iFSI main --squash
```

This pulls the latest `iFSI` code into your project.

---

## 3. Updating the Subtree (pulling upstream changes)

Commit or stash all working-tree changes before updating the subtree. To update
your local subtree with the latest code from the `iFSI` repository, run:

```bash
git subtree pull --prefix=MFEMFSI iFSI main --squash
```

The `--squash` flag keeps your main repo history clean by condensing subtree commits.

---

## 4. Pushing Local Changes Upstream (Developers)

After committing changes inside `MFEMFSI/`, push them to the upstream repository:

```bash
git subtree push --prefix=MFEMFSI iFSI main
```

This extracts only the subtree changes and pushes them to the `main` branch of `iFSI`.
