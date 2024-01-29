Lab Notes
=========

Lab 404 Laboratory Notebook.
Department of Molecular Biophysics, Bogomoletz Institute of Physiology, NAS of Ukraine.


## Diary
Daily log files detailing laboratory activities.

## Plasmids
Notes on plasmid design and cloning procedures.

## Projects
Directories for various projects.
Each project directory should have a unique tag and year in its name, including a head file (__rational.md__), an experiment plan file with a detailed description of the research design, and an individual folder for experiment results.

## Protocols
Documentation of standard protocols and a __Solutions Cookbook__ file with compositions of solutions.

---

#### GitHub cheat-sheet

**git** commands that you may needed for work are present below, enjoy. But first you should create fork of Labnote repo on your GitHub account.


- Clone `master` from GitHub repo: `git clone https://github.com/wisstock/Labnote.git`
- Delete my notes
- Create new `git remote` for you cloned repo:
```
git remote add personal_remote https://github.com/user_name/Labnote.git
```

&nbsp;

- For editing notes go to **Labnote** folder and switch to `master` branch: `git checkout master`
- For saving changes in Labnote folder follow next steps:
```
git add . 
git commit -m "commit message"
git push origin master
```