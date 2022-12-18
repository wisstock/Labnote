Lab Notes
=========

Lab 404 laboratory notebook, Department of Molecular biophysics, Bogomoletz Institute of Physiology.

## Structure of lab notes
This notebook using Markdown markup for text notes and CSV sheets format for storing reagents compositions and chemicals/equipment lists.

There are three main directories:
 1. **diary** 
 1. **projects**
 2. **protocols**
 3. **reagents**

## Diary
Contain daily notes files with everyday lab activity description and CSV sheets with chemicals/expended supplies lists.

## Projects
Contain individual directories for each project.
Each project folders must have unique tag and year in name and contain head file (**rational.md**), experiments plan file with detail description of reserach design and individual folder for experiments results.

## Protocols
Descriptions of standard protocols for laboratory manipulations. Presence of detailed protocol will simplify writing of future work notes (just insert reference to standard protocol).

## Reagents
Composition of solutions and reagents.
Data storing in CSV table format and have next columns:
 - Chemical - name of chemical
 - ID - manufacturer and serial number
 - MW (g/mol) - molecular weight
 - Stock Solution - concentration of srock solution (if you don't use SS just mark as "solid")
 - Concentration - final concentration
 - volume - mass or volume of stock solution for the specified solution volume (there may be several "volume" columns for different final solution volumes)

---

#### GitHub cheatsheet

Main **git** commands that you may needed for work are present below, enjoy. But first you should create fork of Labnote repo on your GitHub account.


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