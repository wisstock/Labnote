Research notes
==============
*Borys Olifirov, 3.05.2021*

Взаимосвязь между паттернами кальций-зависимого встраивания HPCA и PIP2-багатыми регионами цитоплазматической мембраны. 

---

CA1 neurons culture
===================
## Transfection

Эффективность трансфекции Lipofectamine 2000 составляет доли процентов, необходимо опробовать Lipofectamine 3000 или кальциевую трансфекцию.

### Plasmids mass
**Масса плазмид в расчете на одну лунку (12-wells plate):**

|Plasmid|Mass|
|-|-|
|HPCA-TagRFP|???|
|PH-YFP|???|

## MβCD treatment
???


## Calcium imaging
???

## Optical system resolution
*Для выбора оптимального разрешения изображения и величины шага сканирования*

Латеральное разрешение системы определяется диаметром дика Эйри (Airy unit/AU), который зависит от длинной волны излучения и числовой апертуры:

<img src="pic/lateral.png" width="25%">

Аксиальное разрешение конфокальной системы определяется полной шириной на середине высоты (full width half maximum/FWHM) аксиальной проекции PSF также зависит от коэффициента преломления среды (*n*) и диаметра конфокальной апертуры (*D*):

<img src="pic/axial.png" width="70%">

Для использованной оптической системы *n* = 1.33 (вода), *NA* = 0.9.
Латеральное разрешение и аксиальное разрешение для избранных значений *D* для длин волн возбуждение и эмиссии приведены ниже.


##### Fluorescent agents (data from FPbase):
|Name|Exc.|Ems.|
|-|-|-|-|
|HPCA-TFP|456 (453) nm|488 (485) nm|
|EYFP-Mem|513 nm|527 nm|

#### HPCA-TFP
##### 458 nm (exc.)
dxy = 310 nm

|D (um)|dz (um)|
|-|-|
|500|1.553|
|250|1.262|
|100|1.168|

##### 488 nm (ems.)
dxy = 330 nm

|D (um)|dz (um)|
|-|-|
|500|1.610|
|250|1.331|
|100|1.242|


#### EYFP-Mem
##### 515 nm (exc.)
dxy = 349 nm

|D (um)|dz (um)|
|-|-|
|500|1.662|
|250|1.394|
|100|1.309|

##### 527 nm (ems.)
dxy = 357 nm

|D (um)|dz (um)|
|-|-|
|500|1.685|
|250|1.422|
|100|1.339|


#### 405 nm (uncaging)

dxy = 275.5 nm

|D (um)|dz (um)|
|-|-|
|500|1.457|
|250|1.143|
|100|1.037|


*Useful links:*
- http://www.hi.helsinki.fi/amu/AMU%20Cf_tut/Opt_Pinhole.htm
- https://www.leica-microsystems.com/science-lab/confocal-optical-section-thickness/
- https://www.leica-microsystems.com/science-lab/pinhole-effect-in-confocal-microscopes/







