Inhomogeneity in HPCA insertions pattern  and translocation amplitudes in compact cells: morphological properties vs PIP2 distribution
======================================================================================================================================
*Borys Olifirov, Mykyta Bobylyov, 1.08.2023*

Hippocalcin (HPCA) is a neuronal calcium sensor protein responsible for modulating neuronal functions within the hippocampus. Upon binding with calcium ions, HPCA undergoes conformational changes, leading to highly heterogeneous insertions in the cellular membrane. Previous studies have indicated that the pattern of HPCA insertion in neurons corresponds to specific conservative regions (_Dovgan, unpublished_). Furthermore, HPCA exhibits a strong affinity for the minor phospholipid PIP2 with a dissociation constant (Kd) of approximately 50 nM (_[O'Callaghan et al., 2005](doi:10.1042/BJ20051001)_). 

Under normal resting calcium concentrations in the cytoplasm, only a small fraction of HPCA localizes to the plasma membrane, ranging from 0% to 8% (_[Sheremet et al., 2020](DOI10.1007/s11062-020-09845-6)_; _Cherkas, unpublished_). However, in response to an influx of calcium, the insertions at certain plasma membrane sites can reach as high as 50% to 100%(_Olifirov, unpublished_; _Dovgan, unpublished_). 
It is important to note that HPCA's signaling potential targets, such as KCNQ potassium channels, voltage-gated calcium channels, and clathrin-mediated endocytosis machinery, are directly associated with PIP2 (_[Rodríguez-Menchaca et al., 2012](doi:10.3389/fphar.2012.00170)_; _[Dickson et al., 2014](doi/10.1073/pnas.1407133111)_; _[Jung & Haucke, 2007](doi:10.1111/j.1600-0854.2007.00595.x)_). However, the exact mechanisms of molecular signaling in these instances remain unclear. The significant local concentration of HPCA in response to an increase in calcium concentration suggests that potential regulatory pathways may not solely involve direct protein-protein interactions. HPCA could also exert non-direct influences through the modification of local biophysical properties of the lipid bilayer and PIP2 buffering.

Primary dystonia is a neurological movement disorder that is characterized by repeated or prolonged muscle contractions and postural abnormalities. Genetic studies of patients with primary DYT2 dystonia have identified mutations in the HPCA gene as a causative factor in the disease (_[Charlesworth et al., 2015](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00060-9)_). However, the precise mechanisms through which these mutations contribute to the development of movement disorders remain poorly understood. Prior findings from our laboratory indicate that the DYT2-associated mutant HPCA(N75K) exhibits heterogeneous insertions that overlap with HPCA(WT) insertion sites, but it demonstrates a markedly reduced affinity for calcium ions (_[Osypenko et al., 2019](https://doi.org/10.1016/j.nbd.2019.104529)_).

For close interplay between experimental data and molecular dinamic simulations


Звааючи на це важливо характеризувати не тільки процес зв'зя

Саме на характеризуванні поведінки мембранної форми та різницю у взаємодії HPCA(WT) та HPCA(N75K) із 


## Main points
__Interpaly between HPCA insertions pattern anp PIP2 distribution__


__Dose-dependend responce of WT and N75K to calcium ions concentration increase at different mebrane composition__


Если в модельном объекте HEK 293 в котором отсутствуют предполагаемые эндогенные мишени HPCA паттерн транслокации равномерный по всей площади мембраны и уменьшение количества PIP2 не приводит к изменению амплитуды транслокации и ее паттерна, значит наблюдаемая в нейронах картина встраивания обусловлена не гетерогенностью липидного состава внутренней поверхности плазматической мембраны. Открытым вопросом остается роль в формировании наблюдаемой картины встраивания взаимодействия HPCA с мишенями в составе мембраны: различными компонентами каналов в случае sAHP, AP2 в случае LTD (для AP2 показана аффинность по отношению к PIP2).
- **Неравномерности паттерна встраивания обусловлена сродством с PIP2 локализованном в компактных рафтах**
Неизвестно локализован ли пул PIP2 в HEK 293 в составе рафтов или же равномерно распределен в плазматической мембране (что возможно установить совместным измерением PH-domain + FP-Mem), однако если удаление PIP2 из мембраны в результате которого будет наблюдаться снижению амплитуды транслокации HPCA будет прямым доказательством роли PIP2 в процессе кальций-зависимого встраивания HPCA в плазматическую мембрану. В таком случае HPCA может играть двоякую роль в регуляторных процессах: потенциальные мишени заведомо локализованы в областях плазматической мембраны с высоким содержанием PIP2 и таем образом локальный состав мембраны обеспечивает колокализацию мишени и мембранной формы HPCA; кальциевый сенсор выступает промежуточным звеном обеспечивающим кальций-зависимую доставку цитоплазматических мишеней в области плазматической мембраны с высоким содержанием PIP2.


## Experiments plan
### Diffuson part

### Dose-dependence part


- **EYFP-Mem + HPCA-TagRFP**
  Характер распределения сайтов встраивания HPCA
- **EYFP-Mem + PH-CFP**
  Характер распределения PIP2
- **DrVSP-EYFP + PH-CFP**
  Кинетика утилизация мембранного PIP2 за счет активности DRVSP
- **DrVSP-EYFP + PH-CFP + HPCA-TagRFP**
  Влияние утилизация мембранного PIP2 на амплитуды транслокаций HPCA