AM-form loading curve
================
*25.09.2020-???*

Построение кривой загрузки клеток HEK 293 (new) AM NP-EGTA и AM Fluo-4.

Регистрации произведены 25.09.2020, исходные изображения и протокол TC  в дирректории 25_09_2020 на NextCloud.


## Experiment design
Для оптимизации протокола загрузки поставлен калибровочный опыт, одно стекло с покрытием \~90% переведено по стандартному протоколу и затем помещено в стандартный Outside Solution 2 mM Ca2+  с 5 uM AM NP-EGTA и 3 uM AM Fluo-4.

В течении 1.5 h произведена серия снимков с возбуждением Fluo-4, и последующим анкейджингом. Абсолютная интенсивность флуорисценции отражала процесс загрузки Fluo-4, в то время как соотношение F0/deltaF отражало загрузку NP-EGTA.

Для снимка выбрана большая область с не слишком плотно расположенными клетками, каждую клетку на снимке обрабатывали отдельно, набирая статистику с одного поля зрения.


## Time Controller protocol
Цикл из 6 снимков в одной фокальной плоскости с интервалом в 20 секунд (120 s) и интервалом 120 s повторяются 20 раз (полное время регистрации 80'). 

Стимуляция 405 nm всего поля зрения непосредственно перед 3 фреймом в серии снимков, экспозиция 2 us/px.

Файл протокола *25_09_2020_loading_curve.otd*.


## Imaging setup
#### Fluorescent agents
|Name|Ex.|Em.|Count|
|-|-|-|-|
|Fluo-4|488 nm|516 nm|2 uM|
|NP-EGTA|405 nm|-|5 uM|

#### Initial parameters
**Optical system**
C.A.: 250 um
Exposure: 2 us/px
Image size: 128x128 px
Zoom: 6
Size: 0.276 um/px
Scaning speed: L 1.360 ms, F 0.188 s, S 1.157 s

**Uncaging**
Laser power (405 nm): 100%
Scaning mode: rect.
Scaning area diameter: full field of view
Exposure: 10 us/px
Cycles: 1


## Imaging
**Lasers setup**
|Laser|Power|
|-|-|
|488 nm|10 %|

**Registration setup**
|Channel|HV|Pass band|
|-|-|-|
|CHS1|700V|505-540 nm|
|TD1|250 V|-|
