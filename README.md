# Bitter crab disease dynamics in eastern Bering Sea Tanner and snow crab
Data and analyses for [Bitter crab disease dynamics in eastern Bering Sea Tanner and snow crab: An underestimated and emergent stressor] (https://www.sciencedirect.com/science/article/pii/S016578362500044X). This dataset includes Tanner and snow crab sampled for hemolymph in established bitter crab disease monitoring sites on 2014 - 2017 EBS bottom trawl surveys. Infection status was determined via species-specific PCR assays. 2014 data was eliminated from all analyses due to concerns with DNA extraction from leg muscle tissue instead of hemolymph (2015-2017 at-sea protocol change). 

# Project PIs: 
Erin Fedewa, Pam Jensen, Hamish Small, Michael Litzow, Michael Malick, Leah Zacher, Chris Long, and Stan Kotwicki

# Project Objectives:
1) Evaluate testing accuracy of disease detection methods
2) Estimate bitter crab disease prevalence from 2015 – 2017 in eastern Bering Sea monitoring sites
3) Identify factors influencing Hematodinium spp. infection probability 

# METADATA:
================================================

Bibliographic
-------------

| Published     | 02/02/2023   |
| ------------- | ------------- |
| Keywords      | snow crab |
|               |   Tanner crab |
|               |  Bering Sea |
|               | bitter crab disease |
|               | Hematodinium |

 

Coverage
--------

### Temporal

| Begin    | 2014-06-01 |
| ------------- | ------|
| End   | 2017-08-31 |

 
### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |

![Rplot](https://github.com/Erin-Fedewa-NOAA/BCS_2014-2017/blob/master/figs/n_year.png)


Attributes
----------
One master dataset has been produced for further modeling via the "maturity_haul" script. "pcr_haul_master.csv" includes all crab biometric, disease diagnosis and haul level data, with data attributes listed below. 

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
| `year`     |        Year of specimen collection | numeric
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text  
|  `snow_cpue`   |   Station-level snow crab density   |   numeric, crab/nmi^2
|  `snow70under_cpue`   |   Station-level density of snow crab <70 mmm carapace width   |   numeric, crab/nmi^2
|  `snowimm_cpue`   |   Station-level density of all snow crab below size cutoffs specified in at-sea protocols   |   numeric, crab/nmi^2
|  `tanner_cpue`   |   Station-level Tanner crab density   |   numeric, crab/nmi^2
|  `tanner70under_cpue`   |   Station-level density of Tanner crab <70 mmm carapace width   |   numeric, crab/nmi^2
|  `tannerimm_cpue`   |   Station-level density of all Tanner crab below size cutoffs specified in at-sea protocols   |   numeric, crab/nmi^2
|  `spno` |   Unique specimen ID. First four numbers correspond to year of collection   |  numeric
|  `species_name` | Species name of specimen sampled   |  text
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `size`  |  Carapace width of specimen sampled | numeric, in mm
|  `chela`   |   Chela height (males only) used to determine maturity status   | numeric, in mm
 | `shell_cond`  |   Shell condition of specimen sampled. 0=premolt/molting, 1=softshell, 2=newshell, 3=oldshell, 4=very oldshell, 5=graveyard, 9=unk | numeric
 |  `egg_color` | Egg color of clutch (females only). 0=No eggs, 2=Purple, 3=Brown, 4=Orange, 5=Purple-Brown, 6=Pink  |  numeric
 |  `egg_cond` | Egg condition of clutch (females only). 0=No Eggs, 1=Uneyed eggs, 2=Eyed eggs, 3=Dead eggs, 4=Empty eggs cases, 5=Hatching eggs  |  numeric
 |  `clutch` | Size of clutch (females only). 0=Immature: no eggs, 1=Mature: no eggs, 2=Trace to 1/8 full, 3=1/4 full, 4=1/2 full, 5=3/4 full, 6=Full   |  numeric
|  `collection_comments`    |    Notes on datasheet from at-sea samplers   |   text
|  `pcr_result` | Bitter crab disease diagnosis via PCR assay for detection of Hematodinium spp. DNA. 0=uninfected, 1=infected, 3=undetermined   |  numeric
|  `index_site` | Eastern Bering Sea sites established for bitter crab disease monitoring. Sites 1-3=Tanner crab samples, Sites 4-6=Snow crab samples   |  numeric
|  `general_location` | Large marine ecosystem of sampling event. EBS=eastern Bering Sea   |  text
|  `collected_by` | Agency taking hemolymph samples. SAP = NOAA AFSC Shellfish Assessment Program   |  text
|  `dna_plate_no` | Plate number containing hemolymph sample   |  numeric
|  `dna_well_no` | Individual well plate number containing hemolymph sample   |  alpha-numeric
|  `preservative` | Method of preservation for tissue/blood sample   |  text
|  `sample_status` | Disposition of tissue/blood sample   |  text
|  `c_v_h` | Combined string with cruise, vessel and haul for collections   |  numeric 
|  `cruise` | Cruise ID for Bering Sea/GOA/Atlantic bottom trawl surveys. See RACEBASE or AKFIN for additional NOAA cruise metadata   |  numeric
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `visual_positive` | Infection status via visual diagnosis. NA=not recorded, 0=not infected, 1=infected   |  numeric
|  `maturity`  |  Maturity of specimen sampled, as determined by chela (males) or clutch morphology (females). 0=Immature, 1=Mature"  | numeric  
| `start_date`     |        Date of sampling | date, month/day/year
| `mid_latitude`       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C
 | `start_date`     |        Date of sampling | date, month/day/year

:::

Distribution
------------

  File                                              Format    
  ------------------------------------------------- -------- -------------------------------------------------
  `pcr_haul_master`   `csv`    [Download](https://github.com/Erin-Fedewa-NOAA/BCS_2014-2017/tree/master/data)
 ------------------------------------------------- -------- -------------------------------------------------

 # NOAA License
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.

 
