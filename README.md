# SpeciesTreeConsistency
Simulation scripts and data from Wascher and Kubatko (2020)

Data are in the Data directory.

Files in Scripts directory:
1. **comp_lik.R** and **get_lik.R** carry out the likelihood analysis for CIS data
2. **example_simscript_cis.pl** is a perl script to run a complete simulation for CIS and SNP data. Trees are given below.
3. **example_simscript_multilocus.pl** is a perl script to run a complete simulation for multilocus data. Trees are the same as in (2) and are given below.
4. **simscript_anomalous.pl** is a perl script to analyze the anomalous species tree. It uses the "default" prior settings found in file **a01_master.ctl**.
5. **simscript_anomalous_yang.pl** is a perl script to analyze the anomalous species tree for data generated with theta=0.05. It uses the informative priors in file **a01_yang.ctl**.
6. **simscript_anomalous_yang2.pl** is a perl script to analyze the anomalous species tree for data generated with theta=0.01. It uses the informative priors in file **a01_yang2.ctl. 
7. The trees used in the simulation studies are as follows:


Symmetric case:


mysptree4tax ...

_bl1:
((Species1:0.5,Species2:0.5):0.5,(Species3:0.5,Species4:0.5):0.5);   

_bl2:
((Species1:1.0,Species2:1.0):1.0,(Species3:1.0,Species4:1.0):1.0);

_bl3:
((Species1:2.0,Species2:2.0):2.0,(Species3:2.0,Species4:2.0):2.0); 

_bl4:
((Species1:1.0,Species2:1.0):0.5,(Species3:1.0,Species4:1.0):0.5);  

_bl5:
((Species1:0.5,Species2:0.5):1.0,(Species3:0.5,Species4:0.5):1.0);  

_bl6:
((Species1:1.0,Species2:1.0):0.5,(Species3:0.5,Species4:0.5):1.0);






Asymmetric case:

_bl1:
(Species4:1.5,(Species3:1.0,(Species2:0.5,Species1:0.5):0.5):0.5);

_bl2:
(Species4:3.0,(Species3:2.0,(Species2:1.0,Species1:1.0):1.0):1.0);

_bl3:
(Species4:6.0,(Species3:4.0,(Species2:2.0,Species1:2.0):2.0):2.0);

_bl4:
(Species4:2.5,(Species3:1.5,(Species2:0.5,Species1:0.5):1.0):1.0);

_bl5:
(Species4:2.0,(Species3:1.0,(Species2:0.5,Species1:0.5):0.5):1.0);

_bl6:
(Species4:2.5,(Species3:2.0,(Species2:1.0,Species1:1.0):1.0):0.5);



Anomalous case:

mysptree4tax_a1 = 
(Species1:0.48,(Species2:0.44,(Species3:0.4,Species4:0.4):0.04):0.04);

