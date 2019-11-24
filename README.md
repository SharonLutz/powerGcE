## powerGcE
Empirical power analysis for a gene by environment interaction with a binary outcome and a normally distributed environmental exposure.

## Installation
```
install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/powerGcE")
```

## Input
For n cases and controls (input: nCases, nControls), the SNP is generated from a binomial distribution with a specified MAF(input: MAF) and the environmental exposure E is generated from a normal distributions with user specifed mean and vairance (input: meanE, varE). The binary outcome Y is generated from a binomial distirubtion such that:

logit\[P(Y=1)\] = &beta;<sub>0</sub> + &beta;<sub>SNP</sub> SNP + &beta;<sub>E</sub> E + &beta;<sub>I</sub> E*SNP 

See the manpage for more detail regarding the input of the gxeRC function.

```
library(gxeRC)
?gxeRC # For details on this function
```

## Example
For 5,000 subjects, 3 SNPs X with MAF of 0.05, 0.01, and 0.005, respectively, and a normally distributed environmental factor Z, we generated the power for the SNP by environment interaction on the outcome Y for each SNP independently and for the joint interaction using the following commands.

```
library(powerGcE)
powerGcE(nCase=407,nControl=376,MAF=0.49,meanE=0,varE=0.99,beta0=-0.32,betaSNP=0.17,betaE=0.97,betaI=seq(-1,-0.75,by=0.05),nSim=1000,alpha=0.00000005,plot.output=T,plot.name="powerGcE.pdf",seed=1)
```

## Output
For this example, we get the following matrix and corresponding plot which output the type 1 error rate (row 1 of the matrix) and the power (row 2 and 3 of the matrix) to detect the SNP by environment interaction on the outcome. We can see from the plot below that we have the most power in this scenario when the SNP by environment interaction is tested for all 3 SNPs instead of for each SNP individually.

```
      lmX1  lmX2  lmX3 lmAll
[1,] 0.016 0.015 0.020 0.052
[2,] 0.095 0.030 0.017 0.154
[3,] 0.439 0.091 0.045 0.578
```
<img src="https://github.com/SharonLutz/gxeRC/blob/master/gxeRC.png" width="600">

## References
The power analysis used here was implemented in the following manuscript: <br/>

**Lutz SM**, Frederiksen B, Begum F, Cho MH, Hobbs B, McDonald ML, Parker MM, DeMeo DL, Jiang L, Eringher M, Young K, Foreman MG, Kinney GL, Make BJ, Lomas DA, Bakke P, Gulsvik A, Crapo JD, Silverman EK, Beaty TH, Hokanson JE. (2019) Common and Rare Variants Genetic Association Analysis of Cigarettes Per Day Among Ever Smokers in COPD Cases and Controls. *NTR*. 21(6):714-722.

