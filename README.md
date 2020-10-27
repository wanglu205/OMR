# OMR
OMR is implemented as an open source R package for two-sampling Mendelian randomization analysis using omnigenic genetic architecture. 

OMR makes an omnigenic modeling assumption on the SNP effects on the exposure variable and uses genome-wide SNPs as instrumental variables without pre-selection. OMR directly uses summary statistics and imposes a more general modeling assumption of horizontal pleiotropy than that used in the Egger regression. OMR takes LD information into account through joint modeling of all SNPs and relies on a scalable composite likelihood framework to infer causal effect. 

# Installation
OMR is implemented as an R package, which can be installed from GitHub.

####  Install from GitHub
```
library(devtools)
install_github("wanglu205/OMR")
```

# Usage
The main functions is OMR. You can find the instructions and an example by '?OMR'.

# Example
```
data(exampledata)
ny <- exampledata$ny
nx <- exampledata$nx
num.per <- exampledata$num.per
l.j <- exampledata$l.j
Z <- exampledata$Z
res=OMR(ny,nx,num.per,l.j,Z)
```

# Results reproduced
All results from all methods used in the OMR paper can be reproduced at 

 <https://github.com/willa0205/OMRreproduce>.

Details in [here](https://github.com/willa0205/OMRreproduce)

## Our group

 <http://www.xzlab.org>.
