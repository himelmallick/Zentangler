# Zentangler

**Zentangler** is an R package for performing causal multimodal mediation analysis, enabling researchers to untangle complex direct and indirect effects across diverse multimodal data types.

## Installation
You can install the development version of zentangler from GitHub:
```r
# install.packages("devtools")
devtools::install_github("himelmallick/zentangler")
```

## Overview
Mediation analysis helps researchers understand *how* an exposure affects an outcome by decomposing total effects into direct and indirect pathways through intermediate variables (mediators). When data span multiple modalities — such as imaging, genomics, text, or behavioral measures — standard single-mediator or single-view mediation approaches fall short.

**zentangler** provides a unified framework to:
- Estimate direct and indirect effects across multimodal mediators
- Handle high-dimensional multimodal mediators
- Support flexible multimodal fusion, including early, late, and intermediate
- Support multiple ML models for indirect effect estimation
- Visualize complex mediation pathways
