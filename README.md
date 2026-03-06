# MCL explorer

<p align="center">
  <a href="https://qbio.di.unito.it/">
    <img src="./inst/Shiny/www/MCLlogo.png" alt="QBio Logo" width="100">
  </a>
</p>


This application is designed to aid researchers in inspecting and classifying CONNECTOR clusters for Minimal Residual Disease in Mantle Cell Lymphoma.

1. **Clustering Exploration**: Visualize and explore the 'FIL-MCL0208' dataset, which contains key information about connector clusters.

2. **Classification Exploration**: Classify data from the 'MCL Younger' dataset based on clusters identified in the 'FIL-MCL0208' analysis.

3. **Custom Classification**: Upload and classify your own dataset to apply insights gained from the 'FIL-MCL0208' clusters to new data.
            
<p align="center">
    <img src="./inst/Shiny/www/Framework.png" alt="Framework" width="400">
</p>           
                         
## Getting started

To install the package directly from GitHub, use the following command in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("qBioTurin/MCLexplorer", dependencies=TRUE)
```

### Dependencies

To run the application, it is necessary to install the following version of the `connector` package:

```r
devtools::install_github("qBioTurin/connector", ref="Classification", dependencies=TRUE)
```

For more information on the `connector` methodology, please refer to:

[Simone Pernice, et al. 'CONNECTOR, fitting and clustering of longitudinal data to reveal a new risk stratification system', Bioinformatics (2023)](https://academic.oup.com/bioinformatics/article/39/5/btad201/7133735)

## How to run 

Once installed, you can launch the Shiny application with:

```r
library(MCLexplorer)
MCLexplorer::mclexplorer.run()
```


<p align="center">
  <a href="https://qbio.di.unito.it/">
    <img src="./inst/Shiny/www/Logo_QBio.png" alt="QBio Logo" width="200">
  </a>
</p>