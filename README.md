# MCL explorer

<p align="center">
  <a href="https://qbio.di.unito.it/">
    <img src="./inst/Shiny/www/MCLlogo.png" alt="QBio Logo" width="200">
  </a>
</p>

## Getting started

```
devtools::install_gitlab("shinyapps/shinymcl",
                        host = "gitlab.di.unito.it",
                        dependencies=TRUE)
```

To run the application, it is necessary to install the following version of the CONNECTOR package:

```
devtools::install_github("qBioTurin/connector", ref="Classification",dependencies=TRUE)
```

## How to run 

```
ShinyMCL::ShinyMCL.run()
```


<p align="center">
  <a href="https://qbio.di.unito.it/">
    <img src="./inst/Shiny/www/Logo_QBio.png" alt="QBio Logo" width="200">
  </a>
</p>