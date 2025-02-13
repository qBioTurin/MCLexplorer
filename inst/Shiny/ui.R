library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readxl)
library(sas7bdat)
library(survival)
library(survminer)
library(shinybusy)
library(shinyalert)
library(connector)
library(shinydashboard)
library(ggmosaic)
library(purrr)
library(car)
library(cowplot)
library(DT)

qual_vars<-c("Progression-free survival (PFS)" = "PFS","Overall Survival (OS)" = "OS","time to progression (TTP)" = "TTP",
             "TP53_loss_or_mut",
             "TP53" = "TP53", "TP53_delation" = "TP53_loss_array",
             "TP53 disruption" = "TP53_disruption",
             "Ki67_Classes (0 = .. ; 1 = ...; 2 = ....)" = "Ki67_Classes",
             "Morphologic BM infiltration" = "BM_INFILTRATION","
             HIGH LDH" = "HIGH_LDH\n",
             "ECOG" = "PS_ECOG","HISTOLOGY" = "HISTOLOGY","STAGE" = "AA_STAGE",
             "MARKER (0=NO; 1=IGH; 2=BCL1; 3=BOTH)" = "MARKER (0=NO; 1=IGH; 2=BCL1; 3=BOTH)",
             "Ki67_Classes_1 (0 = .. ; 1 = ...)" = "Ki67_Classes_1","MIPI Classes" = "MIPI_Classes",
             "Clinical response pre ASCT" = "CLINICAL_SIT_PRE_ASCT",
             "Clinical response post ASCT" = "CLINICAL_SIT_postASCT")
quant_vars<-c("Age" = "Age","Flow PB" = "flow PB","Flow BM" = "flow BM",
             "% Morphologic infiltration in BM" = "%_BM_INFILTR\n","MIPI" = "MIPI",
             "HB at the diagnosis" = "HB_LEVEL\n",
             "Neutrophils Counts at the diagnosis" = "NEUTRO_COUNT",
             "Lymphocytes Counts at the diagnosis" = "LYMPHO_COUNT", "Platelets" = "PLTS",
             "% KI_67" = "KI_67")

vars_values = c(qual_vars,quant_vars)

# check TP53_loss_or_mut TP53_disruption sono uguali? si, togliere TP53_loss_or_mut
# TP53_disruption = TP53 + TP53_loss_array (TP53_delation)
# Ki67_Classes_1 -> fare come marker
# MIPI_C -> calcolare dalla referenza

pharmagen_muts<-c("ABCB1 1236 C>T","ABCB1 2677 G>T,A","ABCB1 2677 G>T,A _2","ABCB1 3435 C,T",
                  "Aplotype ABCB1","VEGFA -2055 A>C","VEGFA -2055 A>C _2",
                  "ABCG2 421 C>A","FCGR2A 497 A>G","NCF4 -368 G>A","GSTP1 313 A>G",
                  "CRBN_1_rs1714327","CRBN_2_rs1705814") 


ui <- dashboardPage(
  dashboardHeader(title = "MCL explorer",
                  tags$li(
                    class = "dropdown d-flex align-items-center",
                    tags$head(tags$link(rel = "shortcut icon", href = "MCLlogoSMALL.png")),
                    tags$style(".main-header {max-height: 60px;}
                 .icon-container {
        position: relative;
        display: inline-block;
      }
      .icon-container .icon-text {
        visibility: hidden;
        width:400px;
        background-color: #333;
        color: #fff;
        text-align: center;
        border-radius: 6px;
        padding: 5px;
        position: absolute;
        z-index: 2;
        top: 50%;
        left: 110%;
        transform: translateY(-50%);
        opacity: 0;
        transition: opacity 0.3s;
      }
      .icon-container:hover .icon-text {
        visibility: visible;
        opacity: 1;
      }
      h3 {
        padding-top: 0px; /* Adjust the top padding */
        padding-bottom: 0px; /* Adjust the bottom padding */
        margin-top: 0; /* Adjust the top margin */
        margin-bottom: 0; /* Adjust the bottom margin */
      }"
                    )
                  ),
                  tags$li(a(onclick = "window.open('https://github.com/qBioTurin/MCLexplorer')",
                            href = NULL,
                            icon("github"),
                            title = "GitHub",
                            style = "cursor: pointer;"),
                          class = "dropdown")
  ),
  dashboardSidebar(
    sidebarMenu(id ="SideTabs",
                menuItem("Home", tabName = "Home", icon = icon("home")),
                menuItem(HTML('<img src="cluster.png" height="20" width="20"> Clustering Exploration'),
                         tabName = "mcl208"),
                menuItem( HTML('<img src="classification.png" height="20" width="20"> Classification Exploration'),
                          tabName = "youngerMCL_classification"),
                menuItem(  tabName = "Userclassification",
                           HTML('<img src="profile.png" height="20" width="20"> Costum Classification'), 
                           menuSubItem("Upload", tabName = "user_upload"),
                           menuSubItem("Classification", tabName = "user_classification") )
    )
  ),
  dashboardBody(
    tabItems(
      
      tabItem(tabName = "Home",
              fluidRow(
                column(12, align = "center",
                       h1("Welcome to the MCL MRD Connector Cluster Analysis App", style = "color:#1d8fbd; font-weight:bold;"),
                       br(),
                       
                       p(img(src = "MCLlogo.png", height = "30%", width = "30%"), align = "center"),
                       # Brief description
                       p("This application is designed to aid researchers in inspecting and classifying CONNECTOR clusters for Minimal Residual Disease in Mantle Cell Lymphoma.",
                         style = "font-size:18px; color: #333333;"),
                       br(),
                       
                       # Application features
                       div(
                         style = "font-size:16px; color:#666666; padding:10px; text-align:left;",
                         tags$ul(
                           tags$li(strong("Clustering Exploration:"), " Visualize and explore the MCL0208 dataset, which contains key information about connector clusters."),
                           tags$li(strong("Clasification Exploration:"), " Classify data from the Young dataset based on clusters identified in the MCL0208 analysis."),
                           tags$li(strong("Costum Classification:"), " Upload and classify your own dataset to apply insights gained from the MCL0208 clusters to new data.")
                         )
                       )
                       # 
                       # # Visual appeal
                       # div(
                       #   style = "font-size:14px; color: #999999; text-align:center; padding: 10px;",
                       #   "This application is designed to streamline the analysis of connector clusters, providing insights that help in the identification and tracking of MCL MRD disease progress and response to treatments."
                       # ),
                       # column(width = 6,
                       #        h2(em("Inspection of MRD kinetics in Mantle Cell Lymphoma")),
                       #        h4(" Fondazione Italiana Linfomi (FIL) MCL0208 Clinical Trial")
                       # )
                ),
                p(img(src = "Logo_QBio.png", height = "15%", width = "15%", style = "margin:20px 0px"), align = "center"),
                p(img(src = "loghi.png", height = "60%", width = "100%", style = "margin:20px 0px"), align = "center"),
                column(12, align = "center",
                       hr(),
                       div(
                         style = "font-size:14px; color: #666666;",
                         HTML("For more information on the connector methodology, please refer to:"),
                         HTML("Simone Pernice, et al. 'CONNECTOR, fitting and clustering of longitudinal data to reveal a new risk stratification system', Bioinformatics (2023).")
                       ),
                       br(),
                       tags$a(
                         href = "https://qbioturin.github.io/connector/",
                         target = "_blank",
                         icon("github"), " GitHub Repository",
                         style = "color: #1d8fbd; font-size:14px; text-decoration:none;"
                       )
                )
              )
      ),
      tabItem(tabName = "mcl208",
              h2("FIL-MCL208 Exploration"),
              fluidRow(
                column(6,
                       selectInput(inputId = "tissueMCL",
                                   label = "Tissue:",
                                   choices = c("Bone Marrow" = "BM","Peripheral Blood" = "PB"),
                                   selected = "Bone Marrow")
                )
                # column(6,
                #        selectInput(inputId = "ArmMCL",
                #                    label = "Arm:",
                #                    choices = c("No Division" = "NULL","Observational (Arm 0)" = "0","Experimental (Arm 1)" = "1"),
                #                    selected = "No Division")
                # )
              ),
              fluidRow(
                box(width = 12,
                    title = "Clustering Analysis", status = "primary", solidHeader = TRUE,
                    column(12,
                           tabsetPanel(id = "panelsMCL",
                                       tabPanel("Clustering", value = "panel_Plot",
                                                fluidRow(
                                                  column(12, plotOutput("MCL_clusteringPlot",height = "800px"))
                                                )
                                       ),
                                       tabPanel("Survival analysis",value = "panel_survMCL",
                                                fluidRow(
                                                  column(6, selectInput("selectSurvMCL",selected = "Clusters",
                                                                        label = "Survival analysis grouping by:",
                                                                        choices = c("Clusters","Merging Clusters",
                                                                                    "Arm - all clusters","Arm - single cluster") ))
                                                  # conditionalPanel(condition = "input.ArmMCL == 'NULL'",
                                                  #                  column(6, selectInput("selectSurvMCL",selected = "Clusters",
                                                  #                                        label = "Survival analysis grouping by:", choices = c("Clusters","Merging Clusters","Arm - all clusters","Arm - single cluster") ))
                                                  # ) ,
                                                  # conditionalPanel(condition = "input.ArmMCL != 'NULL'",
                                                  #                  column(6, selectInput("selectSurvMCL",selected = "Clusters",
                                                  #                                        label = "Survival analysis grouping by:", choices = c("Clusters","Merging Clusters") ))
                                                  # ) 
                                                ),
                                                fluidRow(
                                                  column(12, uiOutput("MCL_survUIPlot"))
                                                )
                                       ),
                                       tabPanel("Entropy",value = "panel_entropy",
                                                fluidRow(
                                                  column(6,
                                                         sliderInput(inputId = "entropy",
                                                                     label = "Entropy:",
                                                                     min = 0,
                                                                     max = 0,
                                                                     value = c(0,0)  )
                                                  ),
                                                  column(6,
                                                         sliderInput(inputId = "length",
                                                                     label = "Curve Length:",
                                                                     min = 0,
                                                                     max = 0,
                                                                     value = c(0,0) )
                                                  )
                                                ),
                                                fluidRow(
                                                  column(6,
                                                         plotOutput("MCL_linePlot")
                                                  ),
                                                  column(6,
                                                         plotOutput("MCL_scatterPlot")
                                                  )
                                                ),
                                                fluidRow(
                                                  column(10,
                                                         selectInput("SplineID",
                                                                     label = "Select id curve to check the fitting",
                                                                     choices = ""
                                                         ),
                                                         plotOutput("MCL_FittedCurve")
                                                  )
                                                )
                                       ),
                                       tabPanel("Clinical assestment", value = "panel_Clinical",
                                                fluidRow(
                                                  column(3,
                                                         checkboxInput("ClusterCheckMCL",label = "Unify the CONNECTOR clusters.")
                                                  )
                                                ),
                                                box(title = h2("Statistical Analysis"),width = 12,collapsible = T,collapsed = T,
                                                    status = "primary",solidHeader = T,
                                                    fluidRow(
                                                      column(6,
                                                             # conditionalPanel(condition = "input.ClusterCheckMCL",
                                                             #                  selectInput(inputId = "qual_varsMCL",
                                                             #                              label = "Variables:",
                                                             #                              choices = c("MIPI",vars_values) )
                                                             # ),
                                                             conditionalPanel(condition = "!input.ClusterCheckMCL",
                                                                              selectInput(inputId = "varsMCL",
                                                                                          label = "Variables:",
                                                                                          choices = vars_values)
                                                             )
                                                      )
                                                    ),
                                                    fluidRow(
                                                      column(12, plotOutput("MCL_ClinicalPlot",height = "400px"))
                                                    )
                                                ),
                                                box(title = h2("Qualitative Analysis - Gene Mutations"),width = 12,
                                                    collapsible = T,collapsed = T,
                                                    status = "primary",#solidHeader = T,
                                                    fluidRow(
                                                      # column(3,
                                                      #        selectInput(inputId = "cond_typeMCL",
                                                      #                    label = "Types:",
                                                      #                    choices = c("CHIP panel - anyTimeAndTissue","CHIP panel - DIAGNOSIS","CHIP panel - anyTimeAndTissue","CHIP panel - lost","CHIP panel - gained","MCL panel"),
                                                      #                    selected = "All_anyTimeAndTissue")
                                                      # ),
                                                      column(3,
                                                             selectInput(inputId = "var_typeMCL",
                                                                         label = "Variables:",choices = "")
                                                      )
                                                    ),
                                                    fluidRow(
                                                      column(12, plotOutput("MCL_MutationPlot",height = "400px"))
                                                    )
                                                ),
                                                box(title = h2("Pharmacogenomic Analysis"),width = 12,
                                                    collapsible = T,collapsed = T,
                                                    status = "primary",solidHeader = T,
                                                    fluidRow(
                                                      selectInput(inputId = "pharma_varsMCL",
                                                                  label = "Variables:",
                                                                  choices = pharmagen_muts,
                                                                  selected = "ABCB1 1236 C>T")
                                                    ),
                                                    fluidRow(
                                                      column(12, plotOutput("MCL_pharmaPlot",height = "400px"))
                                                    )
                                                )
                                       ),
                                       tabPanel("Table", value = "panel_table",
                                                fluidRow(
                                                  DT::DTOutput("DTtableMCL")
                                                  
                                                )
                                       )
                           )
                    )
                )
              )
      ),
      tabItem(tabName = "youngerMCL_classification",
              h2("YoungerMCL Exploration"),
              fluidRow(
                column(width = 6,
                       selectInput(inputId = "tissueBox",
                                   label = "Tissue:",
                                   choices = c("Bone Marrow" = "BM","Peripheral Blood" = "PB"),
                                   selected = "Bone Marrow")
                )
              ),
              box(width = 12,
                  title = "Classification Analysis", status = "primary", solidHeader = TRUE,
                  column(4,
                         fluidRow(
                           column(width = 6,
                                  conditionalPanel(condition="input.YoungPanels == 'panel_YoungPlot'", 
                                                   sliderInput(inputId = "Cut",
                                                               label = "Truncate the curves at",
                                                               min = 0, max = 1, value = 0),
                                                   actionButton(inputId = "GoClass", label = "Start Classification\n and Landmark "),
                                                   h2("")
                                  ),
                                  selectInput(inputId = "youngMCLselectColor",
                                              label = "Color by:",
                                              choices = c("None", "TTP","Random","INIERG"),
                                              #choices = c("None", "ttpevent","rnd1","INIERG"), 
                                              selected = "None" )
                                  #selectInput(inputId = "youngMCLselectColor", label = "Color by:", choices = c("None") )
                           ),
                           column(width = 12,
                                  conditionalPanel(condition="input.YoungPanels == 'panel_YmclSurv'", 
                                                   verbatimTextOutput("YmclSummary")
                                  ),
                                  conditionalPanel(condition="input.YoungPanels == 'panel_YoungPlot'", 
                                                   verbatimTextOutput("truncationSummary")
                                  )
                           )
                         )
                  ),
                  column(8,
                         tabsetPanel(id = "YoungPanels",
                                     tabPanel("Classification and Survival analysis", value = "panel_YmclSurv",
                                              fluidRow(
                                                column(12, plotOutput("classYoungPlot",height = "800px")),
                                                column(12, plotOutput("survYoungPlot"))
                                              )
                                     ),
                                     tabPanel("Data processing", value = "panel_YoungPlot",
                                              fluidRow(
                                                column(12, plotOutput("linePlot")),
                                                column(12, plotOutput("lineTruncPlot"))
                                              )
                                     ),
                                     tabPanel("Classification & Landmark analysis",value = "panel_Classification",
                                              fluidRow(
                                                column(12, 
                                                       plotOutput("classifiedCurve")
                                                )
                                              ),
                                              fluidRow(
                                                column(2,offset = 9, 
                                                       actionButton(inputId="saveLP_fromClass",
                                                                    label = "Save Landmark Point")
                                                )
                                              )
                                     ),
                                     tabPanel("Saved Landmark results",value = "panel_Summary",
                                              fluidRow(
                                                column(12,
                                                       uiOutput("landmarkPoints")
                                                )
                                              )
                                     )
                         )
                  )
              )
      ),
      tabItem(tabName = "user_upload",
              h2("Classification"),
              # Add your UI components for Classification here
              fluidRow(
                box(title = "Input Excel Data", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    fileInput("file1", accept = c(".xlsx"),
                              div(class = "icon-container",
                                  h4("Choose Excel File:", icon("info-circle")),
                                  div(class = "icon-text", "The first file should be an Excel file with three columns:\n ID, Time, and Observation. The time unit must be in days.")
                              )),
                    uiOutput("column_selectors")
                ),
                box(title = "Input TXT Data", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    fileInput("file2", accept = c(".txt"),
                              div(class = "icon-container",
                                  h4("Choose Excel File:", icon("info-circle")),
                                  div(class = "icon-text", "The second file should be a TXT file with\n an ID column matching the Excel file's ID column,\n and additional feature columns.")
                              )
                    ),
                    checkboxInput("header2", "Header", TRUE),
                    radioButtons("sep2", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = "\t"),
                    radioButtons("quote2", "Quote", choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"'),
                    uiOutput("txt_id_selector")
                ),
                box(title = "Data Preview", status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(6, tableOutput("contents1")),
                      column(6, tableOutput("contents2"))
                    ),
                    fluidRow(
                      column(6, verbatimTextOutput("info1")),
                      column(6, verbatimTextOutput("info2"))
                    )
                ),
                box(title = "Line Plot", status = "primary", solidHeader = TRUE,
                    selectInput("color_col", "Select Color Column:", choices = "None"),
                    plotOutput("linePlotUser")
                ),
                fluidRow(
                  column(3,
                         actionButton(inputId = "GoClass_Tabuser", label = "Go to the Classification")       
                  )
                )
              )
      ),
      tabItem(tabName = "user_classification",
              h2("User Data Exploration"),
              fluidRow(
                column(width = 4,
                       selectInput(inputId = "tissueBox_user", 
                                   label = "Tissue:",
                                   choices = c("Bone Marrow" = "BM","Peripheral Blood" = "PB"),
                                   selected = "Bone Marrow")
                ),
                column(width = 4,
                       selectInput(inputId = "featureTimeSurv_user", 
                                   label = div(class = "icon-container",
                                               h5(tags$b("Feature to use for the time: "), icon("info-circle")),
                                               div(class = "icon-text", "The unit time of the event should be in days")
                                   ), choices = "" )
                ),
                column(width = 4,
                       selectInput(inputId = "featureEventSurv_user", label = "Feature to use for the event:", choices = "" )
                )
              ),
              box(width = 12,
                  title = "Classification Analysis", status = "primary", solidHeader = TRUE,
                  column(4,
                         fluidRow(
                           column(width = 6,
                                  conditionalPanel(condition="input.UserPanels == 'panel_UserPlot'", 
                                                   sliderInput(inputId = "Cut_user",
                                                               label = "Truncate the curves at",
                                                               min = 0, max = 1, value = 0),
                                                   actionButton(inputId = "GoClass_user", label = "Start Classification\n and Landmark "),
                                                   h2("")
                                  ),
                                  selectInput(inputId = "selectColor_user", label = "Color by:", choices = c("None") )
                           ),
                           column(width = 12,
                                  conditionalPanel(condition="input.UserPanels == 'panel_UserSurv'", 
                                                   verbatimTextOutput("UserSummary")
                                  ),
                                  conditionalPanel(condition="input.UserPanels == 'panel_UserPlot'", 
                                                   verbatimTextOutput("UserTruncationSummary")
                                  )
                           )
                         )
                  ),
                  column(8,
                         tabsetPanel(id = "UserPanels",
                                     tabPanel("Classification and Survival analysis", value = "panel_UserSurv",
                                              fluidRow(
                                                column(12, plotOutput("classUserPlot",height = "800px")),
                                                column(12, plotOutput("survUserPlot"))
                                              )
                                     ),
                                     tabPanel("Data processing", value = "panel_UserPlot",
                                              fluidRow(
                                                column(12, plotOutput("linePlot_user")),
                                                column(12, plotOutput("lineTruncPlot_user"))
                                              )
                                     ),
                                     tabPanel("Classification & Landmark analysis",value = "panel_UserClassification",
                                              fluidRow(
                                                column(12, 
                                                       plotOutput("classifiedCurve_user")
                                                )
                                              ),
                                              fluidRow(
                                                column(2,offset = 9, 
                                                       actionButton(inputId="saveLP_fromClass_user",
                                                                    label = "Save Landmark Point")
                                                )
                                              )
                                     ),
                                     tabPanel("Saved Landmark results",value = "panel_UserSummary",
                                              fluidRow(
                                                column(12,
                                                       uiOutput("landmarkPoints_user")
                                                )
                                              )
                                     )
                         )
                  )
              )
      )
    )
  )
)