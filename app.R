library(shiny)
library(shinythemes)
library(DT)
library(stringr)
library(ggplot2)
library(shinydashboard)
setwd("./")
source("Functions/ClusterAnalysis.R")
if(version$os=="linux-gnu") options(bitmapType='cairo')


load("./Data/OntologyAllGenes.RData")


listFiles <- list.files("./Clusterings/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)

load("./Clusterings/CO2NoIronStarv.RData")
clustList <- c("All",unique(cluster[[1]]))
names(clustList) = clustList
clustList <- lapply(split(clustList, names(clustList)), unname)


ui <- dashboardPage(skin="black",
                    
                    dashboardHeader(title = "Coseq clustering visualisation"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Clustering view", tabName = "view", icon = icon("project-diagram")),
                        menuItem("Ontologies", tabName = "Ontologies", icon = icon("table"))
                        #menuItem("Ontology database", tabName = "Ontology database", icon = icon("seedling"))
                      )
                    ),
                    
                    dashboardBody(
                      
                      tabItems(
                        # First tab content
                        tabItem(tabName = "view",
                                
                                fixedRow(
                                  column(
                                    width = 6,
                                    selectInput("select", label = h3("Select list of genes"), width = 700,
                                                choices = files, selected = "CO2NoIronStarv.RData"), hr(),
                                    plotOutput("profiles", height="800px", width="1150px"),
                                    checkboxInput("boxplot", "Profiles as boxplots", TRUE)
                                  ),
                                  column(
                                    width = 5,
                                    selectInput("k", label = h3("Select cluster"), width = 300,
                                            choices = clustList, selected = "All"), hr(),
                                    tabsetPanel(type = "tabs",
                                                tabPanel("Nitrate pathways enrichment", plotlyOutput("rank", height="700px", width="900px")),
                                                tabPanel("GLM coefficients", plotlyOutput("coefsPlot", height="700px", width="900px")),
                                                tabPanel("GLM summary", verbatimTextOutput("summary"))
                                  )
                                )
                                )
                        ),
                        tabItem("Ontologies", 
                                DT::dataTableOutput("Ontologies", width = 1500))
                                
                      )
                    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
   cluster <- reactive({
     load(paste0("./Clusterings/",input$select))
     clustList <- c("All",unique(cluster[[1]]))
     names(clustList) = clustList
     clustList <- lapply(split(clustList, names(clustList)), unname)
     updateSelectInput(session = session, inputId = "k", choices=clustList, selected = "All")
     cluster
   })
  
  output$profiles <- renderPlot({
    
    if(input$k == "All"){plotProfile(cluster(), boxplot = input$boxplot)}
    else{plotProfile(cluster(), input$k, boxplot=input$boxplot)}
  })
  
  output$Ontologies <- DT::renderDataTable({
    
    if(input$k == "All"){findNitrateGenes(cluster())}
    else{findNitrateGenes(cluster(), input$k)}
  })
  
  output$Ontology <- DT::renderDataTable({
    ontologies[input$gene,]
  })
  
  output$rank <- renderPlotly({
    rankClusters(cluster())
  })
  
  glm <- reactive({
    
    if(input$k != "All"){

      glmCluster(DEgenes = names(cluster()[[1]][cluster()[[1]]==input$k]), 
                 normalized.count = data.frame(cluster()[[2]]@tcounts))
    }
    else{print("No cluster selected")}

  })
  
  output$coefsPlot <- renderPlotly({
    print(glm())
    plotGlmCluster(glm())
  })
  
  output$summary <- renderPrint({
    summary(glm())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

