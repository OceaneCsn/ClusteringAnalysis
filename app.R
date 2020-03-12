


library(shiny)
library(shinythemes)
library(DT)
library(stringr)
library(ggplot2)
library(shinydashboard)
setwd("./")
source("Functions/ClusterAnalysis.R")
#options(bitmapType='cairo')

load("./Data/OntologyAllGenes.RData")


listFiles <- list.files("./Clusterings/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)

load("./Clusterings/CO2.RData")
clustList <- c("All",unique(cluster[[1]]))
names(clustList) = clustList
clustList <- lapply(split(clustList, names(clustList)), unname)

# Define UI for application that draws a histogram
ui <- dashboardPage(skin="black",
                    
                    dashboardHeader(title = "Coseq clustering visualisation"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Clustering view", tabName = "view", icon = icon("project-diagram"))
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
                                                choices = files, selected = "CO2.RData")
                                  ),
                                  column(
                                    width = 6,
                                    selectInput("k", label = h3("Select cluster"), width = 700,
                                            choices = clustList, selected = "All")
                                  )
                                ),
                                hr(),
                                fixedRow(
                                    plotOutput("profiles", height="600px", width="1400px"),
                                    
                                    hr(),

                                    tabsetPanel(type = "tabs",
                                                tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
                                                tabPanel("Nitrate pathways enrichment", plotlyOutput("rank",height="400px", width="1500px"))
                                    
                                  )
                                )
                        )
                        # tabItem(tabName = "Ontology database",
                        #         textInput("gene", "Ask me a gene! (AGI)", value = "AT1G08090"),
                        #         hr()
                        #         #DT::dataTableOutput("Ontology")
                        # )
                                
                      )
                    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  

  output$profiles <- renderPlot({
    load(paste0("./Clusterings/",input$select))
    if(input$k == "All"){plotProfile(cluster)}
    else{plotProfile(cluster, input$k)}
  })
  
  output$Ontologies <- DT::renderDataTable({
    load(paste0("./Clusterings/",input$select))
    if(input$k == "All"){findNitrateGenes(cluster)}
    else{findNitrateGenes(cluster, input$k)}
  })
  
  output$Ontology <- DT::renderDataTable({
    ontologies[input$gene,]
  })
  
  output$rank <- renderPlotly({
    load(paste0("./Clusterings/",input$select))
    rankClusters(cluster)
  })
  
  
  output$test<- renderText({
    print("coucou")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

