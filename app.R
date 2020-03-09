

library(shiny)
library(shinythemes)
library(DT)
library(stringr)
library(ggplot2)
library(shinydashboard)
setwd("./")
source("Functions/ClusterAnalysis.R")
#options(bitmapType='cairo')


listFiles <- list.files("./Clusterings/", full.names = F)
names(listFiles) = listFiles
files <- lapply(split(listFiles, names(listFiles)), unname)

load(paste0("./Clusterings/",files[[1]]))
clustList <- c("All",unique(cluster[[1]]))
names(clustList) = clustList
clustList <- lapply(split(clustList, names(clustList)), unname)

# Define UI for application that draws a histogram
ui <- dashboardPage(skin="black",
                    
                    dashboardHeader(title = "Coseq clustering visualisation"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Clustering view", tabName = "view", icon = icon("project-diagram"))
                        #menuItem("D", tabName = "Expression_data", icon = icon("seedling"))
                      )
                    ),
                    
                    dashboardBody(
                      
                      tabItems(
                        # First tab content
                        tabItem(tabName = "view",
                                selectInput("select", label = h3("Select clustering data"), width = 2000,
                                            choices = files),
                                selectInput("k", label = h3("Select cluster"), width = 2000,
                                            choices = clustList, selected = "All"),
                                hr(),
                                fixedRow(
                                  column(
                                    width = 5,
                                    plotOutput("profiles", height="700px")
                                  ),
                                  column(
                                    width = 5,
                                    tabsetPanel(type = "tabs",
                                                tabPanel("Ontologies", DT::dataTableOutput("Ontologies")),
                                                tabPanel("Nitrate pathways enrichment", plotlyOutput("rank",height="700px"))
                                    )
                                  )
                                )
                        )
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
  
  output$rank <- renderPlotly({
    load(paste0("./Clusterings/",input$select))
    rankClusters(cluster)
  })
  
  
  output$test<- renderText({
    print(input$k)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

