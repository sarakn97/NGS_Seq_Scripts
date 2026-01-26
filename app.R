library(shiny)
library(Seurat)
library(ggplot2)
library(bslib)
source("helpers.R")
# Load Seurat Object
setwd("/")

# User interface ----

shinyApp(
  shinyUI(
    navbarPage("My Application",
               tabPanel(title = "Feature Plot", card(plotOutput("ft_umap")), uiOutput('page1')),
               tabPanel(title = "Violin Plot", card(plotOutput("vln")), uiOutput('page2')),
               tabPanel(title = "UMAP", card(plotOutput("umap")), uiOutput('page3'))
    )
  ),
 
shinyServer(function(input, output, session) {
    output$page1 <- renderUI({
      sidebarLayout(
        sidebarPanel(
          selectInput("dataset",label = "Select Tissue.",
                      choices =
                        c("Bone Marrow",
                          "Lung" )),
          textInput("gene", label = NULL, value = "CD34"),
          selectInput("var",label = "Label Clusters?",
                      choices =
                        c("Yes",
                          "No" ),
                      selected = "Yes" ),
          selectInput("var_rep", label = "Repel Labels?",
                      choices = 
                        c("Yes",
                          "No"),
                      selected = "Yes"),
          selectInput("color1", label = "Color 1",
                      choices =
                        c("Grey",
                          "Yellow",
                          "Blue",
                          "Red"),
                      selected = "Grey"),
          selectInput(
            "color2",
            label = "Color 2",
            choices =
              c( "Grey",
                 "Yellow",
                 "Blue",
                 "Red"),
            selected = "Red"),
          selectInput(
            "split_var",
            label = "Split by Sample?",
            choices =
              c("Yes",
                "No"),
            selected = "No"),
          numericInput("size",
                       "Label Size",
                       value = 1.0),
          numericInput("point_size",
                       "Cell Size",
                       value = 0.5
          )
        ),
        mainPanel(
          plotOutput("ft_umap")
        )
      )
    })
    # PAGE 2
    output$page2 <- renderUI({
      sidebarLayout(
        sidebarPanel(
          textInput("vln_gene", label = NULL, value = "CD34",),
          selectInput("split_vln",label = "Split by Sample?",
                      choices =
                        c("Yes",
                          "No"),
                      selected = "No")
        ),
        mainPanel(
          plotOutput("vln")
        )
      )
    })
    # PAGE 3
    output$page3 <- renderUI({
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "group",
            label = "Group by?",
            choices =
              c(
                "Sample",
                "Annotations",
                "Clusters"
              ),
            selected = "Annotations"
          ),
          selectInput("split_map", label = "Split by Sample?",
                      choices =
                        c("Yes",
                          "No"
                        ),
                      selected = "No"
          ),  
          numericInput("pt_size",
                       "Cell Size",
                       value = 0.75
          )
        ),
        mainPanel(
          plotOutput("umap")
        )
      )
    })
    
    combined.all <- reactive({
      dataset <- switch(input$dataset,
                     "Bone Marrow" = readRDS("/home/sara/GITHUB/Seurat-scRNA-RShinyApp/BM_Sample_scRNA.rds"),
                     "Lung" = readRDS("/home/sara/GITHUB/Seurat-scRNA-RShinyApp/Lung_Sample_scRNA.rds"))
    })

    output$ft_umap <- renderPlot({
      answer <- switch(input$var,
                       "Yes" = TRUE,
                       "No" = FALSE)
      rep <- switch(input$var_rep,
                    "Yes" = TRUE,
                    "No" = FALSE)
      split_umap <- switch(input$split_var,
                           "Yes" = TRUE,
                           "No" = FALSE)
      featureplot(seu = combined.all(), gene = input$gene, lab = answer, size = input$size, repell = rep, point = input$point_size,
                  color = c(input$color1, input$color2), split = split_umap)
      
      
    })
    output$vln <- renderPlot({
      split_vln <- switch(input$split_vln,
                          "Yes" = TRUE,
                          "No" = FALSE)
      violin_plot(seu = combined.all(), gene = input$vln_gene, split_vln = split_vln)
      
      
    })
    output$umap <- renderPlot({
      split_map <- switch(input$split_map,
                          "Yes" = TRUE,
                          "No" = FALSE)
      group <- switch(input$group,
                      "Sample" = "sample",
                      "Annotations" = "main",
                      "Clusters" = "integrated_snn_res.0.6")
      umap(seu = combined.all(), grp =  group, pt= input$pt_size, split = split_map)
      
      
    })
  })
)

