#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import cluster
#' @import wesanderson
#' @import shinythemes
#' @noRd

app_ui = function(request){
  tagList(
    golem_add_external_resources(),
    navbarPage(title = strong("Data Report"),
               theme = shinytheme("cerulean"),
               tabPanel("Summary",
                        icon = icon("clipboard"),
                        mainPanel(width =12,strong("General description"),
                                  htmlOutput("general")),
                        sidebarPanel(
                          fileInput("GEXP_file","Gene Expression dataset", accept=".csv", width = NULL, placeholder = "Gexp.csv"),
                          fileInput("CLINICAL_file","Clinical dataset", accept=".csv", width = NULL, placeholder = "Clinical.csv")),
                        mainPanel(
                          plotOutput("heatmap")),
                        tags$footer("sebastien.renaut@gmail.com --- 2022")),

               navbarMenu("Clustering",
                          icon = icon("circle-nodes"),
                          tabPanel("K-means",
                                   mainPanel(width =12,h4("Gene expression clustering (K-means)"),htmlOutput("kmeans")),
                                   sidebarPanel("",
                                                radioButtons("clustering_metrics", "Select clustering metric:",
                                                             c("Silhouette" = "silhouette",
                                                               "Detailed Silhouette" = "detailed",
                                                               "Total Within Sum of Square" = "twss")),
                                                br(),
                                                sliderInput("k_selected","Select a K for kmeans clustering",min = 2, max = 8, value = 2)
                                   ),
                                   mainPanel("",
                                             plotOutput('clustering_plots')
                                   ),
                                   sidebarPanel(
                                     checkboxGroupInput("pc", "Choose two PCs for visualisation",
                                                        choiceNames = paste0("PC",1:8),
                                                        selected = c(1,2),
                                                        choiceValues =  1:8)
                                   ),
                                   mainPanel(plotOutput("pca")),
                                   tags$footer("sebastien.renaut@gmail.com --- 2022")
                          ),

                          tabPanel("Hierarchical clustering",
                                   mainPanel(width =12,
                                             h4("Gene expression clustering (Hierarchical)"),
                                             htmlOutput("dendrogram")
                                   ),
                                   sidebarPanel("",
                                                radioButtons("hc_metrics", "Select clustering metric:",
                                                             c("Silhouette" = "silhouette_hc",
                                                               "Detailed Silhouette" = "detailed_hc",
                                                               "Total Within Sum of Square" = "twss_hc")),
                                                br(),
                                                sliderInput("k_selected_hc","Select a K for Hierarchical clustering",min = 2, max = 8, value = 2)
                                   ),
                                   mainPanel(plotOutput("hc_plots"),
                                             plotOutput("pca_hc")),
                                   tags$footer("sebastien.renaut@gmail.com --- 2022")
                          )
               ),



               tabPanel("Clinical",
                        icon = icon("microscope"),
                        sidebarPanel(
                          radioButtons("variable", "Choose a clinical variable to display:",
                                       selected = "age",
                                       choices= c("kmeans","hierarchical","age","stage","histology","gender","survival"))


                        ),
                        mainPanel(
                          h4("The Clinical dataset"),
                          htmlOutput("clinical"),
                          plotOutput("summary_data")),
                        mainPanel(
                          dataTableOutput("mesotable"),
                          tags$table(tableFooter("sebastien.renaut@gmail.com --- 2022")))
               ),
               navbarMenu("Downloads",icon = icon("floppy-disk"),
                          tabPanel("Clinical/Clustering",
                                   sidebarPanel(
                                     downloadButton('downloadData', 'Clinical/Clustering'))),
                          tabPanel("Gene Expression",
                                   sidebarPanel(
                                     downloadButton('downloadGexp', 'Gene Expression'))),
                          tabPanel("Plots",
                                   sidebarPanel(
                                     downloadButton('downloadHeatmap', 'Heatmap'))),






               )



    )
  )


}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "GexpCluster"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
