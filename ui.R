########################################################################################
## Title: ui.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for HSCs (NK Wilson and FK Hamey)
## Date: 22 November 2017
########################################################################################


library(shiny)
library(Matrix)
library(data.table)

ui <- fluidPage(

  
  
  tags$head(
    tags$style(HTML("
                    
                    h1 {
                    padding-left: 20px;
                    padding-top: 40px;
                    padding-bottom: 40px;
                    
                    background: #007acc; /* For browsers that do not support gradients */
                    background: -webkit-linear-gradient(#007acc, white); /* For Safari 5.1 to 6.0 */
                    background: -o-linear-gradient(#007acc, white); /* For Opera 11.1 to 12.0 */
                    background: -moz-linear-gradient(#007acc, white); /* For Firefox 3.6 to 15 */
                    background: linear-gradient(#007acc, white); /* Standard syntax */
                    background-color: #007acc;
                    color: white;
                    }
                    .explain {
                    
                    padding-left: 20px;
                    padding-top: 20px;
                    padding-bottom: 20px;
                    }
                    
                    .footer{ 
                    text-align: right;    
                    bottom: 0px; 
                    color: gray;
                    margin-right: 20px
                    margin-left: 20px
                    margin-top: 200px
                    margin-bottom: 20px
                    
                    
                    }
                    .downloadBut{ 
                    text-align: center; 
                    
                    
                    }
                    .downloadHeader{
                    align:center;
                    background-color: #f0f0f5;
                    margin-right: 20px;
                    margin-left: 20px;
                    margin-top: 20px;
                    margin-bottom: 20px;
                    padding-left: 20px;
                    padding-top: 10px;
                    padding-bottom: 10px;
                    
                    }
                    
                    "))
                    ),
  
  headerPanel(
    "Single-cell RNA-seq of Haematopoietic Stem and Progenitor Cells"),
  
  #titlePanel("Single-cell RNA-seq data of Haematopoietic Stem Cells"),
  # tags$div(
  #   class="explain",
  #   tags$p("The adult blood system is maintained by cells residing the bone marrow known as haematopoietic stem and progenitor cells (HSPCs). These cells can differentiate and make different mature blood cell types.To better understand the changes in gene expression as these cells differentiate towards the mature blood lineages, we profiled 44,802 individual HSPCs using a droplet-based sequencing method. Cells were collected from mouse bone marrow in two gates: Lin- c-Kit+ (LK) and Lin- Sca-1+ c-Kit+ (LSK). Gene expression patterns in the resulting transcriptional landscape can be visualised in the force-directed graph layout of the single-cell profiles below.")
  # ),
  tabsetPanel(
    tabPanel("LSK and LK compartments",
             fluidRow(
               tags$div(
                 class="explain",
                 tags$p("The adult blood system is maintained by cells residing the bone marrow known as haematopoietic stem and progenitor cells (HSPCs). These cells can differentiate and make different mature blood cell types.To better understand the changes in gene expression as these cells differentiate towards the mature blood lineages, we profiled 44,802 individual HSPCs using a droplet-based sequencing method. Cells were collected from mouse bone marrow in two gates: Lin- c-Kit+ (LK) and Lin- Sca-1+ c-Kit+ (LSK). Gene expression patterns in the resulting transcriptional landscape can be visualised in the force-directed graph layout of the single-cell profiles below."),
                 tags$br(),
                 tags$p("This dataset has been published in", tags$a(href="http://www.bloodjournal.org/content/131/21/e1",
                        "Dahlin JS*, Hamey FK*, et. al., (2018), Blood.")),
                 

                tags$p("You can find a Python notebook with the analysis and the corresponding data in the links below:"),
                 tags$a(href="html_links/notebook_dahlin_analysis_scanpy_0pt4pt2.html",
                        'Check out the Python notebook.'),
                 tags$br(),
                 tags$a(href="html_links/Notebook_input_data.tar.gz",
                        target='blank', 
                        download = 'Notebook_input_data.tar.gz',
                        'Download the input data used in the notebook.')
               
                )
               
             ),
             fluidRow(
               column(4, offset=0.5,
                      textInput("gene", label = h3("Gene of interest"),
                                placeholder="Enter an Ensembl Gene ID or Gene Name",
                                value = "Procr")
               )),#close fluidRow and column for geneID input
             fluidRow(
               column(4, offset=0.5,
                      imageOutput("LSKLK_colLSK"),
                      imageOutput("LSKLK_colLK")
                      
               ),#close column for images
               column(4,
                      tags$br(),
                      fluidRow(
                        div(style = "height:480px;align:center;",
                            plotOutput("WTLK_LSK_plot"))
                      ),##close fluidRow for geneLevels
                      fluidRow(
                        div(style = "height:200px;margin-top:20px;margin-left:120px;align:center;",
                            downloadButton('download_WTLK_LSK_plot', 
                                           'Download plot'))
                        
                      )#close fluidRow for download button geneLevels
                      
                      
                      
               ),#close middle column for expression levels
               
               column(4,
                      tags$br(),
                      tags$br(),
                      tags$br(),
                      tags$br(),
                      tags$br(),
                      
                      fluidRow(
                        div(style = "height:300;align:center;",
                            plotOutput("WTLK_LSK_boxplot"))
                      ),#close fluidRow for boxplot
                      fluidRow(
                        div(style = "height:200px;margin-top:20px;margin-left:150px;",
                            downloadButton('download_WTLK_LSK_boxplot', 
                                           'Download plot'))
                        
                      )#Close fluidRow for download button
               )#close last column for the two fluidRows boxplots
             )#close fluidRow for all but gene ID input.
    ),#close tab panel LSK LK
    
    
    
    tabPanel("LK WT and W41",
             fluidRow(
               tags$div(
                 class="explain",
                 tags$p("To investigate the effect of a signalling mutation on the haematopoietic transcriptional landscape we also profiled LK cells from W41/W41 mice, which have impared c-kit kinase activity. Clustering was performed on WT LK cells, and cells from the mutant mice projected onto the WT clusters. Gene expression can be visualised below in both the force-directed graph layouts and in violin plots of the expression in each cluster.")
               )
             ),
             fluidRow(#create fluidRow for gene input 2.
               column (4,offset=0.5,
                       textInput("gene2", label = h3("Gene of interest"),
                                 placeholder="Enter an Ensembl Gene ID or Gene Name",
                                 value = "Procr")
                       
               )#close column for gene input 2
             ),#close fluidRow gene Input 2.
             
             fluidRow(#create fluidRow for WT LK images
               column(4,#first column will be for static image with clusters
                      imageOutput("LK_WT")
                      #img(src = "images/LK_WT.tiff",width=250,height=250)
               ),#close first column with image coloured by cluster
               column(4,# second column for gene levels WT LK
                      fluidRow(
                        div(style = "height:440px;align:center;",
                            plotOutput("WTLK_plot"))
                      ),##close fluidRow for geneLevels
                      fluidRow(
                        div(style = "height:50px;margin-top:20px;margin-left:120px;align:center;",
                            downloadButton('download_WTLK_plot', 'Download plot'))
                        
                      )#close fluidRow for download button geneLevels
                      
               ),#close second column
               
               column(4,#open last column for boxplot
                      tags$br(),
                      tags$br(),
                      fluidRow(
                        div(style = "height:300;align:center;",
                            plotOutput("WTLK_boxplot"))
                      ),#close fluidRow for boxplot
                      fluidRow(
                        div(style = "height:50px;margin-top:20px;margin-left:150px;",
                            downloadButton('download_WTLK_boxplot',
                                           'Download plot'))
                        
                      )#Close fluidRow for download button
               )
             ),
             
             
             fluidRow(#Create last row for W41
               column(4,#first column for the image
                      imageOutput("LK_W41")
                      
                      #img(src = "images/LK_W41.tiff",width=250,height=250)
               ),#close first column
               column(4,#second column for the gene levels
                      fluidRow(
                        div(style = "height:440px;align:center;",
                            plotOutput("W41LK_plot"))
                      ),##close fluidRow for geneLevels
                      fluidRow(
                        div(style = "height:50px;margin-top:20px;margin-left:120px;align:center;",
                            downloadButton('downloadPlot_W41LK_plot', 
                                           'Download plot'))
                        
                      )#close fluidRow for download button geneLevels
                      
               ),#close second column
               
               column(4,#create third column
                      tags$br(),
                      tags$br(),
                      
                      fluidRow(
                        div(style = "height:300;align:center;",
                            plotOutput("W41LK_boxplot"))
                      ),#close fluidRow for boxplot
                      fluidRow(
                        div(style = "height:50px;margin-top:20px;margin-left:150px;",
                            downloadButton('download_W41LK_boxplot', 
                                           'Download plot'))
                        
                        
                      )#Close fluidRow for download button
               )#close third column.
             )#close last fluid row W41
             
             
    ),#tabPanel Close
    tabPanel("Download data",#create tabPanel for download 
             tags$br(),
             tags$br(),
             
             fluidRow(
               column(4,
                      downloadButton('download_WT_LSK_LK_data', 'Download LSK and LK data')
               
                      
               ),
               column(4,
                      downloadButton('download_WT_LK_data', 'Download LK WT data')
                      
               ), 
               column(4,#open third column
                      downloadButton('download_W41_LK_data', 'Download W41 LK data')
                      
                      
               )#close third column
               
               
             )#close fluid Row
             
             #fluidRow(
             #  div(class="downloadHeader",
            #       
             #      tags$h4("Coordinates from plots")
                   
              # )
             #),#close first fluidrow
             
             
             #fluidRow(
            #   column(4,
             #         downloadButton('download_WT_LSK_LK_plotCoordinates', 
              #                       'Download WT LSK and LK plot coordinates')
              # ),
               #column(4,
              #        downloadButton('download_WT_LK_plotCoordinates', 
              #                       'Download WT LK plot coordinates')
              # ), 
              # column(4,#open third column
              #        downloadButton('download_W41_LK_plotCoordinates', 
              #                       'Download W41 LK plot coordinates')
              # )#close third column
              # 
               
             #),#close fluid Row
             
             
             
            # fluidRow(
            #   div(class="downloadHeader",
                   
            #       tags$h4("Clusters")
                   
            #   )
            # ),#close first fluidrow
             
             
            # fluidRow(
            #   column(4,
            #          downloadButton('download_WT_LSK_LK_clusters', 
            #                         'Download LSK and LK labels')
            #   ),
            #   column(4,
            #          downloadButton('download_WT_LK_clusters', 
            #                         'Download cell labels for WT LK clusters')
            #   ), 
            #   column(4,#open third column
            #          downloadButton('download_W41_LK_clusters', 
            #                         'Download cell labels for W41 LK clusters')
            #   )#close third column
               
               
            # )#close fluid Row
             
             
    )#close tabPanel for download 
  ),#tabset close
  tags$br(),
  tags$br(),
  tags$br(),
  tags$br(),
  
  tags$div(class="footer",
           tags$hr(),
           tags$p("Göttgens Laboratory 2017"))
                    )

