# R Shiny app demo - display PDF in app as reference document
library(shiny)
library(shinyjs)
#install.packages("shinythemes")
library(shinythemes)
library(shinyBS)

# Simple shiny layout for demo sake
shinyUI(fluidPage(theme = shinytheme("flatly"),
  pageWithSidebar(
    headerPanel("DNA Multiple Sequence Alignment"),
    sidebarPanel(#shinythemes::themeSelector(),
                 shinyjs::useShinyjs()
                 ,
                 fluidRow(
                   column(11,align="center",
                          imageOutput("image1", height = 250)
                   )),
                 textAreaInput(inputId = "input_2",
                               label = "STEP 1 - Enter your input sequences",
                               value = "",
                               cols = 100, rows = 6,
                               
                               placeholder = "sequences in any supported format"),
                 fileInput(inputId = "input_1",
                           label = "Or, upload a file: ",
                           accept = c('.fasta')),
                 selectInput("var", "STEP 2 - Choose your Tool",
                             list("Default","ClustalW","Muscle")
                 ),
                 fluidRow(
                   column(6, align="center", offset = 3,
                          actionButton("go",label = "Submit"),
                          tags$style(type='text/css'))),
                 
                 fluidRow(
                   column(6, align="center", offset = 3,
                          shinyjs::hidden(p(id = "text1", "Processing..."),
                          tags$style(type='text/css')))
                 )
                 
                 ,p(HTML("<A HREF=\"javascript:history.go(0)\">Clear</A>"))
    ),
    mainPanel(
      verbatimTextOutput('result'),
        #verbatimTextOutput("stats"),
        tabsetPanel(
          # using iframe along with tags() within tab to display pdf with scroll, height and width could be adjusted
          tabPanel("Allignment with identity",htmlOutput('pdfviewer_1')),
          tabPanel("Allignment with similarity",htmlOutput('pdfviewer_2')),
          
          tabPanel("Distance Matrix",plotOutput("test__"))
          
      
      ))
  )
)
)