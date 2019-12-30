library(shiny)

shinyServer(function(input, output,session){
  
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    #BiocManager::install("msa")
    #BiocManager::install("Biobase")
    #BiocManager::install("DECIPHER")
    
    library(Biostrings)
    library(XVector)
    library(msa)
    library(shinyjs)
    library("DECIPHER")
  
    library(seqinr)
    require(ade4)
    library(Biostrings)
    
    #v <- reactiveValues(dic = NULL)
    #reactive value for example :: input$num
    #reactive function for exxample :: renderPlot({})
    #reactive value run just in reactive function
    
    #Builds a reactive object
    #
    #data <- reactive({
    #  rnorm(input$num)
    #})
    #or
  
  inFile <- "C:/Users/dell/Desktop/R__/Test_print_pdf/text_input.fasta"
  file.remove(inFile, overwrite = TRUE)
  unlink("text_input.fasta", recursive=TRUE, force = TRUE)
  
    observeEvent(input$go, {
      choice_Algorithme <- input$var
      
      inFile <- "C:/Users/dell/Desktop/R__/Test_print_pdf/text_input.fasta"
      file.remove(inFile, overwrite = TRUE)
      unlink("text_input.fasta", recursive=TRUE, force = TRUE)
      
      file_1 <- input$input_1
      file_1 <- file_1$datapath
      text_input<- input$input_2
      print("***")
      print(input$input_1$datapath)
      print(input$input_1$name)
      print("***")
      
      file_2 <- NULL
      if(text_input != ""){
      file_2<-file("text_input.fasta")
      writeLines(text_input, file_2)
      close(file_2)
      file_2 <- "C:/Users/dell/Desktop/R__/Test_print_pdf/text_input.fasta"
      }
      
      print(input$file_1)
      print(input$file_2)
      
      if(is.null(file_1) & is.null(file_2)){
        file <- NULL
      }else if(!is.null(file_1) & !is.null(file_2)){
        file <- file_1
      }else if(is.null(file_1)){
        file <- file_2
      }else if(is.null(file_2)){
        file <- file_1
      }
      
      
      if(!is.null(file)){
      shinyjs::disable("go")
      shinyjs::show("text1")
      withProgress(message = 'Traitement Alignment', value = 0, {
        n <- 10
        i <- 1
        incProgress(1/n, detail = paste("Doing part", i))
      
      #Path_upload <- isolate(input$input_1)
      #Chargement du fichiers / READ Sequences
      #mySeqs <-readDNAStringSet("C:/Users/dell/Desktop/R__/DNA.fasta")
        
      mySeqs <-readDNAStringSet(file)
      
      if(choice_Algorithme == "Default"){
          #Allignement par default
          myAlignment <- msa(mySeqs)
          print(myAlignment)
          
          
      }else if(choice_Algorithme == "ClustalW"){
          #Allignement with ClustalW
          myAlignment<-msaClustalW(mySeqs, gapOpening=1, gapExtension=1, maxiters=16,
                                   cluster="upgma", kimura=FALSE, order="input", maxdiv=23)
      }else if(choice_Algorithme == "Muscle"){
          #Allignement with Muscle
          myAlignment<-msaMuscle(mySeqs)
          print(myAlignment)
      }
      
      
      i <- 2
      incProgress(1/n, detail = paste("Doing part", i))
      #Masquer les seqs
      rowM <- IRanges(start=1, end=2)
      myMaskedAlignment <- myAlignment
      rowmask(myMaskedAlignment) <- rowM
      myMaskedAlignment
      unmasked(myMaskedAlignment)
      
      i <- 3
      incProgress(1/n, detail = paste("Doing part", i))
      ## show resulting LaTeX code with default settings
      msaPrettyPrint(myAlignment, output="asis", askForOverwrite=FALSE)
      
      i <- 4
      incProgress(1/n, detail = paste("Doing part", i))
      ## create PDF file according to some custom settings
      tmpFile <- tempfile(pattern="msa", tmpdir=".", fi=".pdf")
      tmpFile
      myAlignment_Identical <- myAlignment
      myAlignment_Similar <- myAlignment
      
      i <- 5
      incProgress(1/n, detail = paste("Doing part", i))
      msaPrettyPrint(myAlignment_Identical, output="tex",
                     showNames="left", showNumbering="none", showLogo="top",
                     showConsensus="bottom", logoColors="rasmol",shadingMode =c("identical"),
                     verbose=FALSE, askForOverwrite=FALSE)
      
      i <- 6
      incProgress(1/n, detail = paste("Doing part", i))
      tools::texi2pdf("myAlignment_Identical.tex")
      
      i <- 7
      incProgress(1/n, detail = paste("Doing part", i))
      msaPrettyPrint(myAlignment_Similar, output="tex",
                           showNames="left", showNumbering="none", showLogo="top",
                           showConsensus="bottom", logoColors="rasmol",shadingMode =c("similar"),
                           shadingModeArg = 20,verbose=FALSE, askForOverwrite=FALSE)
      
      i <- 8
      incProgress(1/n, detail = paste("Doing part", i))
      tools::texi2pdf("myAlignment_Similar.tex")
      
      i <- 9
      incProgress(1/n, detail = paste("Doing part", i))
      library(filesstrings)
      
      print("finish")
      print("Start moving")
      destDir <-paste(getwd(),"www", sep="/")
      inFile <-paste(getwd(),"myAlignment_Identical.pdf", sep="/")
      file.move(inFile,destDir, overwrite = TRUE)
      inFile <-paste(getwd(),"myAlignment_Similar.pdf", sep="/")
      file.move(inFile,destDir, overwrite = TRUE)
      print("finish moving")
      
        output$pdfviewer_1 <- renderText({
          return(paste('<iframe style="height:600px; width:100%" src="myAlignment_Identical.pdf"></iframe>', sep = ""))
        })

        output$pdfviewer_2 <- renderText({
          return(paste('<iframe style="height:600px; width:100%" src="myAlignment_Similar.pdf""></iframe>', sep = ""))
        })
        
        i <- 10
        incProgress(1/n, detail = paste("Doing part", i))
        shinyjs::enable("go")
        shinyjs::hide("text1")
        
        plotInput <- reactive({
          seqs<-readDNAStringSet(file)
          seqs <- OrientNucleotides(seqs)
          aligned <- AlignSeqs(seqs)
          d<-DistanceMatrix(aligned)
          image(1:dim(d), 1:dim(d), d, axes = FALSE, xlab="", ylab="")
          axis(2, 1:dim(d), rownames(d)[1:dim(d)], cex.axis = 0.9, las=1)
          axis(1, 1:dim(d), colnames(d)[1:dim(d)], cex.axis = 0.9, las=3)
          text(expand.grid(1:dim(d), 1:dim(d)), sprintf("%0.1f", d), cex=0.8)
        })
        
        output$test__ <- renderPlot({
          print(plotInput())})
        
          
        })
    }})
    
    
    observeEvent(input$reset_button, {
                  js$reset()
                 })
    plotInput <- reactive({
      seqs<-readDNAStringSet(file)
      seqs <- OrientNucleotides(seqs)
      aligned <- AlignSeqs(seqs)
      d<-DistanceMatrix(aligned)
      image(1:dim(d), 1:dim(d), d, axes = FALSE, xlab="", ylab="")
      axis(2, 1:dim(d), rownames(d)[1:dim(d)], cex.axis = 0.9, las=1)
      axis(1, 1:dim(d), colnames(d)[1:dim(d)], cex.axis = 0.9, las=3)
      text(expand.grid(1:dim(d), 1:dim(d)), sprintf("%0.1f", d), cex=0.8)
    })
    
    
    
    output$image1 <- renderImage({
      
        return(list(
          src = "C:/Users/dell/Documents/R/Projet/Projet/Image1.png",
          contentType = "image/png"
        ))
     
      
    }, deleteFile = FALSE)
    
    
  })