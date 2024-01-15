rm(list =ls(all = T))

##IMPORTING LIBRARIES##

print("Loading libraries...")

library(shiny) # works with the chr22.vcf.gz in the VariantAnnotation extdata
library(dplyr)
library(DT)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(ReactomePA)
library(org.Hs.eg.db)
library(biomaRt)
library('hexbin')

#setting up datasets

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
target_chromosomes <- seqlevels(txdb)
all_snps <- SNPlocs.Hsapiens.dbSNP150.GRCh38
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
ensembl1 <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#Defining constant

# Transition (Ti)
ti <- c("A>G","G>A","C>T","T>C")
# Transversion (Tv)
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")

#Importing modules

print("Importing R files script...")

path_to_r_file <- "C:/Users/lilol/Documents/VCF_EXPL/shiny/modules/"
to_import <- c("Plots.R", "Variants_treatment.R", "Patterns.R","Annotations.R")
for (script in to_import)
{
  source(paste0(path_to_r_file, script), local = T)
}

##SERVER##

server <-  function(input, output, session) {
  options(shiny.maxRequestSize = 40 * 1024^2)
  
  ##FUNCTIONS##
  
  #import the data 
  data1<- reactive({ 
    validate(need(!is.null(input$fileInput$datapath), "please pick file"))
    raw_vcf_data <- readVcf(input$fileInput$datapath, "hg38")
    # Perform filtering based on the FILTER column
    
    filtered_vcf_data <- subset(raw_vcf_data, FILTER == "PASS")
    return (filtered_vcf_data)
  }) #end of function
  
  #get the sample names
  
  samples <- reactive({
    sampleID <- colnames(Depth())
    return (sampleID)
    
  }) #end of function
  
  ##QUALITY##
  
  #DEPTH DATA
  
  Depth <- reactive ({
    matDP <- geno(data1())$DP
    matDT <-as.data.frame(matDP)
    return(matDT)
    
  }) #end of function
  
  #Filter the vcf by sample
  
  SAMPLE_SELECTED <- reactive({ 
    if (is.null(input$inSelect)){ #if there is no sample selected
      sample <- data1()
      sample <- data1()[,1]
    }
    else {
      sample <- data1()
      sample <- data1()[,as.numeric(input$inSelect)]
    }
    return(sample)
  }) #end of function
  
  #Genotypes
  
  GT <- reactive({
    tbl <- table(geno(data1())$GT)
    tbl_dat <- as.data.frame(tbl)
    return(tbl_dat)
  }) #end of function
  
  
  # Histogram function to plot the SQ in function of the AF BEFORE Filtering
  
  HB <- reactive({
    hb <- hexbin(geno(SAMPLE_SELECTED())$SQ, geno(SAMPLE_SELECTED())$AF, xbins = 20)
    return(hb)
  })
  
  #Filter SQ value :
  
  sampleDT <- reactive({
    sample <- subset(SAMPLE_SELECTED(), geno(SAMPLE_SELECTED())$SQ >= as.numeric(input$SQValue))
    return(sample)
  })
  
  #transform into a VRange object easier to manipulate
  data2 <- reactive({
    myv <- as(sampleDT(), "VRanges")[seq_len(input$numvar)]
    myv <- as.data.frame(myv)%>% dplyr::select(-width)
    return(myv)
  })  #end of function

  #GT By samples
  
  GTsamples <- reactive ({ 
    matGT <- geno(data1())$GT
    matGT <- as.data.frame(matGT)
    return(matGT)
  })#end of function
  
  # Histogram function to plot the SQ in function of the AF AFTER Filtering SQ
  
  HB2 <- reactive({
    hb2 <- hexbin(geno(sampleDT())$SQ, geno(sampleDT())$AF, xbins = 20)
    return(hb2)
  })  
  
  ##PATTERNS##
  
  #Get Variants types 
  
  Variants_tab <- reactive({
    Variants<-Get_Variants(sampleDT())
    return(Variants)
  }) #end of function
  
  #Turn it into a datatable to plot
  
  Mutations <- reactive ({ 
    tbl <- table(Variants_tab()$mutType)
    tbl_dat <- as.data.frame(tbl)
    return(tbl_dat)
  }) #end of function
  
  #GET THE SNPS ONLY
  
  VarX <- reactive({
    varX <- Variants_tab()[Variants_tab()$mutType == "SNP", ] 
    return(varX)
  }) #end of function
  
  #Transitions and transversion datatable
  
  Transition <- reactive({
    tbl <- table(VarX()$nuSub)
    tbl_dat <- as.data.frame(tbl)
    return(tbl_dat)
  }) #end of function
  
  #Ti/tv ratio
  
  Ratio<-reactive({
    tbl <- table(VarX()$TiTv)
    tbl_dat <- as.data.frame(tbl)
    return(tbl_dat)
  }) #end of function
  
  #Turn vcf into a range object
  
  RD <- reactive ({
    rd <- rowRanges(sampleDT()) 
    return(rd)
  }) #end of function
  
  #IS THERE PATTERNS IN THE SAMPLE DATA? APOBEC AND Trinucleotide
  
  PATTERNS <- reactive({
    patterns <- GET_PATTERNS(RD()) #call function
    return(patterns)
  }) #end of function
  
  
  ###ANNOTATION###
  
  #Subset VCF by chromosome (analysis too long for whole exome)
  
  VCF_CHR <- reactive({
    vcf_chr<- sampleDT()[grepl(names(sampleDT()), pattern = input$CHRSelect)]
    return(vcf_chr)
  })
  
  #Get the sample Range
  
  FILTEREDBIS <- reactive({
    rd_filtered <- rowRanges(VCF_CHR())
    return(rd_filtered)
  })
  
  #GET SNPS database for the target chromosome
  
  MY_SNPS<- reactive ({
    my_snps <-GET_SNPS(FILTEREDBIS()) #call function
    return(my_snps)
  }) #end of function
  
  #GET SNP IDS for this chromosome
  
  SNPID <- reactive({
    snp_ID <- data.frame(posIDX = paste0(seqnames(MY_SNPS()), 
                                         ":", pos(MY_SNPS())), 
                         rsID = MY_SNPS()$RefSNP_id,
                         stringsAsFactors = FALSE)
    return(snp_ID)
  }) #end of function
  
  #GET A DATATABLE FOR SNPS
  
  MATS <- reactive ({
    mats <- GET_MATS(data = VCF_CHR(), #call function
                     snp_id = SNPID())
    return(mats)
  }) #end of function
  
  #ANNOTATE 
  
  ANNOTATIONS <- reactive({
    annotations <- GET_ANNOTATIONS(MATS())
    return(annotations)
  })
  
  #REACTOME PATHWAYS
  
  PATHWAYS <- reactive({
    pathways <- GET_PATHWAYS(ANNOTATIONS())
    return(pathways)
  })
  
    FILTERED_ANNOTATIONS <- reactive({
      filtered <- ANNOTATIONS()[!is.na(ANNOTATIONS()$hgnc_symbol), ]
      filtered$AF <- as.numeric(filtered$AF)
      return(filtered)
  })
  
  #VARIATIONS IN DBSNP DATATABLE
  
  VARIATIONS <- reactive ({
    taC2 <- table(!is.na(MATS()$rsID))
    taC2_dat <- as.data.frame(taC2)
    return(taC2_dat)
  })#end of function
  
  #GET THE CODING VARIANTS
  
  MATA <- reactive({
    mata <- GET_MATA(sampleDT()) #call function
    return(mata)
  }) #end of function
  
  #CODING VARIANTS
  
  TACDT<-reactive({
    taC<-GET_TACDT(data=sampleDT(), #call function
                   matA=MATA())
    return(taC)
  }) #end of function
  
  
  
  
  
  
  
  ###DATATABLES###
  
  #SAMPLE DATATABLE (YOUR RAW VCF FILE)
  
  output$table1<- renderDataTable({
    data2()
  })#end of function
  
  #GENOTYPES table
  
  output$table2 <- renderDataTable(
    GT()
  ) #end of function
  
  #VARIANTS datatable
  
  observe({output$table3 <- renderDataTable(Variants_tab())
  }) #end of function
  
  #Mutation types datatable
  
  observe({output$table4 <- renderDataTable(Mutations())
  }) #end of function
  
  #Transition datatable
  
  observe({output$table5 <- renderDataTable(Transition())
  }) #end of function
  
  #Patterns datatable
  
  observe({output$table6 <- renderDataTable(PATTERNS())
  }) #end of function
  
  #Annotation datatable
  
  observe({output$table7 <- renderDataTable(ANNOTATIONS())
  })#end of function
  
  #Reactome pathways datatable
  
  observe({output$table8 <- renderDataTable(PATHWAYS())
  })#end of function
  
  
  
  
  ###PLOTS###
  
  #PLOT DEPTH
  
  observe({
    if (is.null(input$inSelect)) #if no sample selected
      output$plot2 <- renderPlot({
        plot_histogram(Depth(),colnames(Depth()[1]))
      }) 
    else #for the selected sample
      output$plot2 <- renderPlot({ 
        plot_histogram(Depth(),
                       colnames(Depth()[as.numeric(input$inSelect)]))
      })}) #end of function
  
  #plot GENOTYPES
  
  observe({ output$plot1 <- renderPlot({
    Bar_plot1(GT())
  })}) #end of function
  
  #plot GT BY SAMPLE
  
  observe({
    if (is.null(input$inSelect)) #if no sample selected
      output$plot3 <- renderPlot({
        Bar_plot2(GTsamples(),
                  colnames(GTsamples()[1]))
      }) 
    else #for the selected sample
      output$plot3 <- renderPlot({
        Bar_plot2(GTsamples(),
                  colnames(GTsamples()[as.numeric(input$inSelect)]))
      })}) #end of function
  
  #plot the mutations
  
  observe ({ output$plot4 <- renderPlot({
    Bar_plot1(Mutations())
  })}) #end of function
  
  #plot Transitions
  
  observe({output$plot5 <- 
    renderPlot({Bar_plot3(Transition())
    })}) #end of function
  
  #Plot Ti/Tv ratio
  
  observe({output$plot6 <- renderPlot({Bar_plot3(Ratio())
  })}) #end of function
  
  #plot the patterns ratio
  
  observe({output$plot7 <- renderPlot({
    Bar_plot3(PATTERNS())
  }) }) #end of function
  
  #other patterns plot
  
  observe({output$plot8<-renderPlot({
    Bar_plot4(PATTERNS())
  }) }) #end of function
  
  #plot variations
  
  observe({output$plot9<-renderPlot({
    Bar_plot5(VARIATIONS())
  })})#end of function
  
  #plot SNP types
  
  observe({output$plot10<-renderPlot({
    Bar_plot3(TACDT())
  })}) #end of function
  
  #plot SQ before filtering
  
  observe({
    output$plot11 <- renderPlot({
      plot(HB(), main = "2D Histogram: AF vs. SQ", xlab = "SQ", ylab = "AF")
    })})
  
  #plot SQ after filtering
  
  observe({
    output$plot12 <- renderPlot({
      plot(HB2(), main = "2D Histogram: AF vs. SQ", xlab = "SQ", ylab = "AF")
    })})
  
  #plot GENES AF 
  
  observe({
    output$plot13 <- renderPlot({
      Bar_plot6(FILTERED_ANNOTATIONS())
    })})
  
  #lolliplot 
  
  observe({
    output$plot14 <- renderPlot({
      plot_lolliplot(FILTERED_ANNOTATIONS(),"start", "AF", "hgnc_symbol")
    })})
  

  
}#end of server
  
##UI##

ui <- fluidPage(
  tabsetPanel(
    tabPanel("MENU", fluid = TRUE, #MAIN TAB
             h3("VCF XPLORR"),
             fileInput("fileInput", "Choose a VCF file",
                       accept = c(".vcf", ".vcf.gz"),
                       multiple = FALSE),
             p("Please select the sample of your choice"),
             selectInput("inSelect", "Select your sample : 1 is for N and 2 for T",
                         c(c("1", "2"))),
             sliderInput("SQValue", "Pick the value of SQ to filter your VCF:",
                                     min = 0, max = 200,
                                     value = 10),
             selectInput("CHRSelect", "Select a chromosome : ",
                         c(c("chr1:", "chr2:", "chr3:", "chr4:", "chr5:",
                             "chr6:", "chr7:", "chr8:", "chr9:", "chr10:", 
                             "chr11:", "chr12:", "chr13:", "chr14:", "chr15:", 
                             "chr16:", "chr17:", "chr18:", "chr19:", "chr19:",
                             "chr20:","chr21:","chr22:","chrX:","chrY:"
                             )))
                         ), #end of MAIN TAB
    tabPanel("VCF TABLES", fluid = TRUE, #VCF TAB
             p("Write the number of sequence you want to be displayed"),
             numericInput("numvar", "num2chk", value = 50, min = 50, max = 500, step = 10), #end
             h3("Your Raw VCF File"),
             DTOutput("table1"),
             h3("Your sample vcf"),
             DTOutput("table3"),
             DTOutput("table10")
             ), #end of VCF TAB
    tabPanel("QUALITY", fluid = TRUE, 
             h3("Depth Plot"),
             plotOutput("plot2"),
             h3("SQ/AF BEFORE SQ Filtering"),
             plotOutput("plot11"),
             h3("SQ/AF AFTER SQ Filtering"),
             plotOutput("plot12")
             ), #end of QUALITY TAB
    tabPanel("GENOTYPE", fluid = TRUE, #GENOTYPE TAB
             h3("HERE IS DISPLAYED YOUR GENOTYPE DATA"),
             h3("Genotypes table"),
             DTOutput("table2"),
             h3("Genotype Plot"),
             plotOutput("plot1"),
             h3("Sample GT plot"),
             plotOutput("plot3")
             ), #end of GENOTYPE TAB
    tabPanel("PATTERNS", fluid = TRUE, #PATTERNS TAB
             h3("HERE IS DISPLAYED THE PRESENCE OF PATTERNS WITHIN YOUR SAMPLE"),
             h3("APOBEC patterns"),
             DTOutput("table6"),
             h3("APOBEC patterns plots"),
             plotOutput("plot8"),
             h3("Trinucleotide pattern plots"),
             plotOutput("plot7")
             ),#end of PATTERNS TAB
    tabPanel("VARIANT TYPES", fluid = TRUE, #VARIANT TYPES TAB
             h3("HERE IS DISPLAYED THE VARIANT TYPES DATA"),
             h3("Variant types plot"),
             plotOutput("plot4"),
             h3("Variant types datatable"),
             DTOutput("table4")
             ),#end of VARIANT TYPES TAB
    tabPanel("SNP", fluid = TRUE, #SNP TAB
             h3("HERE ARE DISPLAYED SNPs INFORMATIONS"),
             h3("Transition and Transversion"),
             plotOutput("plot5"),
             h3("Transition/Transversion ratio"),
             plotOutput("plot6"),
             h3("Types of transition and Transversion"),
             DTOutput("table5"),
             h3("Variations in dbSNP database"),
             plotOutput("plot9"),
             h3("Types of SNPs"),
             plotOutput("plot10")
             ),#end of SNP TAB
    tabPanel("ANNOTATIONS", fluid = TRUE, #ANNOTATION TAB
             h3("HERE IS DISPLAYED THE ANNOTATION DATA"),
             h3("Genes allele frequency plot for the target chromosome"),
             plotOutput("plot13"),
             h3("Lolliplot for the target chromosome (Allele frequency)"),
             plotOutput("plot14"),
             h3("Annotated dataframe"),
             DTOutput("table7"),
             h3("Reactome Pathways dataframe"),
             DTOutput("table8")
             )#end of ANNNOTATION TAB
  ) #end of fluid page
)#end of UI


##RUNAPP##

shinyApp(ui = ui, server = server)