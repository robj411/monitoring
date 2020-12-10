library( lubridate )
library( shiny )
library(shinyjs)
library(grid)
#library(ape)
library(DT)
library(scales)
library(ggtree)
library( ggplot2 )
library(plotly)
library(dplyr)
library(RColorBrewer)
library(plotrix)

#scale_fill_discrete <- function(...) {
#  scale_fill_brewer(..., palette="Accent")
#}
cllist <- readRDS('sequencesummary.Rds')
genes <- read.csv('Genes.CSV',stringsAs=F)




trends <- lapply(cllist,function(x) 
  lapply(x,function(y){
    xx<-rowSums(y$trend,dims=2); 
    total=colSums(xx[1:5,]); 
    xx[,sort(total,decreasing = T,index.return=T)$ix[1:2]]
  }))
median_growths <- lapply(trends,function(y)lapply(y,function(x){
  n1 <- x[,1]
  n2 <- x[,2]
  ll <- length(n1)
  pos <- n1>0 
  pos2 <- n2>0 
  n1 <- n1[pos]
  n2 <- n2[pos]
  rate <- 0
  if(sum(pos)>4){
    suppressWarnings(l1 <- min(c((ll-13):ll)[pos2[(ll-13):ll]],na.rm=T))
    odds <- approx(x=c(1:ll)[pos],y=n2/n1,xout=c(l1,ll-3))$y
    logodds <- diff(log(odds))
    if(!is.na(logodds)) rate <- logodds
  }
  round(rate,2)
  }
  ))
summary_table <- as.data.frame(do.call(rbind,
                         lapply(1:length(trends),
                                function(g) 
                                  t(sapply(1:length(trends[[g]]),
                                         function(i) {
                                           pos <- gsub('p','',names(trends[[g]])[i])
                                           frst <- colnames(trends[[g]][[i]])[1]
                                           scnd <- colnames(trends[[g]][[i]])[2]
                                           lh <-  cllist[[g]][[i]]$lighthouse[c(frst,scnd),]
                                           lhodds <- lh[2,]/lh[1,]
                                           lhoddsratio <- round(log(lhodds[1]/lhodds[2]),2)
                                           gen <-  cllist[[g]][[i]]$source_sex[c(frst,scnd),c('F','M')]
                                           genodds <- gen[2,]/gen[1,]
                                           genoddsratio <- round(log(genodds[1]/genodds[2]),2)
                                           c(names(trends)[g],pos,frst,scnd,
                                                       paste0(frst,pos,scnd),median_growths[[g]][[i]],
                                             lhoddsratio,genoddsratio)
                                }))
                                )
                         ),stringsAsFactors=F)
colnames(summary_table) <- c('Gene','Position','Original variant','New variant','Label','Mean recent growth','Hospital log odds ratio','Female log odds ratio')
summary_table$`Mean recent growth` <- as.numeric(summary_table$`Mean recent growth`)
for(x in c('Hospital log odds ratio','Female log odds ratio'))
  summary_table[[x]][grepl('Na|Inf',summary_table[[x]])] <- ''


for(g in 1:length(trends)){
  growths <- subset(summary_table,Gene==gene_names[g])$`Mean recent growth`
  neworder <- order(growths,decreasing=T)
  trends[[g]] <- trends[[g]][neworder]
  cllist[[g]] <- cllist[[g]][neworder]
}

regions <- dimnames(cllist[[1]][[1]]$trend)$NUTS1
regions <- regions[regions!='']
cols <- rainbow(length(regions))
reg.options <- gsub('_',' ',regions)
reg.colors <- cols[1:length(regions)]
variant_names <- names(cllist[[2]])
gene_names <- names(cllist)

shiny::shinyServer(function(input, output, session) {
  
  ## initialise checkboxes
  reg.fun <- lapply(reg.options,function(o)
    tags$div(
      HTML(paste(tags$span('―',style = paste0('font-weight: bold; color: ', reg.colors[which(reg.options == o)],';')),o, sep = " "))
      , encoding = 'UTF-8'
    )
  )
  
  updateSelectInput(session, inputId='ti_regions', label = 'Regions', 
                    choices = c('All',regions),#reg.fun(),#
                    #choiceNames = c('All',reg.options),#reg.fun,#names(parms$region_lineage),#
                    #choiceValues = c('All',regions),#names(parms$region_lineage), #reg.colors,#
                    selected = 'All')
  updateCheckboxGroupInput(session, inputId='ti_region2', label = 'Regions', 
                    #choices = c('All',regions),#reg.fun(),#
                    choiceNames = reg.fun,#reg.fun,#names(parms$region_lineage),#
                    choiceValues = regions,#names(parms$region_lineage), #reg.colors,#
                    selected = NULL)
  
  updateSelectInput(session, inputId='ti_genes', label='Gene', choices = gene_names, selected = gene_names[2])
  updateSelectInput(session, inputId='ti_variants', label='Position', choices = variant_names, selected = 'p501')
  
  
  output$legend <- renderPlot({
    par(mar=c(0,0,0,0)); plot.new()#plot(c(0,1),c(0,1)); 
  })
  
  ## observe input events
  
  
  # update variants
  observeEvent({
    input$ti_genes
  },{
    if(input$ti_genes!=''){
      variant_names <<- names(cllist[[which(gene_names==input$ti_genes)]])
      updateSelectInput(session, inputId='ti_variants', label='Position', choices = variant_names, selected = variant_names[1])
    }
  })
  
  ## plot
  observeEvent({
    input$ti_regions
    input$ti_region2
    input$ti_lh
    input$ti_variants
  },{
    if(input$ti_variants%in%variant_names){
      #    if(input$ti_genes!=''){
      g <- 2; i <- 1
      g <- which(gene_names==input$ti_genes)
      i <- which(variant_names==input$ti_variants)
      lab <- input$ti_lh
      if(lab=='Both'){
        subtotal <- rowSums(cllist[[g]][[i]]$trend,dims=3)
      }else{
        subtotal <- cllist[[g]][[i]]$trend[,,,which(c('Hospital','Lighthouse')==lab)]
      }
      reg <- input$ti_regions
      if(reg=='All'){
        total <- rowSums(subtotal,dim=2)
      }else{
        total <- subtotal[,,which(dimnames(subtotal)$NUTS1==reg)]
      }
      varindices <- sort.int(colSums(total),decreasing=T,index.return=T)$ix[1:2]
      total <- total[,varindices]
      new_variant <- which.min(total[1,])
      n1 <- total[,new_variant]
      n2 <- total[,-new_variant]
      value <- names(new_variant)
      variable <- paste0(names(cllist)[g],gsub('p','',names(cllist[[g]])[i]))
      n <- rowSums(total)
      output$trend <- renderPlot({
        suppressWarnings(suppressMessages(ggplot(  aes( x = as.numeric(week), y = variant ) , 
                                                   data= data.frame(week=rownames(total),variant=n1,stringsAsFactors = F) ) + 
                                            geom_point ()  + 
                                            theme_minimal() + 
                                            theme( legend.pos='',
                                                   axis.text=element_text(size=14),axis.title=element_text(size=14)) + 
                                            #geom_smooth(method=stats::loess, method.args=list(span = 1)) + 
                                            xlab('') + 
                                            ggtitle( paste0('Frequency of ', variable, '=', value) )))
      })
      
      ## second plot
      reg2 <- input$ti_region2
      if(is.null(reg2)){
        hide('trend2')
      }else{
        show('trend2')
        if(!reg%in%reg2&reg!='All'){
          reg2 <- c(reg,reg2)
          updateCheckboxGroupInput(session, inputId='ti_region2', label = 'Regions', 
                                   choiceNames = reg.fun,#reg.fun,#names(parms$region_lineage),#
                                   choiceValues = regions,#names(parms$region_lineage), #reg.colors,#
                                   selected = reg2)
        }
        total <- subtotal[,,which(dimnames(subtotal)$NUTS1%in%reg2)]
        if(length(dim(total))>2){
          t1 <- data.frame(week=rownames(total),variant=total[,varindices[new_variant],1],regionnanme=reg2[1],stringsAsFactors=F)
          for(d in 2:dim(total)[3]){
            td <- data.frame(week=rownames(total),variant=total[,varindices[new_variant],d],regionnanme=reg2[d],stringsAsFactors=F)
            t1 <- rbind(t1,td)
          }
        }else{
          t1 <- data.frame(week=rownames(total),variant=total[,varindices[new_variant]],regionnanme=reg2,stringsAsFactors=F)
        }
        colScale <- scale_colour_manual(name = regions,values = reg.colors)
        output$trend2 <- renderPlot({
          suppressWarnings(suppressMessages(ggplot(t1, aes( x = as.numeric(week), y = variant )  ) + 
                                              geom_line (aes( colour=factor(regionnanme) ),size=1.5)  + #colScale +
                                              theme_minimal() + 
                                              scale_colour_manual(values = reg.colors[match(unique(t1$regionnanme),regions)]) +
                                              theme( legend.pos='',
                                                            axis.text=element_text(size=14),axis.title=element_text(size=14)) + 
                                              #geom_smooth(method=stats::loess, method.args=list(span = 1)) + 
                                              xlab('Week') + ylab('Number') + 
                                              ggtitle( paste0('Frequency of ', variable, '=', value) )))
        })
      }
      
      ## log odds
      p = n1  / n
      v = p * (1 - p) / n 
      weight <- 1 / sqrt( v )
      weight = weight / median( weight,na.rm = T )
      logodds = log( n1 )  - log(n2 )
      if(sum(is.na(weight))>3){
        hide('logodds')
      }else{
        suppressMessages(show('logodds'))
        output$logodds <- renderPlot({
          {
            suppressWarnings(suppressMessages(ggplot(  aes( x = as.numeric(week), y = logodds, size=weight, weight=weight ) , 
                                                       data= data.frame(week=rownames(total),weight=weight,logodds=logodds,stringsAsFactors=F) ) + 
                                                geom_point ()  + 
                                                theme_minimal() + 
                                                theme( legend.pos='',
                                                       legend.title=element_text(size=14), 
                                                       legend.text=element_text(size=12),
                                                       axis.text=element_text(size=14),axis.title=element_text(size=14)) + 
                                                geom_smooth(method=stats::loess, method.args=list(span = 1)) + 
                                                xlab('') + 
                                                ggtitle( paste0('Log odds of ', variable, '=', value) )))
          }
      })
      }

      
      output$start <- renderPlot({
        g <- which(gene_names==input$ti_genes)
        if(input$ti_variants%in%variant_names  ){
          i <- which(variant_names==input$ti_variants)
        value <- names(new_variant)
        reg <- input$ti_regions
        if(reg=='All'){
          sampletimes <- rowSums(cllist[[g]][[i]]$cluster_first_sampled,dim=2)
        }else{
          sampletimes <- cllist[[g]][[i]]$cluster_first_sampled[,,which(dimnames(cllist[[g]][[i]]$cluster_first_sampled)$NUTS1==reg)]
        }
        #sampletimes <- cllist[[g]][[i]]$cluster_first_sampled
        if(!ncol(sampletimes)<max(varindices)){
        vrs <- colnames(sampletimes)
        data <- data.frame(Week=rep(rownames(sampletimes),2),
                           Frequency=c(sampletimes[,varindices[-new_variant]],sampletimes[,varindices[new_variant]]),
                           genotype=rep(c(vrs[varindices[-new_variant]],vrs[varindices[new_variant]]),each=nrow(sampletimes)),
                           stringsAsFactors=F)
        ggplot(  aes( x = as.numeric(Week), y = Frequency, col = genotype ) , 
                 data= data ) + 
          geom_line ()  + 
          theme_minimal() + 
          theme( legend.pos='right',
                 legend.title=element_text(size=14), 
                 legend.text=element_text(size=12),
                 axis.text=element_text(size=14),axis.title=element_text(size=14)) + 
          xlab('Week') + 
          ggtitle( paste0('Frequency time of first cluster sample') )
        }}
      })
    }
    #}
  })
  
  output$variants <- DT::renderDataTable({
    datatable(summary_table,options = list("pageLength" = 500),rownames=F)
  })
  
  output$genes <- DT::renderDataTable({
    datatable(genes,options = list("pageLength" = 500),rownames=F)
  })
  
  
  output 
  
}
)
