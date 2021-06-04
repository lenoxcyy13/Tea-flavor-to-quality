#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(pheatmap)
library(rpart)
library(rpart.plot)
library(dplyr)
library(plotly)

ui <- fluidPage(
    
    titlePanel("Gapped Histogram for TEAflavor Data"),
    
    sidebarLayout(
        
        sidebarPanel(
            selectInput("trt", 
                        "treatments to analysis", 
                        choices = c("Origin"="Origin",
                                    "Season"="Season",
                                    "Baking"="Baking",
                                    "Variety"="Variety"
                        )),
            selectInput("var", 
                        "Variables to analysis", 
                        choices = c("MOISTURE" = "MOISTURE",
                                    "T-N" = "T.N",
                                    "TFAA" = "TFAA",
                                    "THEANINE" = "THEANINE",
                                    "NDF-ASH" = "NDF.ASH",
                                    "TANNIN" = "TANNIN",
                                    "CATECHIN" = "CATECHIN",
                                    "CAFFEINE" = "CAFFEINE",
                                    "T.V.C" = "T.V.C"
                        )),
            selectInput("var2", 
                        "Another Variable to analysis in scatter plot", 
                        choices = c("MOISTURE" = "MOISTURE",
                                    "T-N" = "T.N",
                                    "TFAA" = "TFAA",
                                    "THEANINE" = "THEANINE",
                                    "NDF-ASH" = "NDF.ASH",
                                    "TANNIN" = "TANNIN",
                                    "CATECHIN" = "CATECHIN",
                                    "CAFFEINE" = "CAFFEINE",
                                    "T.V.C" = "T.V.C"
                        )),
            selectInput("var3", 
                        "Another Variable to analysis in scatter plot", 
                        choices = c("MOISTURE" = "MOISTURE",
                                    "T-N" = "T.N",
                                    "TFAA" = "TFAA",
                                    "THEANINE" = "THEANINE",
                                    "NDF-ASH" = "NDF.ASH",
                                    "TANNIN" = "TANNIN",
                                    "CATECHIN" = "CATECHIN",
                                    "CAFFEINE" = "CAFFEINE",
                                    "T.V.C" = "T.V.C"
                        )),
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 15,
                        value = 4)
        ),
        
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Plot", plotOutput("distPlot"), plotOutput("boxplot")),
                        tabPanel('3D Scatter', plotlyOutput("scatter"), width = "100%"),
                        tabPanel("Tree", plotOutput("Tree"), width = "100%")
            )
        )
    )
)

server <- function(input, output) {
    tea = read.csv('teaFlavor.csv', header = TRUE)
    
    output$boxplot <- renderPlot({
        sel.var = match(input$var,colnames(tea))
        sel.trt = match(input$trt,colnames(tea))
        
        if (sel.trt==4){
            tea <- tea %>% filter(Variety %in% c('V1', 'V2', 'V3', 'V4'))
        }
        boxplot(tea[,sel.var]~factor(tea[,sel.trt]),
                xlab=paste(colnames(tea)[sel.trt]),
                ylab=paste(colnames(tea)[sel.var]))
    })
    output$distPlot <- renderPlot({
        sel.var = match(input$var,colnames(tea))
        sel.trt = match(input$trt,colnames(tea))

        if (sel.trt==4){
            tea <- tea %>% filter(Variety %in% c('V1', 'V2', 'V3', 'V4'))
        }

        data = tea[, sel.var] 
        h = hclust(dist(data))
        cl = cutree(h, input$bins)
        tcl = table(cl)
        mm = rbind(tapply(data, cl, min), tapply(data, cl, max))
        
        par(mfrow=c(1,2))
        # Uniformity of data in each group
        Fn = ecdf(data)
        plot(Fn, col="grey")
        tdata = cbind(data,cl)
        utdata = unique(tdata)
        mutdata = tapply(utdata[,1],utdata[,2],range)
        for (i in 1:length(mutdata)){
            lines(mutdata[[i]],Fn(mutdata[[i]]),
                  col=2,lwd=2)
        }
        
        # draw the histogram with the specified number of bins
        tcls = table(cl, tea[,sel.trt])
        
        plot(NULL, 
             xlim=c(min(data),max(data)),
             ylim=c(0, max(tcl)), 
             ylab="freq.", xlab="value", 
             main="")
        
        for (n in 1:length(unique(tea[,sel.trt]))){
            for (i in 1:ncol(mm)){
                rect(xleft=mm[1, i], 
                     if (n==1)
                     {ybottom = 0}
                     else 
                     {ybottom = sum(tcls[i, 1:(n-1)])},
                     xright=mm[2, i],
                     ytop=sum(tcls[i, 1:n]),
                     col = n)
                
                
                
            }}
        legend("topright",colnames(tcls),fill=1:n)
        
    })
    output$Tree <- renderPlot({
        sel.var = match(input$var,colnames(tea))
        sel.trt = match(input$trt,colnames(tea))

        if (sel.trt==1){
            tree <- rpart(Origin~., data=tea, cp=.01)
        }
        if (sel.trt==2){
            tree <- rpart(Season~., data=tea, cp=.01)
        }
        if (sel.trt==3){
            tree <- rpart(Baking~., data=tea, cp=.01)
        }
        if (sel.trt==4){
            tea <- tea %>% filter(Variety %in% c('V1', 'V2', 'V3', 'V4'))
            tree <- rpart(Variety~., data=tea, cp=.01)
        }
        
        rpart.plot(tree, box.palette="RdBu", shadow.col="gray", nn=TRUE)
        
    })
    output$scatter <- renderPlotly({
        sel.var = match(input$var,colnames(tea))
        sel.var2 = match(input$var2,colnames(tea))
        sel.var3 = match(input$var3,colnames(tea))
        sel.trt = match(input$trt,colnames(tea))
        
        if (sel.trt==4){
            tea <- tea %>% filter(Variety %in% c('V1', 'V2', 'V3', 'V4'))
        } 
        
        fig <- plot_ly(type = "scatter3d", mode = "markers", data = tea,
                       x = ~tea[, sel.var], y = ~tea[, sel.var2], z = ~tea[, sel.var3],
                       color = ~tea[, sel.trt], width = 900, height = 600)
        fig <- fig %>% layout(scene = list(xaxis = list(title = colnames(tea)[sel.var]),
                                           yaxis = list(title = colnames(tea)[sel.var2]),
                                           zaxis = list(title = colnames(tea)[sel.var3])))
        fig
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



