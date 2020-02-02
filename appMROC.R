library(shiny)
library(boot)
library(shinythemes)
####function to define stepwise procedure
stepwise<-function(z,D,H){
  precision<-function(x,D,H,pF){
    muD<-(colMeans(subset(x, x[,1]==D))[-1])
    sigmaD<-(cov(subset(x, x[,1]==D))[2:ncol(x),2:ncol(x)])
    muH<-(colMeans(subset(x, x[,1]==H))[-1])
    sigmaH<-(cov(subset(x, x[,1]==H))[2:ncol(x),2:ncol(x)])
    sigma<-cov(x[2:ncol(x)])
    
    MROC<-function(x,a){
      b<-solve((a*sigmaD)+((1-a)*sigmaH))%*%(muD-muH)
      AUC<-pnorm((t(b)%*%(muD-muH))*((t(b)%*%(sigmaD+sigmaH)%*%b)^(-1/2)))
      c<-((t(b)%*%muD)*((t(b)%*%(sigmaH)%*%b)^(1/2))+(t(b)%*%muH)*((t(b)%*%(sigmaD)%*%b)^(1/2)))/(((t(b)%*%(sigmaD)%*%b)^(1/2))+((t(b)%*%(sigmaH)%*%b)^(1/2)))
      TPR<-pnorm(((t(b)%*%muD)-c)/((t(b)%*%(sigmaD)%*%b)^(1/2)))
      FPR<-1-pnorm((c-(t(b)%*%muH))/((t(b)%*%(sigmaH)%*%b)^(1/2)))
      J<-TPR-FPR
      return(c(a,AUC,c,TPR,FPR,J))
    }
    Res<-data.frame(rbind(MROC(x,0.1),MROC(x,0.2),MROC(x,0.3),MROC(x,0.4),MROC(x,0.5),MROC(x,0.6),MROC(x,0.7),MROC(x,0.8),MROC(x,0.9)))
    ResultOptimal<-Res[which(Res[,6]==max(Res[,6])),]
    
    topt<-Res[,1][which(Res[,6]==max(Res[,6]))]
    b<-function(x,a){
      solve((a*sigmaD)+((1-a)*sigmaH))%*%(muD-muH)
    }
    
    U<-c(0)
    for(j in 1:nrow(x))
    {
      U[j]<-as.vector(b(x,topt))%*%t(x[j,-1])
    }
    s<-c(0)
    for(j in 1:nrow(x))
    {
      if (x[j,1]==H) s[j]=(-(table(x[,1])[2])/nrow(x)) else s[j]=(table(x[,1])[1]/nrow(x))
    }
    linear<-lm(s~as.matrix(x[-1])-1,x)
    
    V<-c(0)
    for(j in 1:nrow(x))
    {
      V[j]<-as.vector(linear$coefficients)%*%t(x[j,-1])
    }
    rp<-cor(U,V)
    if (pF==0) R<-sqrt(summary(linear)$r.squared)
    if (pF==1) R<-sqrt(1-(det(cor(x[-1]))/det(cor(x[-(1:2)]))))
    r<-sqrt((rp^2*(1-R^2))/(1-(R^2*rp^2)))
    
    lambda2<-sum(s^2)
    Fratio<-((nrow(x)-ncol(x[-1])-1)*(R^2*(1-r^2)))/((ncol(x[-1])-1)*(1-R^2*(1-r^2)))
    pvalue<-1-pf(Fratio,ncol(x[-1])-1,nrow(x)-ncol(x[-1])-1)
    return(c(Fratio,pvalue))
  }
  alphain<-0.05
  alphaout<-0.10
  com<-as.matrix(combn(2:ncol(z),2))
  Result<-matrix(c(0),nrow=ncol(com),ncol=3)
  for(i in 1:ncol(com))
  {
    x<-data.frame(z[,1],z[,com[1:nrow(com),i]])
    m<-paste(com[,i],collapse=",")
    Result[i,]<-c(m,precision(x,D,H,0))
  }
  colnames(Result)=c("subset","F ratio","p-value")
  #noquote(Result)
  ##variables in the model - vim
  vim<-com[,which(Result[,2]==max(as.numeric(Result[,2])))]
  #vim
  
  repeat
  {
    vti<-(2:ncol(z))[! (2:ncol(z)) %in% vim]
    Result1<-matrix(c(0),nrow=length(vti),ncol=3)
    for(i in vti){
      x<-data.frame(z[,1],z[,i],z[,vim])
      Result1[which(vti==i),]<-c(i,precision(x,D,H,1))
    }
    #print(Result1)
    #Result1[which(Result1[,2]==max(subset(Result1,Result1[,3]<alphain)[,2])),1]
    if (nrow(subset(Result1,Result1[,3]<alphain))==0) break
    vim<-c(Result1[which(Result1[,2]==max(subset(Result1,Result1[,3]<alphain)[,2])),1],vim)
    #print(c("vim = ",vim),quote=FALSE)
    
    Result1<-matrix(c(0),nrow=length(vim),ncol=3)
    for(i in vim){
      x<-data.frame(z[,1],z[,i],z[,vim[! vim %in% i]])
      Result1[which(vim==i),]<-c(i,precision(x,D,H,1))
    }
    #print(Result1)
    #Result1[which(Result1[,2]==min(subset(Result1,Result1[,3]>alphaout)[,2])),1]
    if (nrow(subset(Result1,Result1[,3]>alphaout))==0) {vim<-vim} else {
      vim<-vim[! vim %in% (Result1[which(Result1[,2]==min(subset(Result1,Result1[,3]>alphaout)[,2])),1])]
      #print(c("vim = ",vim),quote=FALSE)
    }
  }
  #print(c("vim = ",vim),quote=FALSE)
  
  #### data after variable selection
  x<-data.frame(z[,1],z[,vim])
  return(x)
}

ui<-navbarPage("MROC Model",
               tabPanel("Data Import",
                        sidebarLayout(sidebarPanel(fileInput(inputId="data",label="Choose CSV file",multiple=FALSE,accept=c("text/csv","text/comma-seperated-values,text/plain",".csv"),width=NULL,buttonLabel="Browse",placeholder="No file selected"),
                                                   tags$hr(),
                                                   h5(helpText("Select the header if data contains column titles")),
                                                   checkboxInput(inputId="header",label="Header",value=FALSE),
                                                   radioButtons(inputId = 'sep', label = 'Separator', 
                                                                choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ',')
                        ),
                        mainPanel(uiOutput("dataframe")),position=c("left","right"),fluid=TRUE
                        )),
               tabPanel("MROC model",
                        sidebarLayout(sidebarPanel(
                          
                          uiOutput("status1"),
                          uiOutput("D1"),
                          uiOutput("H1"),
                          selectInput(inputId = "model",label = "Select Model",choices = c("Full Model"="fm","Stepwise Model"="sm","Custom Model"="cm")),
                          uiOutput("Independents"),
                          actionButton(inputId="click",label="Submit")
                        ),
                        mainPanel(tabsetPanel(type="tabs",
                                    tabPanel("MROC measures",tableOutput("measures")),
                                    tabPanel("MROC Coefficients and Precision",tags$h3(helpText("MROC Coefficients")),
                                             verbatimTextOutput("coefficients"),
                                             tags$h3(helpText("Precision")),
                                             verbatimTextOutput("precision")),
                                    tabPanel("MROC Curve",plotOutput("curve",width="400px",height="400px"))),fluid=TRUE
                        ))))

server<-function(input,output){
  file<-reactive({
    x<-input$data
    if(is.null(x)){return()}
    read.table(file=x$datapath,sep=input$sep,header=input$header)
  })
  
  output$x1<-renderTable({
    if(is.null(file())){return()}
    file()
  })
  output$dataframe<-renderUI({
    tableOutput("x1")
  })
  
  output$status1<-renderUI({
    selectInput("status","Select Dependent variable",choices=names(file()),multiple=FALSE)
  })
  output$D1<-renderUI({
    selectInput("D","Select Group 1 value",choices=unique(file()[,input$status]),multiple=FALSE)
  })
  output$H1<-renderUI({
    selectInput("H","Select Group 2 value",choices=unique(file()[,input$status]),multiple=FALSE)
  })
  output$Independents<-renderUI({
    if(input$model=="cm"){
      checkboxGroupInput("Variables","Select Independent Variables",choices=names(file()))
    }
  })
  
  y1<-reactive({
    if(input$model=="fm"){
      y<-file()
      
      file()[,c(input$status,setdiff(names(file()),input$status))]
    }
    else if(input$model=="sm"){
      y<-file()
      z<-file()[,c(input$status,setdiff(names(file()),input$status))]
      #z<-data.frame(cbind(file()[input$status],file()[!names(file()) %in% input$status]))
      stepwise(z,input$D,input$H)
    }
    else if(input$model=="cm"){
      data.frame(cbind(file()[,c(input$status,input$Variables)]))
    }
  })
  
  values<-eventReactive(input$click,{
    y<-(y1())
    muD1<-colMeans(subset(y[-1],y[,1]==input$D))
    muH1<-colMeans(subset(y[-1],y[,1]==input$H))
    sigmaD1<-cov(subset(y[-1],y[,1]==input$D))
    sigmaH1<-cov(subset(y[-1],y[,1]==input$H))
    list(muD1=muD1,muH1=muH1,sigmaD1=sigmaD1,sigmaH1=sigmaH1)
  })
  
  meas<-reactive({
    y<-y1()
    muD<-values()$muD1
    muH<-values()$muH1
    sigmaD<-values()$sigmaD1
    sigmaH<-values()$sigmaH1
    
    MROC2<-function(y,d){
      x=y[d,]
      res<-matrix(c(0),nrow=9,ncol=6)
      
      for(a in seq(0.1,0.9,0.1)){
        MROC<-function(a){
          b<-solve((a*sigmaD)+((1-a)*sigmaH))%*%(muD-muH)
          AUC<-pnorm((t(b)%*%(muD-muH))*((t(b)%*%(sigmaD+sigmaH)%*%b)^(-1/2)))
          c<-((t(b)%*%muD)*((t(b)%*%(sigmaH)%*%b)^(1/2))+(t(b)%*%muH)*((t(b)%*%(sigmaD)%*%b)^(1/2)))/(((t(b)%*%(sigmaD)%*%b)^(1/2))+((t(b)%*%(sigmaH)%*%b)^(1/2)))
          TPR<-pnorm(((t(b)%*%muD)-c)/((t(b)%*%(sigmaD)%*%b)^(1/2)))
          FPR<-1-pnorm((c-(t(b)%*%muH))/((t(b)%*%(sigmaH)%*%b)^(1/2)))
          J<-TPR-FPR
          return(c(a,AUC,TPR,FPR,c,J))
        }
        res[10*a,]<-(MROC(a))
      }
      res
      ResultOpt<-res[which(res[,6]==max(res[,6])),]
      ResultOpt
      return(ResultOpt)
    }
    
    bootMROC<-boot(y,MROC2,100)
    
    Results<-data.frame(rbind(c("t",round(bootMROC$t0[1],4),round(mean(bootMROC$t[,1]),4),round(mean(bootMROC$t[,1])-qnorm(0.975)*sqrt(var(bootMROC$t[,1])),4),round(mean(bootMROC$t[,1])+qnorm(0.975)*sqrt(var(bootMROC$t[,1])),4),"-","-"),
                              c("AUC",round(bootMROC$t0[2],4),round(mean(bootMROC$t[,2]),4),round(mean(bootMROC$t[,2])-qnorm(0.975)*sqrt(var(bootMROC$t[,2])),4),round(mean(bootMROC$t[,2])+qnorm(0.975)*sqrt(var(bootMROC$t[,2])),4),round((mean(bootMROC$t[,2])-0.5)/sqrt(var(bootMROC$t[,2])),4),round(1-pnorm((mean(bootMROC$t[,2])-0.5)/sqrt(var(bootMROC$t[,2]))),4)),
                              c("TPR",round(bootMROC$t0[3],4),round(mean(bootMROC$t[,3]),4),round(mean(bootMROC$t[,3])-qnorm(0.975)*sqrt(var(bootMROC$t[,3])),4),round(mean(bootMROC$t[,3])+qnorm(0.975)*sqrt(var(bootMROC$t[,3])),4),round((mean(bootMROC$t[,3])-0.5)/sqrt(var(bootMROC$t[,3])),4),round(1-pnorm((mean(bootMROC$t[,3])-0.5)/sqrt(var(bootMROC$t[,3]))),4)),
                              c("FPR",round(bootMROC$t0[4],4),round(mean(bootMROC$t[,4]),4),round(mean(bootMROC$t[,4])-qnorm(0.975)*sqrt(var(bootMROC$t[,4])),4),round(mean(bootMROC$t[,4])+qnorm(0.975)*sqrt(var(bootMROC$t[,4])),4),"-","-"),
                              c("c",round(bootMROC$t0[5],4),round(mean(bootMROC$t[,5]),4),round(mean(bootMROC$t[,5])-qnorm(0.975)*sqrt(var(bootMROC$t[,5])),4),round(mean(bootMROC$t[,5])+qnorm(0.975)*sqrt(var(bootMROC$t[,5])),4),"-","-"),
                              c("J",round(bootMROC$t0[6],4),round(mean(bootMROC$t[,6]),4),round(mean(bootMROC$t[,6])-qnorm(0.975)*sqrt(var(bootMROC$t[,6])),4),round(mean(bootMROC$t[,6])+qnorm(0.975)*sqrt(var(bootMROC$t[,6])),4),"-","-")
    ))
    colnames(Results)=c("Measures","Original Value","Bootstrapped Estimate","Lower Limit","Upper Limit","Z value","p-value")
    return(Results)
  })
  
  coeff<-reactive({
    y<-(y1())
    muD<-values()$muD1
    muH<-values()$muH1
    sigmaD<-values()$sigmaD1
    sigmaH<-values()$sigmaH1
    
    topt<-as.numeric(paste(meas()[1,2]))
    bopt<-solve((topt*sigmaD)+((1-topt)*sigmaH))%*%(muD-muH)
    return(bopt)
  })
  
  output$measures<-renderTable(
    meas()[2:nrow(meas()),],rownames=FALSE,quoted=F
  )
  
  output$coefficients<-renderPrint({
    print(round(coeff(),4))
  })
  
  output$precision<-renderPrint({
    y<-(y1())
    
    bopt<-coeff()
    
    #precision of curve
    U<-c(0)
    for(j in 1:nrow(y))
    {
      U[j]<-as.vector(bopt)%*%t(y[j,-1])
    }
    s<-c(0)
    for(j in 1:nrow(y))
    {
      if (y[j,1]==input$H) s[j]=(-(table(y[,1])[2])/nrow(y)) else s[j]=(table(y[,1])[1]/nrow(y))
    }
    linear<-lm(s~as.matrix(y[-1])-1,y)
    
    V<-c(0)
    for(j in 1:nrow(y))
    {
      V[j]<-as.vector(linear$coefficients)%*%t(y[j,-1])
    }
    rp<-cor(U,V)
    R<-sqrt(summary(linear)$r.squared)
    r<-sqrt((rp^2*(1-R^2))/(1-(R^2*rp^2)))
    
    lambda2<-sum(s^2)
    Fratio<-((nrow(y)-ncol(y[-1])-1)*(R^2*(1-r^2)))/((ncol(y[-1])-1)*(1-R^2*(1-r^2)))
    pvalue<-1-pf(Fratio,ncol(y[-1])-1,nrow(y)-ncol(y[-1])-1)
    
    cat(c("F-ratio=",round(Fratio,4),", p-value=",round(pvalue,4)),"\n")
  })
  
  output$curve<-renderPlot({
    
    y<-(y1())
    muD<-values()$muD1
    muH<-values()$muH1
    sigmaD<-values()$sigmaD1
    sigmaH<-values()$sigmaH1
    
    bopt<-coeff()
    
    #MROC curve construction
    ROCTPR<-NULL
    ROCFPR<-NULL
    for(j in 1:nrow(y))
    {
      ROCTPR[j]<-pnorm(((t(bopt)%*%muD)-(t(bopt)%*%t(y[j,2:ncol(y)])))/((t(bopt)%*%(sigmaD)%*%bopt)^(1/2)))
      ROCFPR[j]<-1-pnorm(((t(bopt)%*%t(y[j,2:ncol(y)]))-(t(bopt)%*%muH))/((t(bopt)%*%(sigmaH)%*%bopt)^(1/2)))
    }
    coordinates<-cbind(ROCFPR,ROCTPR)
    
    diag<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
    
    plot(diag,diag,type="l",lwd=2,col="blue",lty=2,0:1,0:1,xlab="1-Specificity",ylab="Sensitivity",main="MROC Curve")
    lines(coordinates[order(coordinates[,1]),],type="l",lwd=2,lty=1,col="red")
    legend(0.55,0.35,c("Diagonal","MROC"),lwd=2,col=c("blue","red"),lty=c(2,1))
  })
}

shinyApp(ui=ui,server=server)