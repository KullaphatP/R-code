#install.packages("outliers")
#library(outliers)
for (n in c(5,10,15,20,30,50,100)) {
  for (alpha in c(0.01,0.05,0.1)) {
    
    powerQ<-0
    powerG<-0
    powerFer<-0
    powerWor<-0
    Gcri<-((n-1)/sqrt(n)) * sqrt( (qt((alpha/(2*n)),n-2)^2)/((n-2)+(qt((alpha/(2*n)),n-2,lower.tail = T)^2)) ) 
    
    Zcri<-qnorm(alpha,lower.tail = F)
    
    Qcri.01<-c(0.7819,0.526,0.4385,0.3922,0.3425,0.296,0.2498)
    Qcri.05<-c(0.6423,0.4122,0.3389,0.3005,0.2594,0.2216,0.1846)
    Qcri.1<-c(0.5578,0.3492,0.2848,0.2511,0.2155,0.1829,0.1512)
    ifelse(alpha==0.01,Qcri<-Qcri.01,ifelse(alpha==0.05,Qcri<-Qcri.05,Qcri<-Qcri.1))
    ifelse(n==5,Qcri<-Qcri[1],ifelse(n==10,Qcri<-Qcri[2],ifelse(n==15,Qcri<-Qcri[3],ifelse(n==20,Qcri<-Qcri[4],Qcri<-Qcri[5]))))
    
    
    
    fer.01<-c(1.3,1.3,1.2,1.1,0.99)
    fer.05<-c(1.0,0.9,0.8,0.8,0.66)
    fer.1<-c(1.0,0.9,0.8,0.8,0.66)
    ifelse(alpha==0.01,Fcri<-fer.01,ifelse(alpha==0.05,Fcri<-fer.05,Fcri<-fer.1))
    ifelse(n==5,Fcri<-Fcri[1],ifelse(n==10,Fcri<-Fcri[2],ifelse(n==15,Fcri<-Fcri[3],ifelse(n==20,Fcri<-Fcri[4],Fcri<-Fcri[5]))))
    
    
    for (i in 1:10000) {
      set.seed(i)
      demo.data <- sort(c(rnorm(n-1,0,1)));sort(demo.data)
      Q3<-demo.data[floor((3*n)/4)]+(((3*(n))/4)-floor((3*(n))/4))*(demo.data[ceiling((3*(n))/4)]-demo.data[floor((3*(n))/4)])
      Q1<-demo.data[floor(n/4)]+((n/4)-floor(n/4))*(demo.data[ceiling(n/4)]-demo.data[floor(n/4)])
      IQR<-Q3-Q1
      Ufence4.5<-Q3+(4.5*IQR)
      Ufence3<-Q3+(3*IQR)
      outlier<-sample((seq(Ufence3,Ufence4.5,0.001)),1)
      dat<-c(demo.data,outlier)
      y<-sort(dat)
      
      Q<-(y[n]-y[n-1])/(max(y)-min(y)) #Dixon test 
      #####
      G<-((y[n])-mean(y))/sd(y) #grubbs test
      #####
      b<-0
      for (a in 1:n) {
        b<-b+(y[a]-mean(y))^3
      }
      c<-b
      e<-0
      for (d in 1:n) {
        e<-e+(y[d]-mean(y))^2
      }
      f<-e^(3/2)
      Fer <-(sqrt(n)*c)/f #ferguson
      #####
      p<-0
      for (q in 1:n) {
        p<-p+(y[q]-mean(y))^3
      }
      up<-p/(n*var(y))
      low<-sqrt((((n-1)/n)^3)*((6*(n-2))/((n+1)*(n+3))) )
      Wor<-up/low # worapan
      #####
      
      ifelse(Q>Qcri,powerQ<-powerQ+1,powerQ<-powerQ+0)
      ifelse(G>Gcri,powerG<-powerG+1,powerG<-powerG+0)
      ifelse(Fer>Fcri,powerFer<-powerFer+1,powerFer<-powerFer+0)
      ifelse(Wor>Zcri,powerWor<-powerWor+1,powerWor<-powerWor+0)
      
    }
    newdf = data.frame(n,alpha,powerQ,powerG,powerFer,powerWor)
    names(newdf)=c("n","alpha", "Dixon","Grubbs","Ferguson","Worrapan")
    print(newdf)
  }
}
