#install.packages("outliers")
#library(outliers)

## set dataframe "df"
n<0;alpha<0;powerQ<0;powerG<0;powerFer<0;powerWor<0
df<-data.frame(n,alpha,powerQ,powerG,powerFer,powerWor)

for (alpha in c(0.01,0.05,0.1)) {
for (n in c(5,10,15,20,30,50,100)) {

  powerQ<-0
  powerG<-0
  powerFer<-0
  powerWor<-0
    
  Gcri<-((n-1)/sqrt(n)) * sqrt( (qt((alpha/(2*n)),n-2)^2)/((n-2)+(qt((alpha/(2*n)),n-2,lower.tail = T)^2)) ) 
  #Qcri<-qdixon(alpha,n)
  Qcri.01<-c(0.7819,0.526,0.4385,0.3922,0.3425,0.296,0.2498)
  Qcri.05<-c(0.6423,0.4122,0.3389,0.3005,0.2594,0.2216,0.1846)
  Qcri.1<-c(0.5578,0.3492,0.2848,0.2511,0.2155,0.1829,0.1512)
  ifelse(alpha==0.01,Qcri<-Qcri.01,ifelse(alpha==0.05,Qcri<-Qcri.05,Qcri<-Qcri.1))
  ifelse(n==5,Qcri<-Qcri[1],ifelse(n==10,Qcri<-Qcri[2],ifelse(n==15,Qcri<-Qcri[3],ifelse(n==20,Qcri<-Qcri[4],Qcri<-Qcri[5]))))
  
  Zcri<-qnorm(alpha,lower.tail = F)
    
    fer.01<-c(1.3,1.3,1.2,1.1,0.99)
    fer.05<-c(1.0,0.9,0.8,0.8,0.66)
    fer.1<-c(1.0,0.9,0.8,0.8,0.66)
    ifelse(alpha==0.01,Fcri<-fer.01,ifelse(alpha==0.05,Fcri<-fer.05,Fcri<-fer.1))
    ifelse(n==5,Fcri<-Fcri[1],ifelse(n==10,Fcri<-Fcri[2],ifelse(n==15,Fcri<-Fcri[3],ifelse(n==20,Fcri<-Fcri[4],Fcri<-Fcri[5]))))

    
    for (i in 1:10000) {
      set.seed(i)
      demo.data <- c(rnorm(n-1,0,1))
      Q3<-quantile(demo.data,0.75)
      IQR<-IQR(demo.data)
      Ufence1.5<-Q3+(1.5*IQR)
      Ufence3<-Q3+(3*IQR)
      outlier<-sample((seq(Ufence1.5,Ufence3,0.0001)),1)
      dat<-c(demo.data,outlier)
      
      y<-sort(dat)
      
      Q<-(y[n]-y[n-1])/(max(y)-min(y)) #Dixon right 
      G<-((y[n])-mean(y))/sd(y) #grubbs right test
      #ferguson
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
      
      ##worapan
      p<-0
      for (q in 1:n) {
        p<-p+(y[q]-mean(y))^3
      }
      up<-p/(n*var(y))
      low<-sqrt((((n-1)/n)^3)*((6*(n-2))/((n+1)*(n+3))) )
      Wor<-up/low
      
      ifelse(Q>Qcri,powerQ<-powerQ+1,powerQ<-powerQ+0)
      ifelse(G>Gcri,powerG<-powerG+1,powerG<-powerG+0)
      ifelse(Fer>Fcri,powerFer<-powerFer+1,powerFer<-powerFer+0)
      ifelse(Wor>Zcri,powerWor<-powerWor+1,powerWor<-powerWor+0)
      
    }
    newdf = data.frame(n,alpha,powerQ,powerG,powerFer,powerWor)
    names(newdf)=c("n","alpha", "Dixon","Grubbs","Ferguson","Worrapan")
    print(newdf)
    
    output <- c(n,alpha,powerQ,powerG,powerFer,powerWor)
    
    
    
    df <- rbind(df,output)

  }
}

colnames(df) <- c("n","alpha", "Dixon","Grubbs","Ferguson","Worrapan")
write_xlsx(df,"C:\\THESIS\\result.xlsx")
print(df)


########################################################
