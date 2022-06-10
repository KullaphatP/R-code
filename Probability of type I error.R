#install.packages("outliers")
#library(outliers)
for (n in c(5,10,15,20,30,50,100)) {
  for (alpha in c(0.01,0.05,0.1)) {
    
    errorQ<-0
    errorG<-0
    errorFer<-0
    errorWor<-0
    errorTM<-0
  
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
    
    set.seed(1234)
    demo.data <- c(rnorm(n,0,1))
    x <-c(demo.data)
    tm = function(x,k){
      k = 0
      n = length(x)
      r = abs(x - mean(x))
      df = data.frame(x,r)
      dfs = df[order(df$r),]
      klarge = c((n-k+1):n)
      subx = dfs$x[-klarge]
      ksub = (subx - mean(subx))**2
      all = (df$x - mean(df$x))**2
      ek = sum(ksub)/sum(all)
    }
    ekstat = tm(x,k)
    test = c(1:10000)
    for (i in 1:10000){
      xx = rnorm(length(x))
      test[i] = tm(xx,k)}
    cri0.01<-quantile(test,0.01)
    cri0.05<-quantile(test,0.05)
    cri0.1<-quantile(test,0.1)
    ifelse(alpha==0.01,TMcri<-cri0.01,ifelse(alpha==0.05,TMcri<-cri0.05,TMcri<-cri0.1))
    
    for (i in 1:10000) {
      set.seed(i)
      demo.data <- c(rnorm(n,0,1))
      
      
      y<-sort(demo.data)
      
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
      
      k<-1
      x <-c(demo.data)
      tm = function(x,k){
        n = length(x)
        r = abs(x - mean(x))
        df = data.frame(x,r)
        dfs = df[order(df$r),]
        klarge = c((n-k+1):n)
        subx = dfs$x[-klarge]
        ksub = (subx - mean(subx))**2
        all = (df$x - mean(df$x))**2
        ek = sum(ksub)/sum(all)
      }
      ekstat = tm(x,k)
      
      
      ifelse(Q>Qcri,errorQ<-errorQ+1,errorQ<-errorQ+0)
      ifelse(G>Gcri,errorG<-errorG+1,errorG<-errorG+0)
      ifelse(Fer>Fcri,errorFer<-errorFer+1,errorFer<-errorFer+0)
      ifelse(Wor>Zcri,errorWor<-errorWor+1,errorWor<-errorWor+0)
      ifelse(ekstat<TMcri,errorTM<-errorTM+1,errorTM<-errorTM+0)
   
    }
    errordf = data.frame(n,alpha,errorQ,errorG,errorFer,errorWor,errorTM)
    names(errordf)=c("n","alpha", "Dixon","Grubbs","Ferguson","Worrapan","Tietjen")
    print(errordf)
  }
}