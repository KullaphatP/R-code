#### Detect k extreme outliers ####


for (alpha in c(0.01,0.05,0.1)) {
  for (n in c(10,20,30,50,100)) {
    for (k in c((n/10),(n/5))) {
      
      powerTM<-0
      powerW<-0
      # critical value for Tietjen Moore
      set.seed(1234)
      demo.data <- sort(c(rnorm(n-k,0,1)))
      Q3<-demo.data[floor((3*n)/4)]+(((3*(n))/4)-floor((3*(n))/4))*(demo.data[ceiling((3*(n))/4)]-demo.data[floor((3*(n))/4)])
      Q1<-demo.data[floor(n/4)]+((n/4)-floor(n/4))*(demo.data[ceiling(n/4)]-demo.data[floor(n/4)])
      IQR<-Q3-Q1
      Ufence4.5<-Q3+(4.5*IQR)
      Ufence3<-Q3+(3*IQR)
      outlier<-sample((seq(Ufence3,Ufence4.5,0.001)),k)
      x <-c(demo.data,outlier)
      x<-sort(x)
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
      test = c(1:10000)
      for (i in 1:10000){
        xx = rnorm(length(x))
        test[i] = tm(xx,k)}
      cri0.01<-quantile(test,0.01)
      cri0.05<-quantile(test,0.05)
      cri0.1<-quantile(test,0.1)
      ifelse(alpha==0.01,TMcri<-cri0.01,ifelse(alpha==0.05,TMcri<-cri0.05,TMcri<-cri0.1))
      ##########################
      
      zcri = qnorm(alpha,lower.tail = F) # critical value for T_w
      
      for (i in 1:10000) {
        set.seed(i)
        demo.data <- sort(c(rnorm(n-k,0,1)))
        Q3<-demo.data[floor((3*n)/4)]+(((3*(n))/4)-floor((3*(n))/4))*(demo.data[ceiling((3*(n))/4)]-demo.data[floor((3*(n))/4)])
        Q1<-demo.data[floor(n/4)]+((n/4)-floor(n/4))*(demo.data[ceiling(n/4)]-demo.data[floor(n/4)])
        IQR<-Q3-Q1
        Ufence4.5<-Q3+(4.5*IQR)
        Ufence3<-Q3+(3*IQR)
        outlier<-sample((seq(Ufence3,Ufence4.5,0.001)),k)
        
        
        x <-c(demo.data,outlier)
        x<-sort(x)
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
        
        
        #Worrapan
        
        r = abs(x - mean(x))
        df = data.frame(x,r)
        dfs = df[order(df$r),]
        out<-sort(((n-k+2):n),decreasing = T)
        subx = dfs$x[-out]
        d<-length(subx)
        ksub = ((subx - mean(subx))^3)/(d*((sd(subx))^3))
        ddown = sqrt((((d-1)/d)^3)*((6*(d-2))/((d+1)*(d+3))))
        Wor<- sum(ksub)/ddown
        ifelse(ekstat<TMcri,powerTM<-powerTM+1,powerTM<-powerTM+0)
        ifelse(Wor>zcri,powerW<-powerW+1,powerW<-powerW+0)
      }

      
      powerdf = data.frame(n,k,alpha,powerTM,powerW)
      names(powerdf)=c("n","k","alpha","powerTM","powerW")
      print(powerdf)
    }
  }
}
