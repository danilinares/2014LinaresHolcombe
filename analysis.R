  library('plyr')
library('ggplot2')
library('gridExtra')
library('PropCIs')
library('png')
library('reshape2')
library('boot')
library('utils')

### setttings
tG<-function(t) textGrob(t,x=0,y=1,hjust=0,vjust=1,gp=gpar(fontsize=14))
theme_dur<-theme_set(theme_bw(8))
theme_dur<-theme_update(strip.background=element_blank())
ylimits<-c(-.25,.35)
sl<-.25
dpl<-expression(Delta~'L'['VA'] ~ '(s)')
singlecolumn<-8.9*0.393701
doublecolumn<-18.3*0.393701
defaultheight<-7

### functions
read.files<-function(path='data',subject='.',exp='.',session='.') {
  subject<-gsub(',',replace='|',subject)
  exp<-gsub(',',replace='|',exp)
  session<-gsub(',',replace='|',session)
  r.exp<-paste('^(',subject,').*','(',exp,').*','[',session,']',sep='')
  print(path)
  names.files<-list.files(pattern=r.exp,path=path,full.names=TRUE)
  print(names.files)
  readfun<-function(x) {
    y<-gsub( paste('(',path,'/)|(\\.txt)',sep="") ,'',x)
    subj<-substr(y,0,2)
    session<-gsub('[a-zA-Z]','',y)
    exp<-gsub( paste('(',session,')|(',subj,')',sep='') ,'',y)
    df<-read.table(x, header = TRUE)
    df$subject<-toupper(subj); df$exp<-exp; df$session<-session
    return(df)
  }
  data<- do.call('rbind', lapply(names.files, readfun))
  data
}
averages<-function(df) {
  n<-length(df$response)
  nyes<-sum(df$response)
  nno<-n-nyes
  y<-nyes/n
  ci<-exactci(nyes,n,.95)$conf.int
  ymin<-ci[[1]]; ymax<-ci[[2]]
  data.frame(nyes,nno,n,y,ymin,ymax)
}
lnorm<-function (p, d) {
  pr<-pnorm(d$x,p[1],p[2])
  -sum(d$nyes*log(pr)+d$nno*log(1-pr))
}
curveLnorm<-function(pIni,d){
  p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d)$par
  xseq<-seq(min(d$x),max(d$x),len=100)
  yseq<-pnorm(xseq,p[1],p[2])
  data.frame(x=xseq,y=yseq)
}
fitLnorm<-function(pIni,ythre,d){
  p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d)$par
  xthre<-qnorm(ythre,p[1],p[2])

  yPred<-pnorm(d$x,p[1],p[2])
  sampling<-function(f){
    nyes<-rbinom(length(d$x),d$n,yPred)
    nno<-d$n-nyes
    data.frame(x=d$x,nyes,nno)
  }
  samples<-ddply(data.frame(sample=1:1000),.(sample),sampling)
  calculateThre<-function(f){
    p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=f)$par
    xthre<-qnorm(ythre,p[1],p[2])
    data.frame(x=xthre)
  }
  samplesThre<-ddply(samples,.(sample),calculateThre)
  quant<-quantile(samplesThre$x,c(.025,.975))
  error<-sum((d$y-yPred)^2)
  data.frame(x=xthre,xmin=quant[[1]],xmax=quant[[2]],error)
}
difLnorm<-function(pIni,ythre,d,condition){
  thre<-ddply(d,condition,function(d2){
    p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d2)$par
    xthre<-qnorm(ythre,p[1],p[2])
    data.frame(x=xthre)
  })
  x=-diff(thre$x)

  ci<-ddply(d,condition,function(d2){
    p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d2)$par
    yPred<-pnorm(d2$x,p[1],p[2])
    sampling<-function(f){
      nyes<-rbinom(length(d2$x),d2$n,yPred)
      nno<-d2$n-nyes
      data.frame(x=d2$x,nyes,nno)
    }
    samples<-ddply(data.frame(sample=1:1000),.(sample),sampling)

    calculateThre<-function(f){
      p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=f)$par
      xthre<-qnorm(ythre,p[1],p[2])
      data.frame(x=xthre)
    }
    samplesThre<-ddply(samples,.(sample),calculateThre)
    data.frame(x=samplesThre$x)
  })
  condi<-unique(d[[condition]])
  diffe<-ci$x[ci[[condition]]==condi[1]]-ci$x[ci[[condition]]==condi[2]]
  quant<-quantile(diffe,c(.025,.975))
  data.frame(x,xmin=quant[[1]],xmax=quant[[2]])
}
difLnormSta<-function(pIni,ythre,d,condition){
  thre<-ddply(d,condition,function(d2){
    p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d2)$par
    xthre<-qnorm(ythre,p[1],p[2])
    if (unique(d2[[condition]])==0) xthreSTa<- xthre-unique(d$sta)
    if (unique(d2[[condition]])==1) xthreSTa<- -xthre+unique(d$sta)
    data.frame(x=xthreSTa)
  })
  x=-diff(thre$x)

  ci<-ddply(d,condition,function(d2){
    p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=d2)$par
    yPred<-pnorm(d2$x,p[1],p[2])
    sampling<-function(f){
      nyes<-rbinom(length(d2$x),d2$n,yPred)
      nno<-d2$n-nyes
      data.frame(x=d2$x,nyes,nno)
    }
    samples<-ddply(data.frame(sample=1:1000),.(sample),sampling)

    calculateThre<-function(f){
      p<-optim(par=c(pIni[1],pIni[2]),lnorm,d=f)$par
      xthre<-qnorm(ythre,p[1],p[2])
      if (unique(d2[[condition]])==0) xthreSTa<- xthre-unique(d$sta)
      if (unique(d2[[condition]])==1) xthreSTa<- -xthre+unique(d$sta)
      data.frame(x=xthreSTa)
    }
    samplesThre<-ddply(samples,.(sample),calculateThre)
    data.frame(x=samplesThre$x)
  })
  condi<-unique(d[[condition]])
  diffe<-ci$x[ci[[condition]]==condi[1]]-ci$x[ci[[condition]]==condi[2]]
  quant<-quantile(diffe,c(.025,.975))
  data.frame(x,xmin=quant[[1]],xmax=quant[[2]])
}
lgauss<-function (p, d) {
  pr<-exp(-0.5*((d$x-p[2])/p[3])^2)/(1+exp(-p[1]))
  -sum(d$nyes*log(pr)+d$nno*log(1-pr))
}
curveLgauss<-function(pIni,d){
  p<-optim(par=c(pIni[1],pIni[2],pIni[3]),lgauss,d=d)$par

  xseq<-seq(min(d$x),max(d$x),len=100)
  yseq<-exp(-0.5*((xseq-p[2])/p[3])^2)/(1+exp(-p[1]))
  data.frame(x=xseq,y=yseq)
}
fitLgauss<-function(pIni,d){
  p<-optim(par=c(pIni[1],pIni[2],pIni[3]),lgauss,d=d)$par

  yPred<-exp(-0.5*((d$x-p[2])/p[3])^2)/(1+exp(-p[1]))
  sampling<-function(f){
    nyes<-rbinom(length(d$x),d$n,yPred)
    nno<-d$n-nyes
    data.frame(x=d$x,nyes,nno)
  }
  samples<-ddply(data.frame(sample=1:1000),.(sample),sampling)
  calculateThre<-function(f){
    p<-optim(par=c(pIni[1],pIni[2],pIni[3]),lgauss,d=f)$par
    data.frame(x=p[2])
  }
  samplesThre<-ddply(samples,.(sample),calculateThre)
  quant<-quantile(samplesThre$x,c(.025,.975))

  data.frame(x=p[2],xmin=quant[[1]],xmax=quant[[2]])
}
lnormFree<-function (p, d) {
  pr<-p[3]+(1-p[3]-p[4])*pnorm(d$x,p[1],p[2])
  -sum(d$nyes*log(pr)+d$nno*log(1-pr))
}

curveLnormFree<-function(pIni,d){
  p<-optim(par=c(pIni[1],pIni[2],pIni[3],pIni[4]),lnormFree,d=d)$par
  xseq<-seq(min(d$x),max(d$x),len=100)
  yseq<-p[3]+(1-p[3]-p[4])*pnorm(xseq,p[1],p[2])
  data.frame(x=xseq,y=yseq)
}

fitLnormFree<-function(pIni,ythre,d){
  p<-optim(par=c(pIni[1],pIni[2],pIni[3],pIni[4]),lnormFree,d=d)$par
  q<-(ythre-p[3])/(1-p[3]-p[4])
  xthre<-qnorm(q,p[1],p[2])

  yPred<-p[3]+(1-p[3]-p[4])*pnorm(d$x,p[1],p[2])

  sampling<-function(f){
    nyes<-rbinom(length(d$x),d$n,yPred)
    nno<-d$n-nyes
    data.frame(x=d$x,nyes,nno)
  }
  samples<-ddply(data.frame(sample=1:1000),.(sample),sampling)

  calculateThre<-function(f){
    p<-optim(par=c(pIni[1],pIni[2],pIni[3],pIni[4]),lnormFree,d=f)$par
    q<-(ythre-p[3])/(1-p[3]-p[4])
    xthre<-qnorm(q,p[1],p[2])
    data.frame(x=xthre)
  }
  samplesThre<-ddply(samples,.(sample),calculateThre)
  quant<-quantile(samplesThre$x,c(.025,.975),na.rm = T)
  error<-sum((d$y-yPred)^2)
  data.frame(x=xthre,xmin=quant[[1]],xmax=quant[[2]])
}
################################################################################
#### Duration ##################################################################
################################################################################
### Figure b
datdur<-read.files(subject='',exp='Dur')
datdur$condition<-paste(datdur$sta,datdur$soundFirst)
datdur$x<-datdur$soa

datdur$condition<-factor(datdur$condition,
                         levels=c('0.6 0','0.6 1','1.2 0','1.2 1'),
                         labels=c('VA (0.6 s)  AV (variable)','AV (0.6 s)  VA (variable)',
                                  'VA (1.2 s)  AV (variable)','AV (1.2 s)  VA (variable)'))

datdur$soundFirst<-factor(datdur$soundFirst,levels=c(1,0))

ntrials<-ddply(datdur,.(subject,sta),function(d) length(d$response))
avdur<-ddply(datdur,.(subject,x,condition),function(d) averages(d))
curvedur<-ddply(avdur,.(subject,condition),function(d) curveLnorm(c(1,1),d))

pdur<-ggplot() +
  facet_wrap(~subject,ncol=2,scales='free_x')+
  geom_vline(data=avdur,aes(xintercept=.6),size=sl,lty=2)+
  geom_vline(data=avdur,aes(xintercept=1.2),size=sl,lty=2)+
  geom_hline(data=avdur,aes(yintercept=.5),size=sl,lty=2)+
  geom_point(data=avdur,aes(x=x,y=y,color=condition,shape=condition),size=3)+
  geom_line(data=curvedur,aes(x=x,y=y,color=condition))+
  scale_x_continuous(breaks=c(.6,1.2,1.8,2.4))+
  scale_y_continuous(breaks=c(0,.5,1)) +
  scale_shape_manual(values=c(16,16,17,17))+
  scale_color_manual(values=c('#D55E00','#0072B2','#D55E00','#0072B2'))+
  xlab('Duration variable (s)')+
  ylab('Proportion variable reported longer')+
  theme(legend.position=c(x=.8,y=0.1),
        legend.direction='vertical',
        legend.key.height=unit(.5,'line'),
        axis.title.x=element_text(hjust=0.12),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=sl),
        legend.title=element_blank())
pdur

### Figure c
avdur2<-ddply(datdur,.(subject,x,sta,soundFirst),function(d) averages(d))
thredur<-ddply(avdur2,.(subject,sta,soundFirst),
               function(d) fitLnorm(c(1,1),.5,d))
pdurbars<-ggplot()+
  facet_grid(.~subject)+
  geom_bar(data=thredur,aes(x=factor(sta),y=x,fill=soundFirst),
           position='dodge',stat="identity")+
  geom_errorbar(data=thredur,size=sl,width=0.5,
                aes(x=factor(sta),ymin=xmin,ymax=xmax,fill=soundFirst),
                position=position_dodge(width=.9),stat='identity')+
  scale_fill_manual(values=c('#0072B2','#D55E00'),
                    labels=c('AV (standard) VA (variable)','VA (standard) AV (variable)'))+
  scale_y_continuous(breaks=c(0,.6,1.2,1.8))+
  xlab('Duration standard (s)')+
  ylab('Duration variable that\nmatches standard (s)')+
  theme(legend.position='top',
        legend.direction='horizontal',
        legend.key.height=unit(.25,'line'),
        legend.title=element_blank())
pdurbars

### Figure d
difdur<-ddply(avdur2,.(subject,sta),function(d)
  difLnorm(c(1,1),.5,d,'soundFirst'))
difdur4<-difdur
difdur4$x<-.25*difdur$x #dividing by 4
difdur4$xmin<-.25*difdur$xmin
difdur4$xmax<-.25*difdur$xmax

pdurbarsdif<-ggplot()+
  facet_grid(.~subject)+
  geom_bar(data=difdur4,aes(x=factor(sta),y=x),fill='grey',
           position='dodge',stat="identity")+
  geom_errorbar(data=difdur4,aes(x=factor(sta),ymin=xmin,ymax=xmax),width=0.5,
                position=position_dodge(width=.9),stat='identity',size=sl)+
  ylim(ylimits)+
  xlab('Duration standard (s)')+
  ylab(dpl)+
  theme(legend.position='top',
        legend.direction='horizontal',
        legend.title=element_blank())
pdurbarsdif

### Panel
image<-rasterGrob(readPNG('durationLow.png'))
p.a<-arrangeGrob(image,left=tG('a'))
p.b<-arrangeGrob(pdur,left=tG('b'))
p.c<-arrangeGrob(pdurbars,left=tG('c'))
p.d<-arrangeGrob(pdurbarsdif,left=tG('d'))
p.acd<-arrangeGrob(p.a,p.c,p.d,heights=c(0.36,.34,.30))
pdf('figures/durationFig.pdf',width=doublecolumn,height=7)
p.ab<-grid.arrange(p.acd,p.b,ncol=2,widths=c(.55,.45))
dev.off()

### statistics

# Checking asymmetries in the pse for individual observers
difdurSta<-ddply(avdur2,.(subject,sta),function(d)
  difLnormSta(c(1,1),.5,d,'soundFirst'))
difdurSta$symmetric<-sign(difdurSta$xmin)*sign(difdurSta$xmax)
length(difdurSta$symmetric[difdurSta$symmetric==1])

# Checking asymmetries in the pse across observers
difdurSta$xAbs<-abs(difdurSta$x)
ddply(difdurSta,.(sta),function(d) {
  print(d$sta)
  print(t.test(d$xAbs))
  data.frame(x=1)
})

# correlation long short
difdurLong<-dcast(difdur,subject~sta,value.var='x')
plot(difdurLong[['0.6']],difdurLong[['1.2']])
cor.test(difdurLong[['0.6']],difdurLong[['1.2']])

funCor<-function(d,i) {
  d2<-d[i,]
  cor(d2[,1],d2[,2])
}
bootCorDur<-boot( difdurLong[,c('0.6','1.2')],funCor,R=1000)
quantile(bootCorDur$t,c(.025,.975),na.rm=T)

# the absolute value of the difference between incrementLVA for each standard
difdursta<-ddply(difdur,.(subject),function(d) data.frame(dif=d$x[1]-d$x[2]))
difdursta$difAbs<-abs(difdursta$dif)
t.test(difdursta$difAbs)

# t-test for the effect of the standard
t.test(difdur4$x[difdur4$sta==.6],difdur4$x[difdur4$sta==1.2],paired=T)

t.test(difdur4$x[difdur4$sta==.6])
t.test(difdur4$x[difdur4$sta==1.2])

# 2-way within subjects anova
#summary(aov(x~(sta*soundFirst)+Error(subject/(sta*soundFirst)),thredur))

#anova
#summary(aov(x~sta+Error(subject/sta),difdur4))



################################################################################
#### SJ ########################################################################
################################################################################
### Figure a
datsj<-read.files(subject='',exp='Synch')
datsj$soa<- -datsj$soa #because the way the data was coded
datsj$x<-datsj$soa

ntrialssj<-ddply(datsj,.(subject),function(d) length(d$response))

avsj<-ddply(datsj,.(subject,x),function(d) averages(d))



curvesj<-ddply(avsj,.(subject),function(d) curveLgauss(c(1,0,1),d))

meansj<-ddply(avsj,.(subject),function(d) fitLgauss(c(1,0,1),d))

pfitsj<-ggplot() +
  facet_wrap(~subject,ncol=4,scales='free_x')+
  geom_pointrange(data=avsj,size=.5,aes(x=x,y=y,ymin=ymin,ymax=ymax),
                  color='#CC79A7')+
  geom_line(data=curvesj,size=0.5,aes(x=x,y=y),color='#CC79A7')+
  geom_vline(data=meansj,aes(xintercept=x),color='#CC79A7',lty=2)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Time V was presented before A (s)')+
  ylab('Proportion of simultaneity reports')+
  theme(#axis.title.x=element_text(hjust=0),
    panel.margin = unit(.5, 'lines'))
pfitsj

### Figure b
pbarsj<-ggplot()+
  geom_bar(data=meansj,aes(x=subject,y=x),fill='#CC79A7',
           position='dodge',stat="identity")+
  geom_errorbar(data=meansj,aes(x=subject,ymin=xmin,ymax=xmax),fill='grey',
                position=position_dodge(width=.9),stat="identity",width=.5,size=sl)+
  ylab(dpl)+
  ylim(ylimits)+
  theme(axis.title.x=element_blank())
pbarsj

#### Panel
pa<-arrangeGrob(pfitsj,left=tG('a'))
pb<-arrangeGrob(pbarsj,left=tG('b'))
pdf('figures/sjFig.pdf',width=doublecolumn,height=3)
grid.arrange(pa,pb,ncol=2, widths=c(.6,.4))
dev.off()

### Statistics
mean(meansj$x)
t.test(meansj$x,mu=0)

################################################################################
#### TOJ #######################################################################
################################################################################
### Figure a
dattoj<-read.files(subject='',exp='TOJ')
dattoj$soa<- -dattoj$soa #because the way the data was encoded
dattoj$x<-dattoj$soa
dattoj$response<- -dattoj$response+1 #because the way the data was encoded

ntrialstoj<-ddply(dattoj,.(subject),function(x) length(x$response))

avtoj<-ddply(dattoj,.(subject,x),function(d) averages(d))
curvetoj<-ddply(avtoj,.(subject),function(d) curveLnormFree(c(0,.25,.01,.01),d))
meantoj<-ddply(avtoj,.(subject),function(d) fitLnormFree(c(0,.25,.01,.01),.5,d))

### Figure a
pfittoj<-ggplot() +
  facet_wrap(~subject,ncol=4,scales='free_x')+
  geom_pointrange(data=avtoj,size=0.5,aes(x=x,y=y,ymin=ymin,ymax=ymax),
                  color='#009E73')+
  geom_line(data=curvetoj,size=.5,aes(x=x,y=y),color='#009E73')+
  geom_vline(data=meantoj,aes(xintercept=x),color='#009E73',lty=2)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Time V was presented before A (s)')+
  ylab('Proportion V reported before A')
pfittoj

### Figure b
pbartoj<-ggplot()+
  geom_bar(data=meantoj,aes(x=subject,y=x),fill='#009E73',
           position='dodge',stat="identity")+
  geom_errorbar(data=meantoj,aes(x=subject,ymin=xmin,ymax=xmax),fill='grey',
                position=position_dodge(width=.9),stat="identity",width=.5,size=sl)+
  coord_cartesian(ylim=ylimits)+
  ylab(dpl)+
  theme(axis.title.x=element_blank())
pbartoj

#### Panel
patoj<-arrangeGrob(pfittoj,left=tG('a'))
pbtoj<-arrangeGrob(pbartoj,left=tG('b'))
pdf('figures/tojFig.pdf',width=doublecolumn,height=3)
grid.arrange(patoj,pbtoj,ncol=2, widths=c(.6,.4))
dev.off()

### Statistics
mean(meantoj$x)
t.test(meantoj$x,mu=0)

################################################################################
#### ALL #######################################################################
################################################################################
meansj$sta<-'SJ'
meantoj$sta<-'TOJ'

latAll<-rbind(difdur4,meansj,meantoj)

pAll<-ggplot()+
  facet_grid(.~subject)+
  geom_bar(data=latAll,aes(x=sta,y=x,fill=sta),position='dodge',stat="identity")+
  geom_errorbar(data=latAll,width=0.5,size=sl,
                aes(x=sta,ymin=xmin,ymax=xmax,fill=sta),
                position=position_dodge(width=.9),stat='identity')+
  ylab(dpl)+
  scale_fill_manual(values=c('grey','grey','#CC79A7','#009E73'))+
  coord_cartesian(ylim=ylimits)+
  theme(axis.title.x=element_blank(),
        panel.margin = unit(1,'lines'),
        legend.position='None')
pAll

pdf('figures/comparisonFig.pdf',width=doublecolumn,height=3)
pAll
dev.off()

#Correlation
latAllLong<-dcast(latAll,subject~sta,value.var='x')
cor.test(latAllLong[['0.6']],latAllLong[['SJ']])
cor.test(latAllLong[['0.6']],latAllLong[['TOJ']])
cor.test(latAllLong[['1.2']],latAllLong[['SJ']])
cor.test(latAllLong[['1.2']],latAllLong[['TOJ']])
cor.test(latAllLong[['SJ']],latAllLong[['TOJ']])

funCor<-function(d,i) {
  d2<-d[i,]
  cor(d2[,1],d2[,2])
}
bootCor<-boot(latAllLong[,c(0.6,'SJ')],funCor,R=1000)
quantile(bootCor$t,c(.025,.975),na.rm=T)

combinations<-data.frame(t(combn(unique(latAll$sta),2)))
correlations<-ddply(combinations,.(X1,X2),function(d) {
  d2<-latAllLong[,c(as.character(d[1,1]),as.character(d[1,2]))]
  cor<-cor.test(d2[,1],d2[,2])
  r<-cor$estimate
  p<-cor$p.value
  bootCor<-boot(d2,funCor,R=1000)
  ci<-quantile(bootCor$t,c(.025,.975),na.rm=T)
  data.frame(r,min=ci[[1]],max=ci[[2]],p)
})
correlations




