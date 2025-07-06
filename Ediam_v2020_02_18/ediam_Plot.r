#This function plots results from run in Ediam

ediamPlot<-function(SimulData)
{
#create function to use data by region
data.by.region<-function(data,id.vars,measure.vars)
{
	data<-data[,c(id.vars,measure.vars)]
	data<-melt(data, id.vars=id.vars, measure.vars=measure.vars,  variable.name="Variable.Region")
	data<-cbind(data,colsplit(data$Variable.Region,("\\."),c("Variable","Region")))
	data$Variable.Region<-NULL
	new.form<-as.formula(paste(id.vars,"Region~Variable",sep="+"))
	data<-dcast(data,new.form, value.var="value")
	data$Region<-ifelse(data$Region=="N","OECD Countries","Non-OECD Countries")
	return(data)
}


#
theme_paper<- theme( panel.border = element_rect(fill = NA, colour = "black", size = 1),
				             #text=element_text(family="Arial"),
				          #style of title
				          plot.title = element_text(size=11),
				          #style axis titles
				          axis.title.x = element_text(size=10),
				          axis.title.y = element_text(size=10),
  			          #style axis
				          axis.text.x = element_text(size=8),
				          axis.text.y = element_text(size=8),
				          panel.grid.major = element_blank(),
				          panel.grid.major.x = element_blank(),
				          panel.grid.minor = element_blank(),
				          panel.grid.minor.x = element_blank(),
                  panel.background = element_blank())
#
Plot0 <- ggplot(SimulData, aes(times+2014,Delta.Temp))+ geom_line()+geom_line(size=1)+
                            geom_hline(yintercept = 2.0,linetype = 2)+geom_hline(yintercept = 6.0,linetype = 2)+
                            ggtitle("Global Temperature Anomaly With Respect to Historic Mean")+
                            scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                            scale_y_continuous(name="Degrees Celsius",limits=c(0, 6.0))+theme_paper
#
Plot1<-ggplot(data.by.region(SimulData,"times",c("Realcetax.N","Realcetax.S")), aes(times+2014,Realcetax,colour =Region)) + geom_line()+geom_line(size=1)+
                         ggtitle("Carbon Price")+
                         #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
                         scale_color_manual(values =c("#E7B800","#00AFBB"))+
                         scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                         scale_y_continuous(name="USD 2010 per tonne of co2")+theme_paper
#
Plot2<-ggplot(data.by.region(SimulData,"times",c("RealRDsubsidy.N","RealRDsubsidy.S")), aes(times+2014,RealRDsubsidy,colour =Region)) + geom_line()+geom_line(size=1)+
                         ggtitle("Renewable Energy R&D Investments")+
                         #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
                         scale_color_manual(values =c("#E7B800","#00AFBB"))+
                         scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                         scale_y_continuous(name="Billion USD 2010")+theme_paper

#
Plot3<-ggplot(data.by.region(SimulData,"times",c("RealTecsubsidy.N","RealTecsubsidy.S")), aes(times+2014,RealTecsubsidy,colour =Region)) + geom_line()+geom_line(size=1)+
                         ggtitle("Renewable Energy Investments")+
                         #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
                         scale_color_manual(values =c("#E7B800","#00AFBB"))+
                         scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                         scale_y_continuous(name="Billion USD 2010")+theme_paper

#
#Graph Scientists
Plot4<-ggplot(data.by.region(SimulData,"times",c("s_re.N","s_re.S")), aes(times+2014,s_re,colour =Region)) + geom_line()+geom_line(size=1)+
                         ggtitle("Entrepreneurs in Renewable Energy")+
                         #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
                         scale_color_manual(values =c("#E7B800","#00AFBB"))+
                         scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                         scale_y_continuous(name="Percent of Total", breaks=seq(0,1,by=0.20),limits=c(0, 1.0))+theme_paper

#
#Graph Rate of decarbonization
Plot5<-ggplot(data.by.region(SimulData,"times",c("DecarbY.N","DecarbY.S")), aes(times+2014,DecarbY,colour =Region)) + geom_line(size=1)+
                          ggtitle("Decarbonization Rate")+
                          #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
                          scale_color_manual(values =c("#E7B800","#00AFBB"))+
                          scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                          scale_y_continuous(name="Percent of Total", breaks=seq(0,1,by=0.20),limits=c(0, 1.0))+theme_paper


#
#World Population
Plot6<-ggplot(data.by.region(SimulData,"times",c("L.N","L.S")), aes(times+2014,L,fill =Region)) +geom_bar(stat="identity")+
                         ggtitle("Population Growth")+
                         scale_fill_manual(values =c("#E7B800","#00AFBB"))+
                         scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
                         scale_y_continuous(name="Million People")+theme_paper



#Putting all graphs into one page

 Plot<-ggarrange(Plot1,Plot4,
           Plot2,Plot5,
           Plot3,Plot6,
           Plot0,
           nrow=4,ncol = 2,
           common.legend = TRUE,
           align="h")

return(Plot)
}
