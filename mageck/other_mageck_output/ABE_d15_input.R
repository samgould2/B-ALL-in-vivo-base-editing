pdf(file='ABE_d15_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_d15_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Vhl","Ctnnb1","Rad21","Smarca4","Arid1b","Atr","Ezh2","Dis3","Sh2b3","Max")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='d15_rep1,d15_rep2,d15_rep3_vs_input_rep1,input_rep2,input_rep3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(68.69105274121348,106.1612857638743,42.98003229851769,67.01567573889473,16.59784970000794,6.8596297246542495),c(61.46041561055943,290.4081453540694,129.81724041184935,152.60196041748318,91.84143500671061,103.65662695033089),c(305.4944187701337,217.5867675160399,140.3429626074047,184.09125383696383,164.87197368674555,132.61950800998216),c(164.49699472237967,409.7299210886719,292.08879092666103,394.82729441348823,328.63742406015723,409.29124023770356),c(2762.1033839098473,1115.1321835197045,1820.949939831076,557.9256859707983,722.5597236070124,843.7344561324727),c(797.1777436546091,569.4105327335077,965.7350114422036,617.6746016898129,497.93549100023824,511.42350502700015),c(947.2134641156807,495.7117888974296,569.2661420762853,198.62477387672413,292.12215472013975,279.72045654978996),c(251.26464029022827,115.81231174240833,144.7286801888861,286.63331189527264,219.09161604010484,233.2274106382445),c(216.91911391962154,314.0970273013802,128.06295337925678,157.44646709740329,194.74810314675986,194.3561755318704),c(200.6501803756499,150.0295856663017,54.382898010369324,97.69755137838871,7.745663193337039,29.725062140168415),c(0.0,0.0,4.385717581481397,0.0,0.0,0.0),c(2284.8813332866803,958.961035867063,1891.9985646510747,905.9227491450589,1393.1128514873333,1102.8760235082998),c(2395.1485495291545,1127.415307492384,1820.0727963147797,1007.657389423381,1067.7949973671775,1095.2542127031286),c(151.84337974373506,301.8139033287005,184.20013842221869,280.17396965537915,395.028822860189,259.9037484563443))
targetgene="Vhl"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1133.4023702300224,899.3001479997616,1242.0352190755316,994.7387049435941,973.7405157337993,1281.2263963493103),c(280.18718881284445,544.8442847881483,553.4775587829523,465.88005905231637,564.32688980027,592.2146995618169),c(1362.9750991282885,747.5158303373629,790.3063081829478,448.92428567259606,468.05936154022396,491.60679693355456),c(312.7250559007877,579.0615587120417,394.7145823333257,222.84730727632464,156.01978718007464,236.27613496031304),c(1446.1274261308101,901.9322459939074,825.392048834799,786.4249177070297,752.4358530670266,785.8086940131701),c(218.72677320228505,417.62621507110885,257.88019379110614,99.31238693836207,119.50451784005718,54.11485671671686),c(659.7956381721822,543.9669187900997,353.4888370674006,166.3280626772568,205.81333628009847,209.5997971422132),c(97.6136012638297,314.0970273013802,64.03147668962839,89.62337357852186,81.88272518670584,108.22971343343372),c(614.6041561055944,532.5611608154686,435.0631840829546,55.71182681908116,158.23283380674238,285.05572411340995))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1717.276318530337,1051.0844656621605,1317.4695614770117,799.3436021868167,1150.7842458672174,724.8342075717991),c(1122.5564145340413,654.5150345442167,649.0862020592467,331.04128979454026,321.9982841801541,449.6868375051119))
targetgene="Rad21"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(668.8339345854997,570.2878987315562,525.4089662614714,256.75885403576535,169.29806694008101,391.76107538580936),c(526.0288512550823,528.1743308252259,695.5748084229496,556.310850410825,492.4028744335689,377.2796348559837),c(435.6458871219066,793.1388622358874,654.3490631570245,845.3664156460576,812.1881119870553,775.1381588859302),c(735.7173280440496,789.6293982436932,925.3864096925747,1073.0582296023024,741.3706199336881,710.3527670419734),c(469.9914134925133,509.7496448662064,850.8292108073911,505.44353027166386,390.6027296068536,341.4571240716782),c(723.0637130654051,487.81549491499266,335.06882322517873,259.1811073757254,182.57634670008736,220.27033226945312),c(1272.592134995113,829.9882341539264,735.9234101725784,495.75451691182366,578.7116928736102,852.8806290986784))
targetgene="Smarca4"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(224.14975105027557,268.4739954028557,81.57434701555398,185.70608939693722,263.35254857345933,229.4165052356588),c(726.6790316307321,428.1546070476914,492.954656158509,445.6946145526493,487.9767811802335,416.15086996235783),c(1370.2057362589426,896.668050005616,1355.1867326777517,863.9370245857514,897.3904071137628,978.6405073840062),c(1017.7121761395576,736.9874383607803,871.0035116822055,771.0839798872828,837.6381481937341,973.3052398203863),c(488.0680063191484,478.16446893645866,406.99459156147367,621.7116905897464,479.12459467356257,394.8097997078779),c(1317.7836170617009,746.6384643393143,1042.0464973599799,935.7972070045662,1004.7231685071474,1265.9827747389677),c(97.6136012638297,33.33990792584483,23.682874939999543,96.08271581841534,4.426093253335451,21.341070254479888),c(260.3029367035458,110.54811575411705,173.6744162266633,29.87445785950729,39.83483928001906,62.49884860240539))
targetgene="Arid1b"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1193.0551265579184,780.8557382632077,1120.989413826645,461.0355523723963,719.2401536670108,654.7135481642223),c(1990.2328702125276,853.6771161012372,1398.1667649762694,646.7416417693335,792.2706923470457,975.5917830619377),c(441.0688649698971,525.5422328310802,654.3490631570245,568.4221171106252,371.7918332801779,474.07663208166036),c(81.34466771985807,402.71099310428355,475.41178583258346,584.5704727103589,234.5829424267789,466.454821276489))
targetgene="Atr"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(889.3683670704482,1010.7256297519273,1270.980955113309,744.4391931477222,1317.8692661806306,867.3620696285041),c(1560.0099609386116,801.9125222163728,902.5806782688715,749.2836998276423,658.3813714336484,588.4037941592312),c(0.0,9.651025978534028,20.174300874814428,0.0,37.621792653351335,2.2865432415514166),c(538.6824662337268,618.5430286242263,660.4890677710983,702.4534685884147,689.3640242069965,594.5012428033683))
targetgene="Ezh2"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(61.46041561055943,142.13329168386477,87.71435162962794,20.18544449966709,27.663082833346568,23.627613496031305),c(668.8339345854997,644.8640085656828,468.3946377022132,596.6817394101591,367.36574002684245,291.91535383806416),c(2263.1894218947177,1386.238276916706,1295.5409735696046,1226.4676077997722,994.7644586871427,909.2820290569466),c(1522.0491160026777,930.0079579314609,668.383359417765,971.3235893239803,958.2491893471251,1055.6207965162373),c(124.72849050378238,485.18339692084703,174.5515597429596,415.820156693142,285.4830148401366,448.1624753440776))
targetgene="Dis3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(477.2220506231674,460.6171489754877,242.0916104977731,212.35087613649776,331.95699400015883,208.83761606169605))
targetgene="Sh2b3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2104.115405020329,902.8096119919559,1251.6837977547907,644.3193884293735,585.3508327536134,787.3330561742044))
targetgene="Max"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Rbm10","Erbb3","Myd88","Pik3r3","Syk","Tcf3","Daxx","Map2k1","Keap1","Smo")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='d15_rep1,d15_rep2,d15_rep3_vs_input_rep1,input_rep2,input_rep3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(180.76592826635127,314.0970273013802,121.04580524888655,352.8415698541807,483.55068792689804,442.8272077804577),c(173.5352911356972,590.4673166866728,311.3859482851792,587.8001438303056,673.8726978203224,539.6242050061343),c(155.4586983090621,250.92667544188473,156.13154590073773,308.4335919549131,308.7200044201477,381.09054025856943),c(83.15232700252159,286.89868136187516,71.92576833629491,299.55199637505956,320.8917608668202,236.27613496031304))
targetgene="Rbm10"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(500.72162129779304,670.3076225090906,580.6690077881369,628.9784506096264,1201.6843182805749,1066.2913316434772),c(714.0254166520875,774.7141762768679,641.1919104125802,663.6974151490539,847.5968580137389,835.3504642467842),c(1193.0551265579184,667.675524514945,871.8806551985017,623.3265261497197,971.5274691071315,923.7634695867723),c(303.6867594874701,278.12502138138973,170.1658421614782,559.5405215307717,515.6398640135801,544.1972914892372),c(57.84509704523241,315.85175929747726,128.94009689555307,654.8158195692004,440.3962787068774,594.5012428033683),c(242.2263438769107,473.7776389462159,264.89734192147637,628.1710328296398,717.027107040343,600.5986914475054),c(128.3438090691094,551.8632127725367,232.44303181851404,523.206721431371,641.7835217336404,490.8446158530374),c(1708.2380221170195,1041.4334396836264,785.0434470851701,835.6774022862174,834.3185782537325,862.7889831454012),c(1131.5947109473589,944.0458139002376,642.0690539288765,447.30945011262264,559.9007965469345,659.2866346473251),c(140.99742404775398,492.2023249052354,399.1002999148071,503.8286947116905,484.6572112402319,357.46292676253813),c(1189.4398079925913,709.7890924212753,855.2149283888724,641.8971350894134,628.505241973634,730.169475135419),c(421.18461286059846,554.4953107666823,389.45172123554806,512.7102902915441,797.8033089137151,842.2100939714384),c(507.95225842844707,475.532370942313,492.07751264221275,700.0312152484546,595.3095425736182,388.7123510637408),c(260.3029367035458,485.18339692084703,511.3746700007309,436.0056011928091,642.8900450469743,814.0093939923042),c(674.2569124334902,489.5702269110898,270.1602030192541,528.8586458912777,549.9420867269298,404.71815375460073),c(694.1411645427888,828.2335021578293,427.16889243628805,1006.0425538634076,739.1575733070204,955.7750749684922),c(115.69019409046481,80.71767182046642,41.225745265925134,160.67613821735003,179.25677676008576,169.96638095532197),c(309.1097373354607,330.7669812643026,140.3429626074047,765.432055427376,746.9032365003574,653.9513670837051))
targetgene="Erbb3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1268.976816429786,754.5347583217513,791.1834516992441,533.7031525711978,591.9899726336166,628.0372103461224),c(258.4952774208823,656.2697665403139,568.3889985599891,891.3892291052986,898.4969304270966,808.6741264286843),c(209.68847678896748,274.61555738919554,90.34578217851679,450.5391212325694,388.38968298018585,364.32255648719234),c(166.30465400504318,384.286307145264,199.9887217155517,385.9456988336347,419.372335753534,380.32835917805227),c(92.19062341583916,321.9933212838171,123.6772357977754,192.16543163683068,272.20473508013026,213.41070254479888),c(274.76421096485393,383.40894114721544,201.74300874814426,478.79874353210334,478.0180713602287,475.60099424269464))
targetgene="Myd88"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(164.49699472237967,248.29457744773907,88.59149514592421,284.2110585553126,407.2005793068615,333.07313218598966),c(43.3838227839243,146.52012167410751,38.59431471703629,242.22533399600505,268.88516514012866,291.15317275754705),c(206.07315822364046,226.36042749652538,66.66290723851723,180.86158271701711,88.52186506670903,86.12646209843669))
targetgene="Pik3r3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(124.72849050378238,78.08557382632077,91.22292569481306,462.65038793236965,505.6811541935753,592.2146995618169))
targetgene="Syk"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(262.1105959862093,464.12661296768187,288.58021686147595,701.646050808428,478.0180713602287,470.26572667907465),c(186.18890611434182,573.7973627237503,309.63166125258664,714.5647352882149,896.2838838004288,544.1972914892372),c(216.91911391962154,493.9570569013325,203.49729578073683,478.79874353210334,535.5572836535896,547.2460158113057))
targetgene="Tcf3"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(620.0271339535849,708.9117264232268,624.5261836029509,687.1125307686676,1079.9667538138501,629.5615725071567),c(117.49785337312832,186.00159158629216,84.20577756444283,681.4606063087609,599.7356358269536,785.8086940131701))
targetgene="Daxx"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1041.2117468141832,771.2047122846736,940.2978494696115,581.3408015904122,600.8421591402874,603.647415769574),c(1064.711317488809,611.5241006398379,555.2318458155448,661.2751618090938,559.9007965469345,587.641613078714),c(535.0671476683998,448.334025002808,258.75733730740245,675.8086818488541,596.4160658869521,573.1601725488885),c(375.99313079401065,543.9669187900997,441.20318869702857,547.4292548309714,284.37649152680274,192.069632290319),c(197.0348618103229,293.9176093462636,159.64011996592285,532.0883170112245,530.0246670869203,479.41189964528036),c(160.88167615705262,125.46333772094236,110.5200830533312,171.9799871371636,356.30050689350384,333.83531326650683),c(332.60930801008635,419.3809470672059,176.30584677555217,304.3965030549797,574.2855996202748,578.4954401125084),c(209.68847678896748,194.77525156677765,129.81724041184935,88.81595579853519,39.83483928001906,54.877037797233996),c(101.22891982915671,291.2855113521179,123.6772357977754,170.36515157719023,231.2633724867773,225.6055998330731),c(207.88081750630397,245.66247945359345,140.3429626074047,502.21385915171714,551.0486100402636,495.4177023361402))
targetgene="Map2k1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1202.093422971236,617.6656626261778,655.2262066733207,438.4278545327692,391.7092529201874,650.9026427616366),c(1646.7776065064602,690.4870404642072,992.049316931092,638.6674639694667,705.9618739070045,682.1520670628393),c(211.49613607163099,253.55877343603038,141.220106123701,582.9556371503855,340.80918050682976,426.82140508959776),c(99.4212605464932,268.4739954028557,94.73149975999817,426.3165878329689,307.6134811068139,450.44901858562906),c(121.11317193845535,302.6912693267491,241.21446698147685,242.22533399600505,225.730755920108,182.9234593241133),c(157.2663575917256,387.7957711374582,165.7801245799968,430.35367673290233,476.91154804689484,701.2065940757677))
targetgene="Keap1"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(189.80422467966883,343.9274712350308,171.04298567777448,443.27236121268925,556.581226606933,458.07082939080044),c(92.19062341583916,201.79417955116602,70.17148130370235,285.82589411528596,292.12215472013975,319.3538727366812),c(484.4526877538214,623.8072246125176,527.1632532940639,685.4976952086943,614.1204389002938,799.5279534624786))
targetgene="Smo"
collabel=c("input_rep1","input_rep2","input_rep3","d15_rep1","d15_rep2","d15_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("ABE_d15_input_summary.Rnw");
library(tools);

texi2dvi("ABE_d15_input_summary.tex",pdf=TRUE);

