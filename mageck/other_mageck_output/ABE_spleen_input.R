pdf(file='ABE_spleen_input.pdf',width=4.5,height=4.5);
gstable=read.table('ABE_spleen_input.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Ctnnb1","non_target","Vhl","Dis3","Epha7","Fbxw7","Pik3cd","Foxa1","Arid1a","Traf7")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='spleen1,spleen2,spleen3,spleen4_vs_input_rep1,input_rep2,input_rep3 neg.'


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
targetmat=list(c(921.6885617223307,726.1453448992164,1037.1622098579965,939.7963143000115,903.5687043936992,1181.5806083002888,1137.6593817661244),c(227.8496444449143,439.9378138365009,462.1817474720309,734.1656677523167,804.1412600878289,516.8638623059096,651.9052829469141),c(1108.3782703965508,603.5861793698853,659.945728165293,575.1232145630839,441.899752470535,367.76851740997415,225.4153007592414),c(254.30960315464628,467.5667586668125,329.60663448877006,44.98170393230824,123.84822010029468,29.819068979187094,58.205828834545905),c(1175.998164876977,728.2706483477019,689.2440956754059,1975.981994169255,919.8492215899821,570.2896942269532,616.9817856461866),c(177.86972243764276,337.21481382636784,215.3430011993298,142.97755892769405,118.0337496730508,342.9192932606516,343.9435340223167),c(536.5491627251207,439.2293793536724,295.18105266438744,266.6772447415417,455.273034453196,262.1593147753532,246.58105669907627),c(79.37987612919595,253.61954485260438,53.469520705956036,427.3261873569283,872.1705640865823,715.6576555004903,2005.4553752993543),c(499.79922007271523,430.01973107690185,363.2997571253999,102.81532327384741,58.72615131516321,38.51629743145,47.622950864628464))
targetgene="Ctnnb1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(468.9292682446946,236.6171172647203,408.7122267660749,440.1781027661592,112.80072628853131,300.67561220680324,164.03460853372027),c(302.81952745582157,213.23877933137965,277.6020321583197,178.3203263030791,35.46826960618768,36.03137501651774,111.12021868413309),c(1396.4978207914103,593.6680966102863,901.6572601237243,486.76629612462136,841.935317864914,550.4103149074951,365.1092899621516),c(743.818839284688,229.53277243643524,291.5187567256233,167.07490032000203,104.66046769038988,188.85410353485162,286.79599298476256),c(241.07962379978028,342.1738552061674,376.4840225049507,94.78287614307808,72.68088034054853,72.06275003303548,224.35701296224966),c(864.3586511845781,534.1596000526919,714.1477080590018,2106.107637687718,1260.5771886264736,991.484043557971,1401.1730432170687),c(435.1193210044815,330.8389034809113,485.62044148012126,398.40937768615873,554.119031716342,244.7648578708274,234.9398909321671),c(1894.8270431580293,674.4296276527357,1024.7104036661985,470.7014018630827,359.3342724036719,214.94578889164032,261.3970858569607),c(379.2594081728251,357.05097934556596,541.2873397493357,88.35691843846261,308.74837968665014,410.0121984638226,338.652095037358),c(145.5297729035259,221.73999312532172,152.35151105258706,141.3710695015402,319.7958734984135,84.48736210769677,146.04371598486063))
targetgene="non_target"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(55.85991283165641,85.72057242224896,35.890500199888294,17.671383687692526,80.23969189596556,27.334146564254837,15.874316954876155),c(49.979922007271526,234.49181381623478,108.4039597874177,69.07904532461623,56.98181018699004,74.54767244796774,135.4608380149432),c(248.4296123302614,175.69175174146895,117.19347004045159,91.56989729077036,133.15137278388488,155.30765093326613,101.5956285112074),c(133.76979125475614,330.8389034809113,243.90890952168985,461.0624653061595,434.9223879578424,282.03869409481126,359.81785097719285),c(2246.156494915026,900.4202276750284,1520.5852737748592,459.4559758800056,423.874894146079,586.4416899240128,620.1566490371619),c(648.2689883884336,459.773979355699,806.4375657158574,244.18639277538762,315.1442971566184,272.09900443508224,659.3132975258563),c(770.27879799442,400.2654827981047,475.36601285158173,163.86192146769432,161.06083083465552,377.70820706970323,276.21311501484513),c(204.32968114737474,93.51335173336251,120.85576597921569,147.79702720615566,104.07902064766549,34.788913809051614,687.8870680446335),c(176.39972473154654,253.61954485260438,106.93904141191207,257.03830818461853,288.97918023402093,130.45842678394354,183.08378887957167),c(163.16974537668057,121.14229656367417,45.412469640674985,1.6064894261538658,0.0,0.0,40.214936285686264),c(0.0,0.0,3.662295938764112,0.0,0.0,0.0,0.0),c(1858.0771005056238,774.3188897315547,1579.9144679828378,1116.5101511769367,1125.1000276716911,777.7807158737968,610.6320588642361),c(1947.7469605774932,910.3383104346275,1519.8528145871064,918.9119517600112,1164.6384265769495,1167.9135350181612,1482.661203585433),c(123.47980731208258,243.70146209300532,153.8164294280927,673.1190695584697,508.1847153411153,634.8976770151919,205.3078326163983))
targetgene="Vhl"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(49.979922007271526,114.76638621821762,73.24591877528223,28.916809670769585,26.746563965321858,33.546452601585486,26.45719492479359),c(543.8991512556019,520.6993448789503,391.13320626000717,106.02830212615514,164.54951309100187,241.03747424842902,221.18214957127444),c(1840.437128032469,1119.326482869036,1081.8422203109187,514.0766163692371,635.5216176977563,291.9783837545403,654.0218585408976),c(1237.7380685330183,750.9405517982141,558.1339010676506,372.70554686769685,409.3387180779693,428.6491165758145,640.264117180005),c(101.42984172063927,391.76426900416266,145.75937836281165,430.53916620923604,402.3613535652766,417.46696570861934,650.8469951499223))
targetgene="Dis3"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1043.698371328317,502.2800483254092,675.3273711081022,369.4925680153891,345.3795433782866,236.0676294185645,683.6539168566665),c(1098.0882864538773,631.9235586830255,646.7614627857422,385.5574622769278,316.3071912420672,699.5056598034306,408.4990896388131),c(1178.9381602891694,697.8079655860763,1050.3464752375473,1755.8929427861754,830.3063770104263,2559.470087380226,714.344262969427),c(1494.987667099857,900.4202276750284,1171.202241216763,411.26129309538965,320.37732054113786,593.8964571688097,486.8123866162021),c(1101.0282818660696,580.2078414365446,700.2309834916982,265.0707553153879,193.04041818449687,188.85410353485162,112.17850648112483))
targetgene="Epha7"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(463.0492774203097,678.6802345497067,773.4769022669805,385.5574622769278,679.1301459020854,431.1340389907468,621.2149368341535),c(263.1295893912236,639.0079035113105,648.9588403490006,652.2347070184695,844.2611060358116,1676.080168871808,486.8123866162021),c(355.7394448752855,716.2272621396174,640.1693300959668,682.758006115393,916.9419863763602,533.0158580029694,661.4298731198398),c(948.1485204320627,558.9548069516895,626.2526055286631,175.10734745077139,242.46341681606987,586.4416899240128,664.6047365108151),c(1073.0983254502414,600.0440069557428,549.3443908146168,77.11149245538556,170.36398351824573,272.09900443508224,295.2622953606965),c(1519.9776281034929,789.9044483537818,983.6926891520404,377.52501514615847,1208.2469547812786,556.6226209448258,1441.387979502755),c(742.3488415785918,648.9259862709096,900.9248009359715,828.9485438953948,739.6006383454218,423.67927174594996,663.5464487138233),c(667.3789585676844,721.8947380022454,656.2834322265288,589.5816193984688,616.9153123305758,470.8927976296629,466.704918473359),c(715.8888828688598,692.1404897234482,572.7830848227071,157.43596376307886,117.45230263032641,63.36552158077258,79.37158477438078),c(567.4191145531414,584.4584483335157,508.3266763004587,165.4684108938482,181.9929243727335,237.31009062603064,185.20036447355514),c(452.75929347763616,534.1596000526919,490.01519660663814,631.3503444784693,209.90238242350412,238.55255183349675,385.2167581049947),c(129.35979813646748,282.66535864857303,161.87348049337373,170.28787917230977,119.77809080122397,85.7298233151629,251.872495684035))
targetgene="Fbxw7"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(887.8786144821177,615.6295655779699,852.5824945442853,355.0341631800043,426.20068231697655,749.2041081020758,318.54462689451486),c(407.1893645886533,427.1859931455878,340.5935223050624,239.366924496926,517.4878680247054,172.70210783779194,86.77959935332298))
targetgene="Pik3cd"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1178.9381602891694,620.5886069577693,798.3805146505764,125.30617524000154,211.0652765089529,213.70332768417418,206.36612041339),c(1302.417967601252,811.157482838637,1026.1753220417042,462.66895473231335,284.9090509349502,585.1992287165467,207.42440821038176),c(723.2388713993408,495.9041379799527,396.2604205742769,1489.2156980446337,521.5579973237762,464.68049159233226,2088.0018234647105),c(429.2393301800966,547.6198552264335,364.03221631315273,403.2288459646203,462.2503989658886,673.4139744466419,429.66484557864794),c(1239.2080662391145,670.1790207557647,706.8231161814736,330.9368217876964,398.8726713089303,231.0977845887,324.8943536764653),c(270.47957792170473,306.04369658191365,157.47872536685682,515.6831057953909,179.08568915911155,228.61286217376772,286.79599298476256),c(107.30983254502415,308.16900003039916,106.93904141191207,263.464265889234,214.55395876529923,257.1894699454887,304.7868855336222),c(29.399954121924424,121.85073104650266,58.59673502022579,114.06074925692447,204.08791199626026,217.43071130657256,83.60473596234775),c(2171.1866119041188,753.774289729528,865.0343007360832,388.7704411292355,596.5646658352223,751.6890305170081,354.52641199223416),c(864.3586511845781,830.2852138750065,422.6289513333785,326.11735350923476,318.05153237024035,397.5875863891613,208.4826960073735),c(1089.2683002173,723.3116069679024,665.8054016673156,689.1839638200084,307.00403855847696,777.7807158737968,620.1566490371619),c(2097.686726599308,949.3022069901951,1195.373394412606,1032.9727010169356,501.2073508284226,708.2028882556936,969.3916220444372),c(1196.578132762324,759.4417655921561,795.4506778995651,404.83533539077416,665.1754168767001,705.7179658407613,383.1001825110112))
targetgene="Foxa1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(961.3784997869287,730.3959517961874,595.4893196430446,719.7072629169319,1173.9415792605398,710.6878106706258,941.8761393226519),c(352.7994494630931,639.0079035113105,604.2788298960785,380.7379939984662,195.36620635539444,638.6250606375903,329.1275048644323),c(458.63928430202105,325.17142761828325,277.6020321583197,33.73627794923118,120.94098488667274,64.6079827882387,16.9326047518679))
targetgene="Arid1a"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(145.5297729035259,58.09162759193732,71.78100039977659,14.458404835384792,11.047493811763376,44.72860346878064,17.99089254885964),c(2047.7068045920362,804.0731380103518,1445.8744366240714,909.273015203088,1055.907829587489,654.77705633465,479.4043720372599))
targetgene="Traf7"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetgenelist=c("Ercc2","Myd88","Smo","Rbm10","Tcf3","Ret","Syk","Keap1","Foxo1","Arid5b")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='spleen1,spleen2,spleen3,spleen4_vs_input_rep1,input_rep2,input_rep3 pos.'


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
targetmat=list(c(0.0,0.0,37.35541857539394,30.523299096923452,100.00889134859477,8.697228452262904,0.0),c(611.519045736028,626.2560828203974,607.9411258348425,612.0724713646229,650.6392408085903,487.0447933267226,600.0491808943186),c(468.9292682446946,608.5452207496849,485.62044148012126,1155.0658974046296,1037.882971263033,1083.4261729104644,901.6612030369656),c(515.9691948397736,467.5667586668125,186.04463368921688,607.2530030861612,726.2273563627608,561.5924657746903,303.72859773663043),c(1827.207148677603,767.9429793860982,1006.3989239723779,1182.3762176492453,468.0648693931325,822.5093193425774,1515.4681252921769),c(1469.9977060962212,848.7045104285476,823.2841270341723,836.9809910261641,537.2570674773347,464.68049159233226,643.4389805709802),c(214.6196650900483,512.9065655678368,226.32988901562211,359.85363145846594,886.1252931119676,277.0688492649468,649.7887073529306),c(858.4786603601932,536.993337984006,733.1916469405752,197.5981994169255,337.82073182286956,212.46086647670805,292.08743196972125),c(195.50969491079744,393.88957245264817,250.50104221146526,372.70554686769685,1070.4440056555986,603.8361468285386,554.5428056236736),c(527.7291764885434,347.1328965859669,270.27744028079144,369.4925680153891,347.7053315491841,140.39811644367256,337.59380724036623),c(1859.54709821172,1141.9963863195483,1262.027180498113,1405.6782478846326,1362.9118681459659,1385.3442463247338,1108.0273234503557),c(292.52954351314804,413.72573797184623,230.72464414213906,345.39522662308116,136.05860799750684,178.91441387512256,125.9362478420175),c(560.0691260226603,643.96694489111,547.1470132513583,1526.1649548461726,1792.6012327192889,1519.5300567310758,1300.635702502853),c(859.9486580662895,592.9596621274577,728.7968918140583,353.42767375385046,384.33649524082057,228.61286217376772,858.2714033603041),c(67.61989448042618,145.22906897984328,71.04854121202376,923.7314200384728,648.3134526376929,1309.5541126693,159.8014573457533),c(680.6089379225505,793.4466207679243,732.4591877528223,1625.7672992677121,1547.230580689597,853.5708495292306,1519.7012764801439),c(424.82933706180796,523.5330828102643,345.7207366193322,128.51915409230926,131.40703165571173,576.5020002642839,184.1420766765634),c(1342.10790566585,894.7527518124003,1278.141282628675,1044.2181270000128,1056.4892766302132,1334.4033368186226,1057.229509194752),c(1787.517210613005,587.2921862648296,1034.2323731069853,395.196398833851,573.3067841262467,149.09534489593548,743.9763212851958))
targetgene="Ercc2"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(1031.9383896795473,609.2536552325133,660.6781873530458,435.35863448769766,234.90460526065283,1268.5528928229178,420.14025540572226),c(210.20967197175963,529.9089931557209,474.6335536638289,1486.002719192326,1051.2562532456939,2438.951350256011,1528.1675788560779),c(170.51973390716168,221.73999312532172,75.4432963385407,819.3096073384716,786.1164017633728,1222.581828146671,1085.803279713529),c(135.23978896085237,310.2943034788847,167.0006948076435,345.39522662308116,526.7910207082957,208.73348285430967,226.47358855623315),c(74.96988301090728,259.99545519806094,103.27674547314795,383.95097285077395,345.3795433782866,341.6768320531855,530.2021862928636),c(223.43965132662564,309.5858689960562,168.46561318314915,565.4842780061608,695.4106630983682,966.6348194086484,460.3551916914085))
targetgene="Myd88"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(154.34975914010323,277.7063172687735,142.82954161180035,1539.0168702554035,1042.534547604828,1408.9510092665903,702.7030972025178),c(74.96988301090728,162.93993105055588,58.59673502022579,282.7421390030804,640.1731940395514,457.22572434753545,329.1275048644323),c(393.9593852337873,503.69691729106626,440.20797183944626,428.9326767830822,537.2570674773347,575.2595390568177,491.04553780416904))
targetgene="Smo"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(146.99977060962212,253.61954485260438,101.07936790988948,560.6648097276992,544.8158790327517,455.9832631400693,479.4043720372599),c(141.11977978523726,476.7764069435831,260.02301165225197,464.2754441584672,343.05375520738903,552.8952373224274,345.0018218193084),c(126.41980272427503,202.6122620889521,130.37773542000238,448.21054989692857,380.8478129844743,242.27993545589516,353.4681241952424),c(67.61989448042618,231.65807588492075,60.06165339573143,102.81532327384741,192.45897114177248,242.27993545589516,158.74316954876156))
targetgene="Rbm10"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(213.14966738395208,374.7618414162785,240.97907277067856,965.5001451184734,393.0582008816864,593.8964571688097,390.5081970899534),c(151.40976372791079,463.3161517698415,258.5580932767463,711.6748157861625,869.2633288729603,729.3247287826177,448.7140259244993),c(176.39972473154654,398.84861383244765,169.9305315586548,318.0849063784654,514.5806328110835,387.64789672943226,398.9744994658874))
targetgene="Tcf3"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(460.10928200811725,522.1162138446074,530.3004519330434,507.6506586646216,441.899752470535,408.76973725635645,427.54826998466444),c(146.99977060962212,269.91353795765997,184.57971531371123,1209.686537893861,985.552737417838,1468.5891472249646,1006.4316949391482),c(470.3992659507908,312.4196069273702,398.45779813753535,509.25714809077544,1028.5798185794426,890.8446857532144,557.7176690146489),c(161.69974767058434,201.1953931232951,194.10168475449794,764.6889668492402,992.5301019305306,473.37772004459515,1302.7522780968366))
targetgene="Ret"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(101.42984172063927,63.05066897173684,76.17575552629353,846.6199275830872,629.7071472705125,539.2281640402999,1523.9344276681109))
targetgene="Syk"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(977.5484745539871,498.7378759112667,547.1470132513583,255.43181875846466,372.70755438633284,161.51995697059678,303.72859773663043),c(1339.1679102536575,557.5379379860325,828.4113413484421,277.9226707246188,529.6982559219176,283.2811553022774,917.5355199918417),c(171.9897316132579,204.7375655374376,117.9259292282044,364.67309973692755,1462.3393124518363,583.9567675090806,498.4535523831113),c(80.84987383529217,216.78095174552217,79.10559227730482,375.9185257200046,689.0147456284,426.16419416088223,115.35336987210006),c(98.48984630844683,244.40989657583384,201.42627663202614,313.26543810000385,359.9157194463963,134.18581040634194,83.60473596234775),c(127.88980043037125,313.1280414101987,138.43478648528344,735.7721571784706,794.8381074042386,1198.9750652048144,996.9071047662226))
targetgene="Keap1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(27.929956415828205,111.93264828690361,43.947551265169345,117.27372810923221,218.62408806436997,296.94822858440483,115.35336987210006),c(161.69974767058434,211.11347588289414,102.54428628539513,318.0849063784654,108.7305969894606,79.51751727783225,142.8688525938854),c(30.869951828020646,194.11104829501005,81.30296984056328,151.01000605846338,213.39106467985047,145.36796127353708,214.83242278932397),c(85.25986695358084,179.23392415561148,90.82493928134997,380.7379939984662,170.36398351824573,241.03747424842902,245.52276890208452))
targetgene="Foxo1"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
targetmat=list(c(19.109970179250876,57.38319310910881,47.609847203933455,263.464265889234,150.59478406561655,246.00731907829353,107.94535529315786),c(149.93976602181456,238.0339862303773,112.06625572618182,501.22470096000615,384.33649524082057,539.2281640402999,1512.2932619012017),c(73.49988530481106,154.43871725661384,131.84265379550803,571.9102357107762,152.9205722365141,100.63935780475644,222.24043736826619))
targetgene="Arid5b"
collabel=c("input_rep1","input_rep2","input_rep3","spleen1","spleen2","spleen3","spleen4")

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
Sweave("ABE_spleen_input_summary.Rnw");
library(tools);

texi2dvi("ABE_spleen_input_summary.tex",pdf=TRUE);

