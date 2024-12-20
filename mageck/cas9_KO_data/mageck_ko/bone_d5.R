pdf(file='bone_d5.pdf',width=4.5,height=4.5);
gstable=read.table('bone_d5.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Trp53","Rara","Inppl1","Dis3","Erbb2","B2m","Irs1","Esr1","Cd79a","Ccnd3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='B1,B2,B3,B4_vs_d5_rep2,d5_rep3 neg.'


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
targetmat=list(c(304.09123474206314,211.40101829618703,728.8470935181421,1311.7889019009494,636.3180152639188,2618.388632934809),c(86.7885659067762,124.5463363605443,161.5252112384685,0.0,2838.5544384612676,376.08588297685776),c(1912.6609906325411,2241.8340544897974,658.002702624077,356.1577394946126,206.21417161330703,555.3826411402434),c(1016.9499898236753,929.1812199530082,190.4297227232471,29.97367114558621,356.128601135362,108.23401864740964),c(1040.8002827446214,1233.9919905196034,1712.1672391277662,574.7892231447709,604.2402552351822,733.5861263879987),c(451.1680410878976,213.03978587987842,358.7559954875458,329.7103826014483,542.703327833116,403.41770586761777),c(174.9021480869383,73.74454126611175,437.53495816174626,708.7891647368033,501.46049351045457,569.5951890434386),c(323.9664788428516,165.51552595282863,397.8620992610698,68.76312792222718,228.47220918426717,205.5353081385153),c(738.0340642759442,481.7976696052635,1556.309579160823,442.55243867894933,849.0786685157436,1789.6877628869654),c(319.9914300226939,152.40538528329765,434.1344273988311,555.3944947564504,250.07559859137552,188.04294148842888),c(203.38999796473504,183.54196937343372,454.5376119763219,262.7104118054321,309.64858150188644,294.09041430457773),c(157.67693653292164,406.41436075546034,90.11406521725085,14.105257009687628,349.5821194968443,144.31202486321286),c(127.86407038173897,170.43182870390274,0.0,178.0788697473063,0.0,263.4787726669265),c(208.0275549215857,149.1278501159149,651.7683962253992,664.7102365815294,783.6138521305668,112.60711030993124),c(641.3078763187738,475.24259927049803,653.4686616068568,551.8681805040285,369.87621257624915,642.8444743906755),c(494.8935781096322,299.8944678155212,329.28472887561475,227.447269281213,1044.1638213435706,965.3599845016436),c(325.9540032529304,471.96506410311525,1839.6871427370834,257.42094042679923,110.63553969094886,135.56584153816965),c(407.44250406616305,337.5861222404227,1453.7269011462165,144.5788843492982,105.3983543801347,781.6901346757363),c(301.44120219529134,181.90320178974233,218.76747908087313,74.05259930086005,72.6659461875463,57.94346452841122),c(402.8049471093124,534.2382322833874,256.74007260009205,357.92089662082356,251.38489491907904,185.85639565716806),c(659.8581041461763,965.2341067942184,81.61273830996304,24.68419976695335,103.43440988857941,5.466364578152002),c(682.3833807937365,647.3131955580922,319.64989171402186,313.8419684655497,547.2858649800784,5154.781797197338),c(329.92905207308814,683.3660823993023,457.3713876120845,223.9209550287911,851.6972611711507,382.6455204706401),c(1010.3249084567458,1155.3311465024176,554.2865143551657,426.6840245430508,1361.0135326478264,1133.7240135087252),c(2337.9912143894135,3362.751081734696,595.0928835101471,4131.077146712264,175.44570791227392,432.93607458963857),c(1188.5395972271488,1548.6353665883469,999.1892891698947,44.078928155273836,1135.814564282818,1046.2621802582933),c(773.8095036573633,855.4366786868964,2295.9250200948627,504.26293809633273,26.840574717922504,97.30128949110563),c(1128.2513567880906,927.5424523693168,1729.7366480694943,516.6050379798094,510.62556780437933,706.2543034972387),c(588.9697335200308,394.9429876696208,3721.880920010606,1937.709681705838,1198.0061398487362,3196.730005303291),c(1091.1509011332855,860.3529814379705,129.2201689907748,290.92092582480734,147.95048503049966,206.62858105414568),c(456.4681061814412,562.0972812061408,358.18924036039334,151.63151285414202,756.1186292487924,741.2390367974115),c(80.82599267653967,26.22028133906196,13.602123051660506,0.0,0.6546481638517684,1.0932729156304004),c(375.64211350490154,96.68728743779097,0.0,7.052628504843814,40.58818615880964,0.0),c(386.90475182868164,311.36584090136074,301.5137276451412,791.6575496687182,1229.429251713621,2024.7414397475015),c(109.97635069102938,111.43619569101332,5.667551271525211,0.0,32.732408192588416,0.0),c(204.052506101428,342.50242499149687,69.71088063976009,72.2894421746491,97.5425764139135,7.652910409412803),c(810.9099593121684,796.441045674007,75.94518703843782,33.499985398008114,696.5456463382816,135.56584153816965),c(639.9828600453878,553.9034432876839,183.0619060702643,389.65772489262076,421.5934175205388,916.1627032982756),c(218.62768510867286,88.4934495193341,403.529650532595,419.63139603820696,170.86317076531154,119.16674780371365),c(268.9783034973369,111.43619569101332,388.22726209947695,558.9208090088723,1240.5582704991011,974.1061678266867),c(247.77804312316258,273.67418647645917,1014.4916776030127,1888.3412821719312,1671.3167623135646,3010.8736096461225),c(571.7445219660142,986.5380853822062,488.5429196054732,990.8943049305559,2159.6842925469837,752.1717659537155),c(284.21599064127474,478.52013443788076,358.18924036039334,306.78933996070595,1095.2263781240085,53.57037286588962),c(433.94282953388097,701.3925258199074,136.58798564375758,283.8682973199635,225.8536165288601,111.51383739430084),c(115.93892392126591,114.71373085839608,493.64371574984585,61.710499417383375,104.74370621628294,142.12547903195204),c(13.91267087055191,4.9163027510741175,54.97524733379454,119.89468458234484,43.20677881421671,3.279818746891201),c(446.530484131047,416.2469662576086,601.3271899088248,3.526314252421907,20.94874124325659,499.62572244309297),c(2106.775874683575,3230.010907455695,899.440386791051,638.2628796883652,72.6659461875463,1042.982361511402),c(172.2521155401665,254.00897547216272,108.81698441328405,121.6578417085558,66.77411271288038,126.81965821312644),c(253.07810821670617,372.00024149794154,321.35015709547946,364.9735251256674,720.767628400797,432.93607458963857),c(527.3564768075867,530.9606971160047,58.94253322386219,109.31574182507912,6329.793096282749,142.12547903195204),c(131.83911920189666,77.0220764334945,339.4863211643601,262.7104118054321,56.29974209125208,14.212547903195205),c(162.31449348977227,195.01334245927333,335.5190352742925,79.3420706794929,104.74370621628294,329.0751476047505),c(382.267194871831,296.6169326481384,652.9019064797043,491.92083821285604,453.6711775492755,364.0598809049233),c(109.97635069102938,101.60359018886508,616.0628232147905,93.44732768918054,542.0486796692642,1888.0823252937016),c(297.4661533751337,188.45827212450783,1089.870109514298,3.526314252421907,206.86881977715882,67.78292076908482),c(345.16673921702596,88.4934495193341,69.71088063976009,137.52625584445437,243.52911695285783,414.35043502392176),c(98.71371236724926,62.27316818027215,82.74624856426807,0.0,20.94874124325659,2.186545831260801),c(191.464851504262,186.81950454081647,184.76217145172188,96.97364194160245,98.19722457776525,83.08874158791043),c(769.8344548372056,665.3396389786972,1675.3281558628523,368.4998393780893,469.38273348171793,69.96946660034563),c(47.03807770519931,36.052886841210196,32.3050422476937,3.526314252421907,35.35100084799549,59.03673744404162),c(736.7090480025582,979.9830150474407,621.1636193591631,892.1575058627425,197.04909731938227,533.5171828276355),c(1758.959102919777,2990.750840236755,1996.1115578311792,1729.6571408129455,1936.449268673531,3438.3433196576093),c(816.872532542405,799.7185808413898,719.7790114837018,495.44715246527795,1961.9805470637498,1493.410802751127))
targetgene="Trp53"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(402.8049471093124,221.2336237983353,78.77896267420043,86.39469918433673,192.46656017241992,24.05200414386881),c(2917.6858339957435,3815.050934833515,188.72945734178953,267.9998831840649,1040.2359323604599,237.2402226917969),c(189.47732709418315,131.1014066953098,28.904511484778574,29.97367114558621,14.402259604738905,0.0),c(958.6492737946959,1253.6572015238999,98.61539212453867,125.1841559609777,868.7181134312966,1151.2163801588117))
targetgene="Rara"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(493.5685618362463,771.8595319186364,35.70557301060883,126.94731308718866,92.30539110309934,56.85019161278082),c(738.0340642759442,1004.5645288028113,242.57119442127902,879.8154059792658,53.02650127199324,79.80892284101922))
targetgene="Inppl1"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(106.66381000756463,65.5507033476549,1.1335102543050422,199.23675526183774,51.06255678043794,24.05200414386881),c(243.14048616631194,385.1103821674725,104.8496985232164,15.868414135898583,27.49522288177427,854.9394200229731),c(1282.6157526375475,2097.622507124957,567.3218822796736,740.5259930086005,1810.7568212139913,447.1486224928338),c(538.6191151313668,693.1986879014505,325.8841981126996,398.4735105236755,163.00739279909033,322.51551011096814),c(100.03872864063516,208.1234831288043,29.471266611931096,12.342099883476674,32.732408192588416,417.63025377081294))
targetgene="Dis3"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(1766.2466924233995,2358.1865529318848,265.8081546345324,458.42085281484793,539.4300870138571,88.55510616606243),c(404.79247151939126,619.4541466353388,239.1706636583639,211.57885514531444,335.1798598921054,237.2402226917969),c(1079.8882628095052,1243.8245960217516,44.20689991789664,342.052482484925,398.68073178572695,2198.5718333327354),c(1105.0635720038374,1443.754241232099,2035.2176616047032,770.4996641541867,1170.5109169669618,2533.113345515638),c(221.27771765544466,260.5640458069282,211.96641755504288,548.3418662516066,127.00174378724307,648.3108389688274),c(36.437947518112146,4.9163027510741175,51.0079614437269,238.02621203847872,28.80451920947781,89.64837908169284),c(255.72814076347794,88.4934495193341,1862.9241029503369,701.7365362319595,375.7680460509151,1141.3769239181381),c(1515.1561086167721,1560.1067396741867,304.91425840805636,481.3418954555903,70.70200169599099,160.71111859766887))
targetgene="Erbb2"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(548.556737181761,445.7447827640533,1021.292739128843,692.9207506009047,470.69202980942146,1358.9382341285877),c(1194.5021704573853,1561.745507257878,2179.1734639014435,516.6050379798094,951.2037820766194,3733.5270068778173),c(1124.2763079679328,1191.3840333436278,29.471266611931096,359.6840537470345,197.70374548323406,178.20348524775525),c(394.85484946899703,598.1501680473509,495.3439811313034,204.5262266404706,89.68679844769227,44.82418954084642),c(127.86407038173897,93.40975227040823,34.57206275630379,0.0,0.0,0.0))
targetgene="B2m"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(637.332827498616,955.4015012920702,836.5305676771211,130.47362733961057,877.8831877252214,401.2311600363569),c(1584.0569548328388,2497.4817975456517,1070.0336800639598,1955.3412529679474,758.0825737403478,2752.861201557348),c(292.82859641828304,126.18510394423568,67.44386013115,236.26305491226776,64.81016822132507,397.95134128946574),c(2093.525711949716,1896.054094330918,885.8382637393904,276.8156688151197,3656.2099951121263,1002.5312636330772),c(1863.635388517263,2543.36728988901,2731.7597128751518,1763.1571262109535,133.54822542576076,809.0219575664963),c(1037.4877420611567,1181.5514278414796,234.63662264114373,1512.788814288998,227.16291285656362,13754.466551546067))
targetgene="Irs1"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(606.8574532107405,1071.7539997341576,366.69056726768116,423.15771029062887,102.77976172472763,1.0932729156304004),c(1255.4529190331366,1258.573504274974,494.7772260041509,1481.051986017201,369.22156441239736,1296.6216779376548),c(1242.2027562992776,2174.6445835584514,649.5013757167892,105.78942757265722,691.9631091913192,2652.280093319351),c(963.9493388882394,866.908051772736,340.61983141866517,486.6313668342232,7.2011298023694525,1235.3983946623525),c(1480.0431773720459,1255.2959691075912,385.39348646371434,749.3417786396552,1887.3506563846483,936.9348886952531))
targetgene="Esr1"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(1018.2750060970611,1451.948079150556,303.7807481537513,497.2103095914889,147.29583686664787,1466.0789798603669),c(1638.3826220416606,2962.8917913140012,79.34571780135295,296.2103972034402,186.574726697754,112.60711030993124),c(147.73931448252742,144.21154736484078,37.97259351921891,158.6841413589858,46.48001963347556,1.0932729156304004),c(280.240941821117,188.45827212450783,70.84439089406513,220.3946407763692,221.92572754574948,2008.3423460130455))
targetgene="Cd79a"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(361.72944263434965,1224.1593850174552,137.1547407709101,183.36834112593917,604.894903399034,261.2922268356657),c(680.3958563836577,888.2120303607238,106.54996390467396,282.10514019375256,153.1876703413138,209.90839980103686),c(594.9323067502673,1155.3311465024176,141.68878178813026,973.2627336684463,912.5795404093651,68.87619368471523),c(1568.156759552208,2430.2923266143052,241.43768416697398,881.5785631054767,475.27456695638386,153.05820818825606))
targetgene="Ccnd3"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetgenelist=c("ST","Ptprd","Cdk12","Eif1ad19","Map2k4","Fbxw7","Braf","Erg","Atr","Arid2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='B1,B2,B3,B4_vs_d5_rep2,d5_rep3 pos.'


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
targetmat=list(c(110.63885882772233,72.10577368242039,333.2520147656824,63.47365654359433,927.6364481779558,62.31655619093282),c(1.3250162733858961,3.277535167382745,100.31565750599623,44.078928155273836,68.73805720443568,348.75406008609775),c(23.187784784253182,3.277535167382745,218.76747908087313,19.394728388320488,0.0,7.652910409412803),c(172.91462367685943,506.3791833606341,151.89037407687564,84.63154205812577,185.92007853390223,181.48330399464646),c(599.569863707118,863.6305166053533,578.0902296955715,163.97361273761868,279.5347659647051,92.92819782858403),c(180.20221318048186,180.26443420605096,900.0071419182035,172.78939836867346,57.60903841895562,781.6901346757363),c(25.837817331024976,6.55507033476549,234.0698675139912,1.7631571262109536,110.63553969094886,1.0932729156304004),c(127.86407038173897,45.88549234335843,13.035367924507984,10.578942757265722,385.58776850869157,52.477099950259216),c(11.925146460473066,3.277535167382745,847.298915093019,190.420969630783,623.8797001507353,349.84733300172815),c(5.962573230236533,1.6387675836913724,101.44916776030128,88.15785631054767,62.84622372976976,2.186545831260801),c(425.33022375687267,208.1234831288043,981.053125101014,1465.1835718813024,845.8054276964848,651.5906577157186),c(178.21468877040303,101.60359018886508,2158.2035241968,1278.2889165029412,1160.0365463453336,558.6624598871346),c(25.837817331024976,96.68728743779097,176.26084454443406,158.6841413589858,136.8214662450196,40.451097878324816),c(2.6500325467717922,3.277535167382745,53.841737079489505,3.526314252421907,45.17072330577202,20.77218539697761),c(257.05315703686387,163.87675836913724,159.82494585701093,393.1840391450426,1577.7020748827617,131.19274987564805),c(328.6040357997022,136.01770944638392,986.7206763725392,423.15771029062887,155.80626299672087,556.4759140558738),c(66.2508136692948,27.85904892275333,1147.6791324838553,121.6578417085558,621.2611074953282,75.43583117849762),c(49.6881102519711,11.471373085839607,2898.385720257993,148.1051986017201,167.5899299460527,774.0372242663235),c(127.20156224504603,167.15429353652,221.03449958948323,167.4999269900406,202.9409307940482,101.67438115362724),c(151.71436330268511,103.24235777255646,473.8072862995076,495.44715246527795,553.1776984547442,203.34876230725447),c(231.87784784253182,294.978165064447,1164.1150311712784,405.5261390285193,1203.2433251595503,385.92533921753136),c(3.9750488201576886,0.0,36.83908326491387,10.578942757265722,22.912685734811895,0.0),c(34.4504231080333,13.11014066953098,228.96907136961852,17.631571262109535,148.60513319435142,26.238549975129608),c(0.0,0.0,1.1335102543050422,0.0,2.6185926554070735,1.0932729156304004),c(208.0275549215857,201.5684127940388,2289.690713696185,610.05236566899,1526.6395181023238,999.2514448861859),c(285.5410069146606,132.74017427900117,1673.6278904813948,248.60515479574445,96.23328008620996,114.79365614119205))
targetgene="ST"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(129.18908665512487,26.22028133906196,775.3210139446488,375.5524678829331,828.129927272487,2190.9189229233225),c(66.2508136692948,18.026443420605098,74.24492165698027,518.3681951060204,92.96003926695111,980.6658053204692),c(27.16283360441087,4.9163027510741175,95.78161648877607,144.5788843492982,1093.917081796305,166.17748317582087),c(0.0,0.0,17.569408941728152,65.23681366980529,13.747611440887136,1.0932729156304004))
targetgene="Ptprd"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(102.026253050714,45.88549234335843,965.1839815407434,661.1839223291075,1414.0400339198197,1402.6691507538037),c(35.7754393814192,24.581513755370587,2000.0788437212468,88.15785631054767,89.68679844769227,2018.181802253719))
targetgene="Cdk12"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(704.2461493046038,314.6433760687435,1164.1150311712784,396.71035339746453,854.9705019904095,873.5250595886899),c(919.561293729812,788.2472077555501,15165.23369234716,1680.2887412790387,591.1472919581469,4255.018187633518),c(593.6072904768814,390.02668491854666,468.1397350279824,601.2365800379351,876.5738913975179,2783.4728431949993),c(117.92644833134476,227.78869413310076,1104.038987693111,1172.499488930284,9.819722457776527,894.2972449856676),c(93.41364727370568,85.21591435195137,816.1273830996304,211.57885514531444,1823.1951363271749,693.1350285096738),c(46.375569568506364,40.96918959228431,60.642798605319754,345.5787967373469,142.7132997196855,332.35496635164174),c(78.83846826646082,19.66521100429647,570.7224130425888,195.71044100941583,7626.651108873101,1380.8036924411956))
targetgene="Eif1ad19"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(38.42547192819099,11.471373085839607,315.1158506968017,192.18412675699395,39.278889831106106,374.99261006122737),c(55.65068348220764,3.277535167382745,256.1733174729395,19.394728388320488,339.10774887521603,86.36856033480163),c(105.33879373417874,16.387675836913726,181.36164068880674,611.8155227952009,155.15161483286911,258.0124080887745),c(94.73866354709158,22.942746171679214,269.20868539744754,29.97367114558621,385.58776850869157,60.13001035967202),c(17.22521155401665,8.193837918456863,536.150350286285,121.6578417085558,85.10426130072989,416.53698085518255))
targetgene="Map2k4"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(276.26589300095935,168.79306112021135,680.6729077101778,405.5261390285193,559.0695319294102,444.96207666157295),c(180.20221318048186,119.63003360947019,1096.1044159129758,1214.815259959347,1667.388873330454,2138.441822973063),c(260.36569772032857,131.1014066953098,712.9779499578715,229.21042640742397,651.3749230325095,681.1090264377394),c(757.9093083767326,1043.8949508114042,533.3165746505223,269.76304031027587,754.1546847572372,260.1989539200353),c(593.6072904768814,635.8418224722525,913.609264969864,911.5522342510629,857.5890946458165,165.08421026019047),c(1497.2683889260627,1456.8643819016302,1189.0522567659891,1394.6572868328642,311.61252599344175,165.08421026019047),c(259.70318958363566,227.78869413310076,1680.9957071343774,1354.1046729300124,666.4318308011002,2452.211149758988),c(537.2940988579809,703.0312934035987,998.0557789155896,613.5786799214119,2627.7577297009984,1152.309653074442),c(604.8699288006616,283.50679197860745,527.0822682518447,1364.683615687278,498.8419008550475,2584.4971725502664),c(305.41625101544906,298.2557002318298,1228.7251156666657,1283.5783878815741,703.746776140651,1773.2886691525096),c(151.05185516599215,49.16302751074117,276.00974692327776,262.7104118054321,404.5725652603929,406.69752461450895),c(205.3775223748139,267.1191161416937,1141.4448260851775,195.71044100941583,26.840574717922504,5431.379844851829))
targetgene="Fbxw7"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(132.5016273385896,101.60359018886508,563.9213515167585,91.68417056296958,342.38098969447486,554.289368224613),c(148.40182261922035,32.77535167382745,1678.7286866257675,1264.1836594932538,692.617757355171,956.6138011766003),c(80.82599267653967,27.85904892275333,498.74451189421853,453.13138143621507,197.70374548323406,309.3962351234033),c(319.9914300226939,214.67855346356978,328.1512186213097,375.5524678829331,466.76414082631084,437.30916625216014),c(483.6309397858521,350.6962629099537,440.9354889246614,800.4733352997729,655.3028120156201,910.6963387201235),c(165.627034173237,52.44056267812392,901.707407299661,618.8681513000447,2172.1226076601674,3632.9458986398204),c(2420.142223339339,2359.825320515576,1904.8639823596234,4280.945502440195,655.3028120156201,1494.5040756667574),c(378.29214605167334,447.3835503477447,608.1282514346551,599.4734229117242,32.732408192588416,298.4635059670993),c(355.10436126742013,232.70499688417488,1233.259156683886,1955.3412529679474,4212.660934386129,12743.189104587947),c(151.71436330268511,37.691654424901564,398.4288543882223,366.73668225187834,616.023922184514,740.145763881781),c(180.20221318048186,129.46263911161842,327.0177083670047,192.18412675699395,89.0321502838405,501.8122682743538))
targetgene="Braf"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(51.013126525357,63.91193576396353,190.9964778503996,225.68411215500205,260.5499692130038,315.9558726171857),c(21.200260374174338,0.0,241.43768416697398,273.2893545626978,80.52172415376751,178.20348524775525),c(59.62573230236533,27.85904892275333,773.6207485631912,181.6051839997282,2057.5591789861082,545.5431848995698))
targetgene="Erg"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(617.4575833978275,857.0754462705878,339.4863211643601,61.710499417383375,134.85752175346428,1495.5973485823877),c(568.4319812825495,1288.0713207814188,1052.4642711222316,507.78925234875464,285.426599439371,144.31202486321286),c(330.5915602097811,195.01334245927333,1034.328107053351,833.973320697781,145.3318923750926,851.659601276082),c(16.5627034173237,60.63440059658078,2905.7535369109755,389.65772489262076,2167.540070513205,3074.283438752686))
targetgene="Atr"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
targetmat=list(c(144.4267737990627,39.33042200859294,105.41645365036892,112.84205607750103,118.49131765717007,220.84112895734089),c(60.95074857575122,24.581513755370587,1046.2299647235538,373.78931075672216,1520.747684627658,882.2712429137331),c(95.40117168378453,249.0926727210886,966.3174917950485,350.86826811597973,818.9648529785622,378.27242880811855),c(49.025602115278154,24.581513755370587,600.1936796545198,160.44729848519677,675.596905095025,2202.9449249952568),c(642.6328925921596,489.9915075237204,718.6455012293967,343.81563961113596,339.76239703906776,275.50477473886093))
targetgene="Arid2"
collabel=c("d5_rep2","d5_rep3","B1","B2","B3","B4")

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
Sweave("bone_d5_summary.Rnw");
library(tools);

texi2dvi("bone_d5_summary.tex",pdf=TRUE);

