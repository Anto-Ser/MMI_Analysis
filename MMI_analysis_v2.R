# Section0: Libraries ----------------------------------------

Packages = c("dplyr", "stringr", "readxl", "tidyr", "heatmaply", "ggplot2",
             "reshape2", "lares", "corrr", "corrplot", "viridis", "grid", "gridExtra", "gridGraphics",
             "circlize", "kableExtra", "janitor", "RColorBrewer", "randomcoloR", "packcircles", "ggbeeswarm", "ggpubr",
             "scales", "networkD3", "gridGraphics", "grid", "ggcorrplot", "vegan", "plyr", "tibble","ggthemes", "circlize",
             "writexl")
package.check = lapply(Packages,FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE) }})

lapply(Packages, library, character.only = TRUE)  #package not installed:"ComplexHeatmap"

# Section1: Functions -----------------------------------------------------

# --- Merge PCR to merge replicate together ---
merge_PCR = function(data,min.rep=2,CPM=TRUE,show.corr=FALSE){
  #Code to generate filtered data by replicate
  #min.rep: is the minimal number of replicates
  #CPM: normalize to counts per million (CPM) prior and after merging
  #note:remove samples with zero counts before running this
  samples <- substr(names(data),1,nchar(names(data))-1) 
  merged = data.frame(row.names = row.names(data))
  
  for(s in unique(samples)){
    print(s)
    d = data[,grep(s,names(data))]
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)>=2 & show.corr==TRUE) {corr=cor(d,use="pair",method="p"); print(heatmap4(as.data.frame(corr)));}
    if(ncol(d)<=1){merged=cbind(merged,0)
    } else {merged=cbind(merged,(rowSums(d>0)>=min.rep)*rowSums(d))}
  }
  if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
  names(merged)=unique(samples)
  if(sum(is.na(colSums(merged)))>0){merged = merged[,-which(is.na(colSums(merged))==TRUE)]} # Remove column with NA values
  return(merged)
}

# --- Merge Tumour samples ---
merge_tumor = function(data,min.rep=1,CPM=TRUE, Remove_NA_col=FALSE){
  tum.pieces = names(select(data,contains("Tum")))
  tum.samples <- paste0(apply(str_split_fixed(tum.pieces,"_",3)[,1:2],1, paste,collapse="_"),"_Tum")
  tum.merged = data.frame(row.names = row.names(data))
  
  for(s in unique(tum.samples)){
    print(s)
    d = data[,grep(s,names(data))]
    if(Remove_NA_col==TRUE){d=d[!is.na(colSums(d))]}
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)<=1){tum.merged=cbind(tum.merged,0)}
    else {tum.merged=cbind(tum.merged,(rowSums(d>0)>=min.rep)*rowSums(d))}
  }
  names(tum.merged)=unique(tum.samples)
  #tum.merged=tum.merged[,colSums(tum.merged)>0]
  if(CPM==TRUE) tum.merged=sweep(tum.merged,2,colSums(tum.merged)/1e6,`/`)
  
  ##now merge this with the non.tumor samples
  return(merge(tum.merged,select(data,!contains("Tum")),by="row.names",all=TRUE)%>%tibble::column_to_rownames(var="Row.names") )
}

#--- Merge experiment or mouse as single sample by row or col
merge_exp = function(data,grp=c("exp","mouse"),by = c("col", "row"), min.rep=1,CPM=TRUE){
  exp = str_split_fixed(names(data),"_",3)[,1]
  mice = str_split_fixed(names(data),"_",3)[,2]
  mice = mice[nchar(mice)>0]
  mice = paste(exp, mice, sep = "_")
  
  df.merged = data.frame(ID= 1:nrow(data), row.names = row.names(data))
  exp.merged = data.frame()
  if(by == "row") {df1 = data.frame()}
  if(by == "col") {df1 <- as.data.frame(matrix(NA, nrow = nrow(data)))}
  
  if(grp == "exp"){
    for(s in unique(exp)){
      print(s)
      d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
      
      if( by == "row"){
        if(ncol(d)<=1){next}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = rbind(df1, df.merged)
          #df.merged=cbind(df.merged,(rowSums(d>0)>=min.rep)*rowSums(d))
        }
      }
      if( by == "col"){
        if(ncol(d)<=1){df1=cbind(df1,0)}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = cbind(df1, df.merged$n)
        }
      }
    }
    if (by == "col"){
      df1 = df1[,-1]
      names(df1)= unique(exp)}
  }
  if(grp == "mouse"){
    for(s in unique(mice)){
      print(s)
      d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
      
      if( by == "row"){
        if(ncol(d)<=1){next}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = rbind(df1, df.merged)
          #df.merged=cbind(df.merged,(rowSums(d>0)>=min.rep)*rowSums(d))
        }
      }
      if( by == "col"){
        if(ncol(d)<=1){df1=cbind(df1,0)}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = cbind(df1, df.merged$n)
        }
      }
    }
    if (by == "col"){ df1 = df1[,-1]
    names(df1)= unique(mice)}
  }
  if(by == "col") {df1=df1[,colSums(df1)>0]}
  return(df1)
}

#--- Function to filter data per mouse by presence of barcode in primary tumour:
Filter_by_Tum = function(data, CPM=TRUE, filter=TRUE){
  exp = str_split_fixed(names(data),"_",3)[,1]
  mice = str_split_fixed(names(data),"_",3)[,2]
  mice = mice[nchar(mice)>0]
  mice = paste(exp, mice, sep = "_")
  df.merged = data.frame(row.names = row.names(data))
  for(s in unique(mice)){
    print(s) #s = "Exp57_339"
    d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
    DA = which(d[,grep("Tum", names(d)), drop=F]>0)
    if(filter==TRUE) d[-DA,]=0
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    df.merged = cbind(df.merged,d)
  }
  return(df.merged)
}

#--- Merge samples I (ex:Lung LungI) as a single sample
merge_sampleI = function(data, CPM=TRUE){
  #data=df.PCR.Tum.Filt
  #CPM=TRUE
  Mouse.nb = str_split(names(data), "_", simplify = T)[,2]
  Mega.merged = data.frame(row.names = row.names(data))
  for (i in unique(Mouse.nb)) {
    #i="333"
    print(i)
    TEST.name = paste("_", i, "_", sep="")
    Ms.x = data[,grep(TEST.name, names(data)), drop=F]
    id.name <- paste(str_split(names(Ms.x[1]), "_", simplify = T)[,c(1,2)], collapse = "_")
    colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
    names(Ms.x)[which(endsWith(names(Ms.x), 'I'))] = str_sub(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))],1,nchar(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))])-1)
    
    samples <- names(Ms.x)
    merged = data.frame(row.names = row.names(Ms.x))
    for(s in unique(samples)){
      print(s)
      #s="Lung"
      d = Ms.x[,grep(s,names(Ms.x)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      d = d[,!is.na(colSums(d)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      if(ncol(d)==2){d=sweep(d,2,colSums(d)/1e6,`/`)}
      merged=cbind(merged,(rowSums(d>0)*rowSums(d)))
    }
    if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
    names(merged) = paste(id.name,unique(samples), sep = "_")
    Mega.merged = cbind(Mega.merged,merged)
  }
  return(Mega.merged)
}

#--- Stacked histogram function: needs a dataframe with one column ID
Stacked_histo = function(data, PCT=TRUE, angle.x=0, my_col = "Default"){
  if(PCT==TRUE) {data= data/10000}
  if(sum(names(data)=="ID")==0){data$ID = 1:nrow(data)}
  
  combo_color = c(brewer.pal(11, "Spectral"), brewer.pal(11, "RdYlGn"), brewer.pal(11, "RdBu"), brewer.pal(11, "PuOr"), brewer.pal(11, "PiYG"))
  getPalette = colorRampPalette(combo_color)
  if(my_col=="Default"){COLOR = sample(getPalette(2608), 2608)} else {
    COLOR = my_col
  }
  
  data$COLOR = COLOR
  COLOR.to.save = COLOR
  data = data[rowSums(data[1:(ncol(data)-2)])>0,] #remove row w/o bc (ignoring last 2 col (ID, colours))
  colnames(data) =gsub("^[^_]*_","",names(data))
  
  COLOR= data$COLOR
  data = data[,-ncol(data)]
  
  mA = melt(data, id=c("ID")) 
  mA$ID = as.factor(mA$ID) #transform barcode number in factor
  
  mA = mA %>% group_by(variable) %>% mutate(pos = sum(value)-(cumsum(value)-0.5*value))
  #Add new colunm in mA table for position of each barcode (to later on be labelled on graph)
  
  #colourCount = length(unique(mA$ID))
  
  
  #write.table(COLOR, file=paste("COLOR_exp51.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
  
  insert_minor <- function(major_labs, n_minor) {labs <- 
    c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
  labs[1:(length(labs)-n_minor)]}
  
  #PLot:
  p123 = ggplot(data =mA, aes(y = value, x = variable, fill = ID, label= ID)) + 
    geom_bar(stat="identity",alpha=0.9) +
    #geom_text(data = subset(mA, value>5), aes(y = pos, label = ID), size = 3) +
    theme_bw()+
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "grey10", face="bold",size=10,angle = 0),
          axis.text.x = element_text(color = "black", face="bold",size=10,angle = angle.x,
                                     vjust = 0.3),
          #hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "grey10", size = 15))+
    xlab("")+ylab("Number of reads(%)")+ ggtitle("")+
    scale_fill_manual(values = COLOR )+
    scale_y_continuous(breaks = seq(0,100, 10),
                       labels = insert_minor (seq(0, 100, by=20),1),
                       limits = c(0,100.001),
                       expand = c(0,0))
  my_list = list(Stack.plot = p123, Colour= COLOR.to.save)
  return(my_list)
}
Bubble_plot = function(data, PCT=TRUE, Y="avg", data2, angle.x=0, my_col=FALSE){
  
  if(PCT==TRUE){data = data/10000} #CPM to percent
  if(Y=="avg"){data$Y = rowSums(data)/ncol(data)} else if (Y=="Random"){#Average all sample for Y order or Random order
    data$Y= sample(length(data[,1]), length(data[,1]))
  } else if(Y=="input") {data$Y = data2 #manual input of previously saved Y
  } else{data$Y = data[,1]}#Take first column as Y ref = Tum collapsed
  colnames(data) =gsub("^[^_]*_","",names(data))
  data$ID = as.numeric(rownames(data))
  
  data1 = data[rowSums(data[1:(ncol(data)-2)])>0,1:(ncol(data)-2)] #remove empty rows
  data1$Y = data[rownames(data1),"Y"]
  data1$ID = data[rownames(data1),"ID"]
  
  mG =melt(data1, id=c("ID","Y"))
  mG$ID = as.factor(mG$ID) 
  mG[mG==0] = NA 
  
  Pallett = unname(distinctColorPalette(500))
  Max.color = sample(Pallett, length(unique(data$ID)),replace = TRUE)
  
  if(my_col==FALSE){col.bbp = Max.color[data1$ID]} else {
    col.bbp = my_col[data1$ID]
  }
  
  Plot.bbp = ggplot()+
    #geom_point(aes(x=variable,y=Y,color=ID,size=value), data=subset(mG,variable=!"X25k_P0_1M", stat="identity"),alpha=0.5)+
    geom_quasirandom(data=mG,aes(x=variable,y=Y,color=ID,size=value),alpha=0.8, width = 0.2)+
    scale_size_continuous(range = c(1.5,15))+
    #geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)+
    scale_color_manual(values = col.bbp,guide="none")+ #breaks=length(unique(mG$ID)) removed
    xlab("")+ylab("Frequency in Tumours")+ ggtitle("")+labs(size="Frq %")+
    #scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
    #                  breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
    #                 labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    theme_minimal()+
    theme( panel.grid.major.x = element_blank() ,
           panel.grid.minor.y = element_blank(),
           axis.text.x = element_text(angle = angle.x))
  
  if(Y=="avg"){Plot.bbp = Plot.bbp +
    scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
                       breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                       labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)} else if(Y=="Random"){
      
      Plot.bbp = Plot.bbp +
        scale_y_continuous(limits = c(0,length(data[,1])))}
  Plot.bbp
  
  
  my_list = list(BBP = Plot.bbp, Colour = Max.color, Y.order = data$Y)
  return(my_list)
}

biomass_fun = function(data){
  #data=Ms.x
  x = matrix(data=NA, nrow = ncol(data), ncol = ncol(data))
  colnames(x) = gsub("^[^_]*_","",names(data))
  rownames(x) = gsub("^[^_]*_","",names(data))
  
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      #print(c(i,j))
      if(sum(data[i])==0){x[i,j] =0} else {x[i,j] = sum(data[data[i]>0,][j])}
    }
  }
  return(x)
}
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

look.samples = function(data, sample.to.look){
  data[,grep(paste(sample.to.look, collapse = "|"), names(data))]
}
#Select the specific organ and order the sample by mode of injection
Select_and_order_Organs = function(data, Model, Organ){
  MI.tum = data[,grep(paste("_", Organ, sep = ""), names(data)), drop=F]
  
  order.tum.MDA= vector()
  for (i in unique(df.bc.info$MI)) {
    print(i)
    #i="IMFPc"
    order.temp = df.bc.info %>% 
      filter(Org == Organ) %>%
      filter(MI == i & Cells == Model) %>% arrange(Mouse) %>% .$Full_name
    order.tum.MDA = c(order.tum.MDA, order.temp)
  }
  Tum.order = MI.tum[,order.tum.MDA]
  newdf <- t(na.omit(t(Tum.order)))
  newdf = as.data.frame(newdf)
  Tum.order = newdf
}
# Section2: Pre - Figures analysis -----------------------------------------------------

info = read_xlsx("Info_All_mice_PhD.xlsx")
weight = read_xlsx("Info_All_mice_PhD.xlsx", sheet=2)

df.MI = read.csv("Data_MMI.csv")
df.MI = df.MI[,-1]

df.MI[df.MI<=10]=0

df.PCR = merge_PCR(df.MI)

df.PCR.Tum = merge_tumor(df.PCR)

df.PCR.Tum$id <- as.integer(row.names(df.PCR.Tum))#Reorder by id name
df.PCR.Tum = df.PCR.Tum[order(df.PCR.Tum$id), ]#Reorder by id name
df.PCR.Tum = df.PCR.Tum[,-ncol(df.PCR.Tum)]#Reorder by id name

df.PCR.Tum.Filt = Filter_by_Tum(df.PCR.Tum)
names(df.PCR.Tum.Filt)[grep("data", names(df.PCR.Tum.Filt))] = setdiff(names(df.PCR.Tum), names(df.PCR.Tum.Filt))

df.MI.vitro = df.PCR.Tum.Filt[,grep(paste(c("Exp38_","Exp1009", "Exp1011", "Exp72", "Exp1012", "BC09"), collapse="|"), names(df.PCR.Tum.Filt)), drop=F]
names(df.MI) = str_replace(names(df.MI), "_Blood", "_CTC")

#merge_sampleI create shift in vitro sample data from previous code for vitro analysis
df.PCR.Tum.Filt.I = merge_sampleI(df.PCR.Tum.Filt)
df.PCR.Tum.Filt.I = df.PCR.Tum.Filt.I[,-which(names(df.PCR.Tum.Filt.I)=="Exp72_587_cfDNA1and2")] #exclusion of cfDNA ms587

df.MI = df.PCR.Tum.Filt.I[,grep(paste(c("Exp38_","Exp1009", "Exp1011", "Exp72", "Exp1012", "BC09"), collapse="|"), names(df.PCR.Tum.Filt.I)), drop=F]
names(df.MI) = str_replace(names(df.MI), "_Blood", "_CTC")

## --- Collapse Multiple integration in Exp38 (barcodes ID : 915, 1312, 1864 added as one) ----
Exp38= df.MI[,grep("Exp38_", names(df.MI))]
bc.collapse = colSums(Exp38[c(915,1312,1864),])
Exp38[c(915,1312,1864),]=0
Exp38[915,]=bc.collapse

df.MI.temp= df.MI[,-grep("Exp38_", names(df.MI))]
df.MI = cbind(df.MI.temp, Exp38)

# Section3: Figures ----------------------------------------------------

MI.info = info[,1:6]
names(MI.info)[2] = "Mouse"
MI.info$Mouse = as.character(MI.info$Mouse)

## --- Tumours number of barcode csv export ----

names(df.MI) = str_replace(names(df.MI), "_Lung1", "_Lung")
names(df.MI) = str_replace(names(df.MI), "_Liver1", "_Liver")
names(df.MI) = str_replace(names(df.MI), "_Blood", "_CTC")

MI.tum = df.MI[,grep("_Tum", names(df.MI)), drop=F]

df.tum.bc = data.frame(nam=names(MI.tum), bc =colSums(MI.tum>0), shannon= diversity(t(MI.tum)))
df.tum.bc = cbind(df.tum.bc, str_split(df.tum.bc$nam, "_", simplify = T))
names(df.tum.bc) = c("Full_name","bc","Shannon", "Exp","Mouse","Org")

df.tum.bc= left_join(df.tum.bc,MI.info, by="Mouse")
#write.csv(df.tum.bc, file="Nb_barcode_Tum_MI.csv")

## --- Number of barcode all organs csv export  ----

df.bc = data.frame(nam=names(df.MI), bc =colSums(df.MI>0), shannon= diversity(t(df.MI)))
df.bc = cbind(df.bc, str_split(df.bc$nam, "_", simplify = T))
names(df.bc) = c("Full_name","bc","Shannon", "Exp","Mouse","Org")

unique(df.bc$Org)

Number_barcode = df.bc %>% 
  filter(Org %in% c("Tum","CTC","Liver","Lung","Liver1","Lung1","Liver2","SNLiver","SNLung","Blood","Ov","Heart","PE")) %>%
  select(Mouse,Org,bc)%>%
  pivot_wider(names_from = Mouse, values_from = bc)
Number_barcode = as.data.frame(Number_barcode)
rownames(Number_barcode) = Number_barcode$Org
Number_barcode = as.data.frame(t(Number_barcode[,-1]))
Number_barcode$Mouse = rownames(Number_barcode)

Number_barcode= left_join(Number_barcode,MI.info, by="Mouse")

write.csv(Number_barcode, file="Nb_barcode_ALL_MI.csv")

## --- Shannon diversity all organs csv export  ----

Shannon.df = df.bc %>% 
  filter(Org %in% c("Tum","CTC","Liver","Lung","Liver1","Lung1","Liver2","SNLiver","SNLung","Blood","Ov","Heart","PE")) %>%
  select(Mouse,Org,Shannon)%>%
  pivot_wider(names_from = Mouse, values_from = Shannon)
Shannon.df = as.data.frame(Shannon.df)
rownames(Shannon.df) = Shannon.df$Org
Shannon.df = as.data.frame(t(Shannon.df[,-1]))
Shannon.df$Mouse = rownames(Shannon.df)

Shannon.df= left_join(Shannon.df,MI.info, by="Mouse")

write.csv(Shannon.df, file="Nb_Shannon_ALL_MI.csv")


## --- Tumours and organs histo and bbp ----

df.bc.info= left_join(df.bc,MI.info, by="Mouse")
df.bc.info$Org = str_replace(df.bc.info$Org, "Blood", "CTC")
df.bc.info$Full_name = str_replace(df.bc.info$Full_name, "_Blood", "_CTC")

#write.csv(MDA.tum.histo$Colour, file="Colour_MDA_Histo.csv")
MDA.color = read.csv("Colour_MDA_Histo.csv")
MDA.color = MDA.color$x

## --- MDA-231 Histo ----

MDA.color = read.csv("Colour_MDA_Histo.csv")
MDA.color$x

MDA.tum = Select_and_order_Organs(df.MI, "MDA-231", "Tum")

MDA.tum.histo = Stacked_histo(MDA.tum, angle.x = 90, my_col = MDA.color$x)
MDA.tum.bbp = Bubble_plot(MDA.tum, Y = "Random", angle.x = 90)

MDA.tum.histo$Stack.plot #5x12
MDA.tum.bbp$BBP

MDA.tum.exp1= MDA.tum[,grep("Exp38_", names(MDA.tum))]
MDA.tum.exp2= MDA.tum[,grep("Exp1009_", names(MDA.tum))]
MDA.tum.exp3= MDA.tum[,grep("Exp1011_", names(MDA.tum))]

MDA.tum.exp.all = cbind(MDA.tum.exp1,MDA.tum.exp2,MDA.tum.exp3)
Bubble_plot(MDA.tum.exp.all, Y = "Random", angle.x = 90)

MDA.tum.histo1 = Stacked_histo(MDA.tum.exp1, angle.x = 90, my_col = MDA.color$x)
MDA.tum.histo2 = Stacked_histo(MDA.tum.exp2, angle.x = 90, my_col = MDA.color$x)
MDA.tum.histo3 = Stacked_histo(MDA.tum.exp3, angle.x = 90, my_col = MDA.color$x)

MDA.tum.histo1$Stack.plot #5x8
MDA.tum.histo2$Stack.plot #5x11
MDA.tum.histo3$Stack.plot #5x8


MDA.lung = Select_and_order_Organs(df.MI, "MDA-231", "Lung")

MDA.lung.histo = Stacked_histo(MDA.lung, angle.x = 90, my_col = MDA.color)
MDA.lung.bbp = Bubble_plot(MDA.lung, angle.x = 90)

MDA.lung.histo$Stack.plot #5x12

MDA.liver = Select_and_order_Organs(df.MI, "MDA-231", "Liver")

MDA.liver.histo = Stacked_histo(MDA.liver, angle.x = 90, my_col = MDA.color)
MDA.liver.bbp = Bubble_plot(MDA.liver, angle.x = 90)

MDA.liver.histo$Stack.plot #5x12

MDA.CTC = Select_and_order_Organs(df.MI, "MDA-231", "CTC")

MDA.CTC.histo = Stacked_histo(MDA.CTC, angle.x = 90, my_col = MDA.color)
MDA.CTC.bbp = Bubble_plot(MDA.CTC, angle.x = 90)

MDA.CTC.histo$Stack.plot #5x12

ggarrange(MDA.tum.histo$Stack.plot,
          MDA.lung.histo$Stack.plot,
          MDA.CTC.histo$Stack.plot,
          MDA.liver.histo$Stack.plot,
          nrow = 4 )


## --- PDX-412 Histo ----
PDX.tum = Select_and_order_Organs(df.MI, "PDX-T412", "Tum")

PDX.tum.histo = Stacked_histo(PDX.tum, angle.x = 90, my_col = PDX.color)
PDX.tum.bbp = Bubble_plot(PDX.tum, angle.x = 90)

PDX.tum.histo$Stack.plot #5x12
PDX.tum.bbp$BBP

write.csv(PDX.tum.histo$Colour, file="Colour_PDX_Histo.csv")
PDX.color = read.csv("Colour_PDX_Histo.csv")
PDX.color = PDX.color$x

PDX.lung = Select_and_order_Organs(df.MI, "PDX-T412", "Lung")
PDX.lung.histo = Stacked_histo(PDX.lung, angle.x = 90, my_col = PDX.color)
PDX.lung.bbp = Bubble_plot(PDX.lung, angle.x = 90)
PDX.lung.histo$Stack.plot #5x12
PDX.lung.bbp$BBP

PDX.CTC = Select_and_order_Organs(df.MI, "PDX-T412", "CTC")
PDX.CTC.histo = Stacked_histo(PDX.CTC, angle.x = 90, my_col = PDX.color)
PDX.CTC.bbp = Bubble_plot(PDX.CTC, angle.x = 90)
PDX.CTC.histo$Stack.plot #5x12
PDX.CTC.bbp$BBP


## --- Scatter plots -----

df.bc.info= left_join(df.bc,MI.info, by="Mouse")

df.scatter=data.frame()

Mouse.ID = unique(str_split(names(df.MI), "_", simplify = T)[,2])
Mouse.ID = Mouse.ID[!Mouse.ID %in% c("Ctrl", "Orga",   "vitroA", "vitroB", "vitroC", "vitroD")]

names(df.MI) = str_replace(names(df.MI),"_Blood","_CTC")

for (i in Mouse.ID) {
  #i="358"
  print(i)
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  TEST.name = paste("_", i, "_", sep="")
  Ms.x.Tum = df.MI[,grep(TEST.name, names(df.MI)),drop=F]
  colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
  
  Orga.MI = c("Tum", "Lung", "Liver", "CTC")
  Ms.df.temp = data.frame(ID=1:2608)
  for (j in Orga.MI) {
    #j="Liver"
    print(j)
    
    if(sum(names(Ms.x.Tum) == j)==0){Ms.df.temp[,which(j==Orga.MI)+1] = NA
    } else {Ms.df.temp[,which(j==Orga.MI)+1] = Ms.x.Tum[,j]
    }
  }
  names(Ms.df.temp) = c("ID", Orga.MI)
  columns.to.check = 1:5

  Ms.df.temp = Ms.df.temp[rowSums(Ms.df.temp[,columns.to.check[-c(1,which(is.na(colSums(Ms.df.temp))==TRUE))],drop = FALSE])>0,]
  
  Ms.df.temp$Mouse = i
  
  df.scatter = rbind(df.scatter, Ms.df.temp)


}
head(df.scatter)

df.scatter= left_join(df.scatter,MI.info, by="Mouse")
unique(df.scatter$MI)

df.scatter$facet = factor(df.scatter$MI, levels = c("IMFPc", "IMFP", "SC", "ID", "IV"))

head(df.scatter)
## --- Scatter Tum vs Lung Log -----
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(Cells~facet)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#5x14

df.scatter %>% filter(Tum>0 & Lung>0) %>% select(Mouse, Exp, Cells, facet) %>% group_by(Cells, facet) %>%
  dplyr::summarise(n_Ms = length(unique(Mouse)), n_Exp =length(unique(Exp)))

## --- Scatter Tum vs Lung NO Log -----
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter,Tum>0&Lung==0),aes(x=Tum/10000 , y=0, color=MI ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(Cells~MI)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#5x12

## --- Scatter MDA Tumour VS Lung, Liver, CTC Log scale -----

df.scatter.MDA = df.scatter %>% filter(Cells == "MDA-231")
#Tum vs Lung:
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter.MDA,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter.MDA,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#4x14
#Tum vs Liver
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA,Tum>0&Liver>0),aes(x=Tum/10000, y=Liver/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter.MDA,Tum==0&Liver>0),aes(x=0.0001, y=Liver/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter.MDA,Tum>0&Liver==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Tum>0&Liver>0),aes(x=Tum/10000, y=Liver/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Liver")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Liver (%)")#4x14
#Tum vs CTC
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter.MDA,Tum==0&CTC>0),aes(x=0.0001, y=CTC/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter.MDA,Tum>0&CTC==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs CTC")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in CTC (%)")#4x14

## --- Scatter PDX Tumour VS Lung, Liver, CTC Log scale -----

df.scatter.PDX = df.scatter %>% filter(Cells == "PDX-T412")
#Tum vs Lung
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.PDX,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter.PDX,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter.PDX,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.PDX,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#4x14

# Tum vs CTC
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.PDX,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= subset(df.scatter.PDX,Tum==0&CTC>0),aes(x=0.0001, y=CTC/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df.scatter.PDX,Tum>0&CTC==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.PDX,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs CTC")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in CTC (%)")#4x10

#Lung vs CTC
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.PDX,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.PDX,Lung==0&CTC>0),aes(x=0.001, y=CTC/10000),width=0.5,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.PDX,Lung>0&CTC==0),aes(x=Lung/10000 , y=0.05),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.PDX,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs CTC")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in CTC (%)")#4x10


## --- Scatter MDA Organs VS Lung, Liver, CTC Log scale -----

df.scatter.MDA = df.scatter %>% filter(Cells == "MDA-231")
# Lung vs CTC
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung==0&CTC>0),aes(x=0.0005, y=CTC/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung>0&CTC==0),aes(x=Lung/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs CTC")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in CTC (%)")#4x16


r2<-ddply(df.scatter.MDA,.(MI),function(x) summary(lm(x$Liver ~ x$Lung))$r.squared)
names(r2)<-c("MI","r2")

# Lung vs Liver
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  #stat_cor(data= subset(df.scatter.MDA,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=MI), alpha=0.9)+
  #geom_text(data=r2,aes(x=0.001,y=100, label = paste("R^2: ", round(r2,digit=4),sep="")),parse=T)+
  geom_point(data= subset(df.scatter.MDA,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=MI), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#4x16

ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  #stat_cor(data= subset(df.scatter.MDA,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_text(data=r2,aes(x=0.001,y=100, label = paste("R^2: ", round(r2,digit=4),sep="")),parse=T)+
  geom_point(data= subset(df.scatter.MDA,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=MI), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA,Liver==0&CTC>0),aes(x=0.0005, y=CTC/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA,Liver>0&CTC==0),aes(x=Liver/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Liver vs CTC")+
  xlab("Barcode frequencies in Liver (%)")+ylab("Barcode frequencies in CTC (%)")#4x16


length(unique(df.scatter.MDA$Mouse))
few.col= c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
few.col= brewer.pal(8, "Set2")
New.pal = colorRampPalette(few.col)
Few.COL = sample(New.pal(51), 51)
Few.COL = New.pal(51)

ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000, color=Mouse ),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005, color=Mouse ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(Cells~MI)+
  theme_classic()+
  #scale_color_manual(values = Few.COL )+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#4x18

## --- MDA Scatter plots IV and IMFP ----
df.scatter.MDA.IV = df.scatter.MDA %>% filter(MI == "IV")

scatter.IV = ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA.IV,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA.IV,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000, color=Mouse ),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA.IV,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005, color=Mouse ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA.IV,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse),method = "lm", se=F)+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #facet_grid(Cells~MI)+
  scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  scale_color_grey()+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#5x6

df.scatter.MDA.IMFP = df.scatter.MDA %>% filter(MI == "IMFP")

scatter.IMFP = ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df.scatter.MDA.IMFP,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= subset(df.scatter.MDA.IMFP,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000, color=Mouse ),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.scatter.MDA.IMFP,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005, color=Mouse ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df.scatter.MDA.IMFP,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse),method = "lm", se=F)+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #facet_grid(Cells~MI)+
  theme_classic()+
  scale_color_grey()+
  scale_x_continuous(limits = c(0,100))+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#5x6

ggarrange(scatter.IMFP, scatter.IV, nrow=1)

## --- Correlation MDA ----
cor.test(df.scatter.MDA.IV$Lung, df.scatter.MDA.IV$Liver, method = "pearson")
cor(df.scatter.MDA.IV[,2:5])

head(df.scatter.MDA)

df.scatter.MDA %>% dplyr::group_by(MI) %>% dplyr::summarize(cor(.[,2]))

test.cor = group_by(df.scatter.MDA, MI)
dplyr::summarize(test.cor, cor(Tum, Lung))
dplyr::summarize(test.cor, cor(Tum, Liver))
dplyr::summarize(test.cor, cor(Tum, CTC))
dplyr::summarize(test.cor, cor(Liver, Lung))

head(df.scatter.MDA)

df.cor = df.scatter.MDA[,c(2:5,11)]
df.cor %>% group_by(MI) %>% do(data.frame(Cor=t(cor(.[,1:4], .[,1:4]))))

cor(df.cor[1:4])

cor.IMFPc = df.cor %>% filter(MI == "IMFPc")
cor(cor.IMFPc[1:4], use="complete.obs")

cor.IMFP = df.cor %>% filter(MI == "IMFP")
cor(cor.IMFP[1:4], use="complete.obs")

cor.SC = df.cor %>% filter(MI == "SC")
cor(cor.SC[1:4], use="complete.obs")

cor.ID = df.cor %>% filter(MI == "ID")
cor(cor.ID[1:4], use="complete.obs")

cor.IV = df.cor %>% filter(MI == "IV")
cor(cor.IV[1:4])

Correlation.TUM = data.frame(IMFPc =cor(cor.IMFPc[1:4], use="complete.obs")[1,],
                             IMFP = cor(cor.IMFP[1:4], use="complete.obs")[1,], 
                             SC = cor(cor.SC[1:4], use="complete.obs")[1,],
                             ID = cor(cor.ID[1:4], use="complete.obs")[1,],
                             IV = cor(cor.IV[1:4])[1,])

Correlation.TUM

Correlation.Lung = data.frame(IMFPc =cor(cor.IMFPc[1:4], use="complete.obs")[2,],
                              IMFP = cor(cor.IMFP[1:4], use="complete.obs")[2,], 
                              SC = cor(cor.SC[1:4], use="complete.obs")[2,],
                              ID = cor(cor.ID[1:4], use="complete.obs")[2,],
                              IV = cor(cor.IV[1:4])[2,])


Correlation.Lung

Correlation.Liver = data.frame(IMFPc =cor(cor.IMFPc[1:4], use="complete.obs")[3,],
                               IMFP = cor(cor.IMFP[1:4], use="complete.obs")[3,], 
                               SC = cor(cor.SC[1:4], use="complete.obs")[3,],
                               ID = cor(cor.ID[1:4], use="complete.obs")[3,],
                               IV = cor(cor.IV[1:4])[3,])


Correlation.Liver

Correlation.CTC = data.frame(IMFPc =cor(cor.IMFPc[1:4], use="complete.obs")[4,],
                             IMFP = cor(cor.IMFP[1:4], use="complete.obs")[4,], 
                             SC = cor(cor.SC[1:4], use="complete.obs")[4,],
                             ID = cor(cor.ID[1:4], use="complete.obs")[4,],
                             IV = cor(cor.IV[1:4])[4,])

Correlation.TUM$name = rownames(Correlation.TUM)
Correlation.Lung$name = rownames(Correlation.TUM)
Correlation.CTC$name = rownames(Correlation.TUM)
Correlation.Liver$name = rownames(Correlation.TUM)

install.packages("writexl")
library(writexl)

write_xlsx(list(Tum = Correlation.TUM, Lung =Correlation.Lung, CTC=Correlation.CTC, Liver =Correlation.Liver),
           path="Correlation_MDA.xlsx")

## --- Correlation PDX ----

df.cor = df.scatter.PDX[,c(2:5,11)]
df.cor = df.cor[,-3]
head(df.cor)
cor.IMFPc = df.cor %>% filter(MI == "IMFPc")
cor(cor.IMFPc[1:3], use="complete.obs")

cor.IMFP = df.cor %>% filter(MI == "IMFP")
cor(cor.IMFP[1:3], use="complete.obs")

cor.SC = df.cor %>% filter(MI == "SC")
cor(cor.SC[1:3], use="complete.obs")

cor.ID = df.cor %>% filter(MI == "ID")
cor.ID.PDX = as.data.frame(cbind(cor(cor.ID[1:2], use="complete.obs"), CTC=c(NA,NA)))
cor.ID.PDX = rbind(cor.ID.PDX, CTC= c(NA,NA,NA))
cor.ID.PDX
as.matrix(cor.ID.PDX)[1,]

cor.IV = df.cor %>% filter(MI == "IV")
cor(cor.IV[1:3])

Correlation.TUM = data.frame(IMFPc =cor(cor.IMFPc[1:3], use="complete.obs")[1,],
                             IMFP = cor(cor.IMFP[1:3], use="complete.obs")[1,], 
                             SC = cor(cor.SC[1:3], use="complete.obs")[1,],
                             ID = as.matrix(cor.ID.PDX)[1,],
                             IV = cor(cor.IV[1:3])[1,])

Correlation.TUM

Correlation.Lung = data.frame(IMFPc =cor(cor.IMFPc[1:3], use="complete.obs")[2,],
                              IMFP = cor(cor.IMFP[1:3], use="complete.obs")[2,], 
                              SC = cor(cor.SC[1:3], use="complete.obs")[2,],
                              ID = as.matrix(cor.ID.PDX)[2,],
                              IV = cor(cor.IV[1:3])[2,])


Correlation.Lung

Correlation.CTC = data.frame(IMFPc =cor(cor.IMFPc[1:3], use="complete.obs")[3,],
                             IMFP = cor(cor.IMFP[1:3], use="complete.obs")[3,], 
                             SC = cor(cor.SC[1:3], use="complete.obs")[3,],
                             ID = as.matrix(cor.ID.PDX)[3,],
                             IV = cor(cor.IV[1:3])[3,])
Correlation.CTC

Correlation.TUM$name = rownames(Correlation.TUM)
Correlation.Lung$name = rownames(Correlation.TUM)
Correlation.CTC$name = rownames(Correlation.TUM)
Correlation.Liver$name = rownames(Correlation.TUM)

write_xlsx(list(Tum = Correlation.TUM, Lung =Correlation.Lung, CTC=Correlation.CTC, Liver =Correlation.Liver),
           path="Correlation_PDX.xlsx")

## --- Histo per mouse PDF with tumour and organs  -----
unique(df.tum.bc$Mouse)
pdf(file= "Histo_bbp_invidual_mice.pdf", height = 15, width = 15)
for (i in unique(df.tum.bc$Mouse)) {
#i=584
print(i)
TEST.name = paste("_", i, "_", sep="")
Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop=F]
id.name <- paste(str_split(names(Ms.x[1]), "_", simplify = T)[,c(1,2)], collapse = " ")
colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
if(ncol(Ms.x)==1){next}
pH = Stacked_histo(Ms.x, angle.x = 90)
pB = Bubble_plot(Ms.x,angle.x = 90, Y="No")
figure <- ggarrange(pH$Stack.plot, pB$BBP, nrow=2)
figure

Cells = df.tum.bc %>% filter(Mouse == i) %>% .$Cells
Mode_i = df.tum.bc %>% filter(Mouse == i) %>% .$MI
print(annotate_figure(figure, top = text_grob(paste(id.name,Cells,Mode_i), color = "Black", face = "bold", size = 14)))
}
dev.off()

## --- Specific example histo ----
Examples =c(365,305,324, 342, 456, 291)
i=859
Exp= paste("_", i, "_", sep = "")
str_split(Exp, "_", simplify=T)[,2]
Ms.Tum.exp = df.PCR.Tum.Filt.I[,grep(Exp, names(df.PCR.Tum.Filt.I)), drop=F]
colnames(Ms.Tum.exp) <- str_split(names(Ms.Tum.exp), "_", simplify = T)[,3]
Ms.Tum.exp = Ms.Tum.exp[,c("Tum", "Lung")]
Ms.Tum.exp$Blood = 0

p2 = Stacked_histo(Ms.Tum.exp)
p2$Stack.plot + labs(title = paste("Mouse",str_split(Exp, "_", simplify=T)[,2])) #4x4

p1 = Bubble_plot(Ms.Tum.exp)
p1$BBP +labs(title = paste("Mouse",str_split(Exp, "_", simplify=T)[,2]))

#ggarrange(p1, p2, nrow=2, ncol=1) 


## --- Bubble plot for all organs groupir for MDA ----

Bubble_plot = function(data, PCT=TRUE, Y="avg", data2, angle.x=0, my_col=FALSE){

  if(PCT==TRUE){data = data/10000} #CPM to percent
  if(Y=="avg"){data$Y = rowSums(data)/ncol(data)} else if (Y=="Random"){#Average all sample for Y order or Random order
    data$Y= sample(length(data[,1]), length(data[,1]))
  } else if(Y=="input") {data$Y = data2 #manual input of previously saved Y
  } else{data$Y = data[,1]}#Take first column as Y ref = Tum collapsed
  colnames(data) =gsub("^[^_]*_","",names(data))
  data$ID = as.numeric(rownames(data))
  
  data1 = data[rowSums(data[1:(ncol(data)-2)])>0,1:(ncol(data)-2)] #remove empty rows
  data1$Y = data[rownames(data1),"Y"]
  data1$ID = data[rownames(data1),"ID"]
  
  mG =melt(data1, id=c("ID","Y"))
  mG$ID = as.factor(mG$ID) 
  mG[mG==0] = NA 
  
  Pallett = unname(distinctColorPalette(500))
  Max.color = sample(Pallett, length(unique(data$ID)),replace = TRUE)

  if(my_col==FALSE){col.bbp = Max.color[data1$ID]} else {
    col.bbp = my_col[data1$ID]
  }
  
  Plot.bbp = ggplot()+
    #geom_point(aes(x=variable,y=Y,color=ID,size=value), data=subset(mG,variable=!"X25k_P0_1M", stat="identity"),alpha=0.5)+
    geom_quasirandom(data=mG,aes(x=variable,y=Y,color=ID,size=value),alpha=0.8, width = 0.2)+
    scale_size_continuous(range = c(1.5,15))+
    #geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)+
    scale_color_manual(values = col.bbp,guide="none")+ #breaks=length(unique(mG$ID)) removed
    xlab("")+ylab("Frequency in Tumours")+ ggtitle("")+labs(size="Frq %")+
    #scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
     #                  breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
      #                 labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    theme_minimal()+
    theme( panel.grid.major.x = element_blank() ,
           panel.grid.minor.y = element_blank(),
           axis.text.x = element_text(angle = angle.x))
  
  if(Y=="avg"){Plot.bbp = Plot.bbp +
    scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
                       breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                       labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)} else if(Y=="Random"){
      
      Plot.bbp = Plot.bbp +
        scale_y_continuous(limits = c(0,length(data[,1])))}
  Plot.bbp
    

  my_list = list(BBP = Plot.bbp, Colour = Max.color, Y.order = data$Y)
  return(my_list)
}

MDA.tum = Select_and_order_Organs(df.MI, "MDA-231", "Lung")
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "Random")
MDA.tum.bbp$BBP #6x20

Col.MDA.bbp = MDA.tum.bbp$Colour
Y.MDA.bbp = MDA.tum.bbp$Y.order

MDA.tum = Select_and_order_Organs(df.MI, "MDA-231", "Liver")
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

MDA.tum = Select_and_order_Organs(df.MI, "MDA-231", "CTC")
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x17

write.csv(data.frame(Y.Position = Y.MDA.bbp, Colours = Col.MDA.bbp ), file="Y.Position_Colours_MDA_bubbles.csv")

## --- Bubble plot for all organs groupir for PDX ----

MDA.tum = Select_and_order_Organs(df.MI, "PDX-T412", "Lung")
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "Random")
MDA.tum.bbp$BBP #6x15

Col.MDA.bbp = MDA.tum.bbp$Colour
Y.MDA.bbp = MDA.tum.bbp$Y.order

MDA.tum = Select_and_order_Organs(df.MI, "PDX-T412", "CTC")
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x8

write.csv(data.frame(Y.Position = Y.MDA.bbp, Colours = Col.MDA.bbp ), file="Y.Position_Colours_PDX_bubbles.csv")

## --- Histo per mouse PDF with tumour and organs  -----

pdf(file= "Histo_bbp_invidual_micev2.pdf", height = 15, width = 15)
Ms.list = unique(df.bc.info$Mouse)[!unique(df.bc.info$Mouse) %in% c("Ctrl","Orga","vitroA","vitroB","vitroC","vitroD")]
for (i in Ms.list) {
  #i=366
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop=F]
  id.name <- paste(str_split(names(Ms.x[1]), "_", simplify = T)[,c(1,2)], collapse = " ")
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  if(ncol(Ms.x)==1){next}
  pH = Stacked_histo(Ms.x, angle.x = 90)
  pB = Bubble_plot(Ms.x,angle.x = 90, Y="No")
  figure <- ggarrange(pH$Stack.plot, pB$BBP, nrow=2)
  figure
  
  Cells = df.bc.info %>% filter(Mouse == i) %>% .$Cells
  Mode_i = df.bc.info %>% filter(Mouse == i) %>% .$MI
  print(annotate_figure(figure, top = text_grob(paste(id.name,Cells,Mode_i), color = "Black", face = "bold", size = 14)))
}
dev.off()


## --- Specific example histo organs ----
Examples =c(365,305,324, 342, 456, 291)
i=857
Exp= paste("_", i, "_", sep = "")
str_split(Exp, "_", simplify=T)[,2]
Ms.Tum.exp = df.PCR.Tum.Filt.I[,grep(Exp, names(df.PCR.Tum.Filt.I)), drop=F]
colnames(Ms.Tum.exp) <- str_split(names(Ms.Tum.exp), "_", simplify = T)[,3]
Ms.Tum.exp$Blood = 0
Ms.Tum.exp = Ms.Tum.exp[,c("Lung", "Blood")]
names(Ms.Tum.exp)
names(Ms.Tum.exp) = c("Lung", "CTC")

p2 = Stacked_histo(Ms.Tum.exp)
p2$Stack.plot + labs(title = paste("Mouse",str_split(Exp, "_", simplify=T)[,2])) #4x4(MDA) #4x3(PDX)


p1 = Bubble_plot(Ms.Tum.exp)
p1$BBP +labs(title = paste("Mouse",str_split(Exp, "_", simplify=T)[,2]))

#ggarrange(p1, p2, nrow=2, ncol=1) 
## --- Vitro analysis -----

names(df.MI)

vitro0 = df.MI.vitro[,grep("BC09", names(df.MI.vitro)),drop=F]
vitro1 = df.MI.vitro[,grep("Amp", names(df.MI.vitro))]
vitro2 = df.MI.vitro[,grep("vitro", names(df.MI.vitro))]
vitroCtrl = df.MI.vitro[,grep(paste("Ctrl3", "Ctrl4", sep = "|"), names(df.MI.vitro))]

vitro3 = cbind(vitro0, vitro1,vitro2, vitroCtrl)

Stacked_histo(vitro3, angle.x = 90)

P1 =  vitro3 %>% select(c("BC09_25k_P0.1M_",
                          "Exp1009_cltrvitro_500k1_", "Exp1009_cltrvitro_500k2_",
                          "Exp1011_Ctrl_Ctrl3_", "Exp1011_Ctrl_Ctrl4_",
                          "Exp1009_Amp1_P1_", "Exp1009_Amp2_P1_", "Exp1009_Amp3_P1_", "Exp1009_Amp4_P1_",
                          "Exp1011_vitroA_P1_","Exp1011_vitroB_P1_","Exp1011_vitroC_P1_","Exp1011_vitroD_P1_"))
names(P1) = c("P0", "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
              "Flask1 P1", "Flask2 P1", "Flask3 P1", "Flask4 P1", "Flask5 P1", "Flask6 P1", "Flask7 P1", "Flask8 P1")

Vitro.p1 = Stacked_histo(P1, angle.x = 90, my_col = vitro.col)

Vitro.p1$Stack.plot #6x10
#vitro.col = Vitro.p1$Colour

cor(as.matrix(P1))
corrplot(as.matrix(cor(P1)),is.corr = F,cor.lim =c(0,1),
         tl.col = 'black',col = viridis(50),
         #type="upper",
         method = "color") #7x7



df.sum.vitro = data.frame(nam=names(vitro3), bc =colSums(vitro3>0), shannon= diversity(t(vitro3)))
df.sum.vitro = cbind(df.sum.vitro, str_split(df.sum.vitro$nam, "_", simplify = T))
names(df.sum.vitro) = c("Full_name","bc","Shannon", "Exp","Flask","Passage")

#write.csv(df.sum.vitro, file="Vitro_barcode_Shannon.csv")


#Bubble_plot(vitro3, Y = "Random", angle.x = 90)

a = df.sum.vitro$Flask
a = str_replace(a, "Amp1", "Flask1")
a = str_replace(a, "Amp2", "Flask2")
a = str_replace(a, "Amp3", "Flask3")
a = str_replace(a, "Amp4", "Flask4")
a = str_replace(a, "vitroA", "Flask5")
a = str_replace(a, "vitroB", "Flask6")
a = str_replace(a, "vitroC", "Flask7")
a = str_replace(a, "vitroD", "Flask8")
a = str_replace(a, "25k", "Pop")
a = str_replace(a, "cltrvitro", "Ctrl")

df.sum.vitro$new.f = a

b = df.sum.vitro$Passage
b = str_remove(b, "P")
b = str_replace(b, "0.1M", "0")
b
df.sum.vitro$new.p = b

df.sum.vitro = df.sum.vitro[,-7]

c = df.sum.vitro %>% filter(new.f != "Ctrl")
head(c)

#P8 and P13 amp1 bis to replace the original P8 and P13
c= c %>% filter(!(Flask == "Amp1" & Passage %in% c("P8", "P13")))
c$new.p = str_remove(c$new.p, "bis")

str(c)
c$new.p = as.numeric(c$new.p)

c = c %>% arrange(new.f, new.p)

c$new.name = paste(c$Exp, c$new.f, paste("P",c$new.p, sep = "") , sep="_")

vitro.order = vitro3[,c$Full_name]
names(vitro.order) =c$new.name
names(vitro.order)

#Move last column to first position:
vitro.order = vitro.order %>% select("BC09_Pop_P0", everything())

vitro.SH = Stacked_histo(vitro.order, angle.x = 90, my_col = vitro.col)
vitro.SH$Stack.plot #6x25 and #10x25

vitro.bbp = Bubble_plot(vitro.order, Y = "Random", angle.x = 90)
vitro.bbp$BBP


P10 = vitro.order[,grep("P10", names(vitro.order))]
MDA.tum = Select_and_order_Organs(df.PCR.Tum.Filt, "MDA-231", "Tum")
MDA.tum.vitro = MDA.tum[,grep(paste("Exp1009", "Exp1011", sep="|"), names(MDA.tum))]

MDA.tum.vitro.P10 = cbind(MDA.tum.vitro, P10)

cor(as.matrix(log(MDA.tum.vitro.P10+1)))
corrplot(as.matrix(cor(log(MDA.tum.vitro.P10+1))),is.corr = F,cor.lim =c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color") #7x7
Bubble_plot(MDA.tum.vitro.P10, angle.x = 90) #6x15

MDA.tum.vitro.P10$Pop_P0 = vitro.order$BC09_Pop_P0
MDA.tum.vitro.P10 = MDA.tum.vitro.P10 %>% select("Pop_P0", everything())


names(P1)
ctrl.avg = P1[,c("Ctrl1","Ctrl2","Ctrl3","Ctrl4")]
ctrl.avg$Y = rowSums(ctrl.avg)/ncol(ctrl.avg)
ctrl.avg$ID = row.names(ctrl.avg)
ctrl.avg = ctrl.avg %>% arrange(desc(Y))

ctrl.avg$ID[1:10]

#Function to have barcode ID of the top X barcode in avg in population
TOP_X_bc = function(data, X_bc=10){
  data$Y = rowSums(data)/ncol(data)
  data$ID = row.names(data)
  data = data %>% arrange(desc(Y))
  return(data$ID[1:X_bc])
}

IMFP.Ms.vitro = df.bc.info %>% filter(MI=="IMFP" & Org =="Tum") %>% filter(Exp.x %in% c("Exp1009", "Exp1011")) %>% .$Mouse

MDA.tum.IMFP = MDA.tum[, grep(paste(IMFP.Ms.vitro, collapse = "|"), names(MDA.tum))]
names(MDA.tum.IMFP) = c(1:10)

P5 = vitro.order[,grep("P5", names(vitro.order))]
mean(colSums(P5[ctrl.avg$ID[1:10],]))/10000

P10 = vitro.order[,grep("P10", names(vitro.order))]
mean(colSums(P10[ctrl.avg$ID[1:10],]))/10000

P15 = vitro.order[,grep("P15", names(vitro.order))]
mean(colSums(P15[ctrl.avg$ID[1:10],]))/10000

P25 = vitro.order[,grep("P25", names(vitro.order))]
mean(colSums(P25[ctrl.avg$ID[1:10],]))/10000


vitro.all.P.tum = cbind(vitro.order, MDA.tum)

df.vitro.X_bc = rbind(data.frame(x= "P1", y = colSums(P1[TOP_X_bc(ctrl.avg),6:13])/10000),
                      data.frame(x= "P5", y = colSums(P5[TOP_X_bc(P5),])/10000),
                      data.frame(x= "P10", y = colSums(P10[TOP_X_bc(P10),])/10000),
                      data.frame(x= "P15", y = colSums(P15[TOP_X_bc(P15),])/10000),
                      data.frame(x= "P25", y = colSums(P25[TOP_X_bc(P25),])/10000))
df.vitro.X_bc$x = factor(df.vitro.X_bc$x, levels = c("P1", "P5", "P10", "P15", "P25"))

df.tum.X_bc = rbind(data.frame(x= "P1", y = colSums(MDA.tum.IMFP[TOP_X_bc(ctrl.avg),])/10000),
                    data.frame(x= "P5", y = colSums(MDA.tum.IMFP[TOP_X_bc(P5),])/10000),
                    data.frame(x= "P10", y = colSums(MDA.tum.IMFP[TOP_X_bc(P10),])/10000),
                    data.frame(x= "P15", y = colSums(MDA.tum.IMFP[TOP_X_bc(P15),])/10000),
                    data.frame(x= "P25", y = colSums(MDA.tum.IMFP[TOP_X_bc(P25),])/10000))
df.tum.X_bc$x = factor(df.tum.X_bc$x, levels = c("P1", "P5", "P10", "P15", "P25")) 


X_vitro = ggplot(df.vitro.X_bc)+
  geom_bar(aes(x=x ,y=y), stat= "summary", fun="mean")+
  geom_point(aes(x=x ,y=y))+
  coord_flip()+
  labs(x="", y="Frequency in vitro (%)")

X_vivo =ggplot(df.tum.X_bc)+
  geom_bar(aes(x=x ,y=y), stat= "summary", fun="mean")+
  geom_point(aes(x=x ,y=y))+
  scale_y_reverse(limits=c(100,0))+
  coord_flip()+
  labs(x="", y="Frequency in tumors (%)")
  #theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

annotate_figure(ggarrange(X_vivo, X_vitro, nrow = 1),
                top = text_grob("Top 10 barcodes in vitro contribution", face = "bold", size = 12))
