#---------------Native--------------------------------#
rm(list=ls(all=TRUE))
library("readxl")
getwd()
setwd("Path")
Data<-read_excel("MSdata.xlsx")
class(Data)
Data<-as.data.frame(Data)
Data_out<-c("File","m/z", "MaxQunatMass","Entry","Composition","Neu5Gc","Neu5Ac","HexNAc","Hex","Fuc","Red_HexNAc","deviation") 

  Neu5Gc<-307.09033
  Neu5Ac<-291.09541
  HexNAc<-203.07937
  Hex<-162.05282
  Fuc<-146.05791

  for (i in 1:nrow(Data)) {
  print(i)
  q<-Data$Mass[i]
  s<-q%/%Neu5Gc
  t<-q%/%Neu5Ac
  j<-q%/%HexNAc
  m<-q%/%Hex
  n<-q%/%Fuc
  
      for(Neu5Gc_n in 0:s){
       for(Neu5Ac_n in 0:t){
        for( HexNAc_n in 0:j){
          for( Hex_n in 0:m){
            for( Fuc_n in 0:n){
                total_1<-sum(325.1009*Neu5Gc_n,309.10598*Neu5Ac_n,221.08994*HexNAc_n,180.06339*Hex_n,164.06848*Fuc_n,223.10559)
                total_2<-total_1-sum(Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n)*18.01057
                deviation<-abs(total_2-q)/total_2*1000000

                if(deviation<5)
                {  
                   Entry<-paste ("N",as.character(Neu5Gc_n),as.character(Neu5Ac_n),as.character(HexNAc_n),as.character(Hex_n),as.character(Fuc_n),"1", sep = "")
                   Composition<-paste("Neu5Gc",as.character(Neu5Gc_n),"Neu5Ac", as.character(Neu5Ac_n),"HexNAc",as.character(HexNAc_n),"Hex",as.character(Hex_n),"Fuc",as.character(Fuc_n),"Red-HexNAc","1", sep="")
                   MS1<-c(Data$`Raw file`[i], Data$`m/z`[i], Data$Mass[i],Entry,Composition,Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n,"1",deviation) 
                   Data_out<-rbind(Data_out, MS1)
                }  else {}
                }
  }
  }
  }
  }
  }
  
 colnames(Data_out)<-Data_out[1,]
 Data_out<-Data_out[-1,]

write.table(Data_out,row.names=F,col.names=T,file="Data_out.csv",sep=",")



#---------------Acetohydrazide(Ah)-based derivatization--------------------------------#
rm(list=ls(all=TRUE))
library("readxl")
getwd()
setwd("Path")
Data<-read_excel("MSdata.xlsx")
class(Data)
Data<-as.data.frame(Data)
Data_out<-c("File","m/z", "MaxQunatMass","Entry","Composition","Neu5Gc","Neu5Ac","HexNAc","Hex","Fuc","Red_HexNAc","deviation") 

Neu5Gc<-363.12778
Neu5Ac<-347.13286
HexNAc<-203.07937
Hex<-162.05282
Fuc<-146.05791

for (i in 1:nrow(Data)) {
  print(i)
  q<-Data$Mass[i]
  s<-q%/%Neu5Gc
  t<-q%/%Neu5Ac
  j<-q%/%HexNAc
  m<-q%/%Hex
  n<-q%/%Fuc
  
  for(Neu5Gc_n in 0:s){
    for(Neu5Ac_n in 0:t){
      for( HexNAc_n in 0:j){
        for( Hex_n in 0:m){
          for( Fuc_n in 0:n){
            total_1<-sum(381.13835*Neu5Gc_n,365.14343*Neu5Ac_n,221.08994*HexNAc_n,180.06339*Hex_n,164.06848*Fuc_n,223.10559)
            total_2<-total_1-sum(Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n)*18.01057
            deviation<-abs(total_2-q)/total_2*1000000
            
            if(deviation<5)
            {  
              Entry<-paste ("N",as.character(Neu5Gc_n),as.character(Neu5Ac_n),as.character(HexNAc_n),as.character(Hex_n),as.character(Fuc_n),"1", sep = "")
              Composition<-paste("Neu5Gc",as.character(Neu5Gc_n),"Neu5Ac", as.character(Neu5Ac_n),"HexNAc",as.character(HexNAc_n),"Hex",as.character(Hex_n),"Fuc",as.character(Fuc_n),"Red-HexNAc","1", sep="")
              MS1<-c(Data$`Raw file`[i], Data$`m/z`[i], Data$Mass[i],Entry,Composition,Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n,"1",deviation) 
              Data_out<-rbind(Data_out, MS1)
            }  else {}
          }
        }
      }
    }
  }
}

colnames(Data_out)<-Data_out[1,]
Data_out<-Data_out[-1,]

write.table(Data_out,row.names=F,col.names=T,file="Data_out.csv",sep=",")



#---------------Permethylation--------------------------------#
rm(list=ls(all=TRUE))
library("readxl")
getwd()
setwd("Path")
Data<-read_excel("MSdata.xlsx")
class(Data)
Data<-as.data.frame(Data)
Data_out<-c("File","m/z", "MaxQunatMass","Entry","Composition","Neu5Gc","Neu5Ac","HexNAc","Hex","Fuc","Red_HexNAc","deviation") 

Neu5Gc<-391.18423
Neu5Ac<-361.17366
HexNAc<-245.12632
Hex<-204.09977
Fuc<-174.08921

for (i in 1:nrow(Data)) {
  print(i)
  q<-Data$Mass[i]
  s<-q%/%Neu5Gc
  t<-q%/%Neu5Ac
  j<-q%/%HexNAc
  m<-q%/%Hex
  n<-q%/%Fuc
  
  for(Neu5Gc_n in 0:s){
    for(Neu5Ac_n in 0:t){
      for( HexNAc_n in 0:j){
        for( Hex_n in 0:m){
          for( Fuc_n in 0:n){
            total_1<-sum(437.2261*Neu5Gc_n,407.21553*Neu5Ac_n,291.16819*HexNAc_n,250.14164*Hex_n,220.13108*Fuc_n,307.19949)
            total_2<-total_1-sum(Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n)*46.04187
            deviation<-abs(total_2-q)/total_2*1000000
            
            if(deviation<5)
            {  
              Entry<-paste ("N",as.character(Neu5Gc_n),as.character(Neu5Ac_n),as.character(HexNAc_n),as.character(Hex_n),as.character(Fuc_n),"1", sep = "")
              Composition<-paste("Neu5Gc",as.character(Neu5Gc_n),"Neu5Ac", as.character(Neu5Ac_n),"HexNAc",as.character(HexNAc_n),"Hex",as.character(Hex_n),"Fuc",as.character(Fuc_n),"Red-HexNAc","1", sep="")
              MS1<-c(Data$`Raw file`[i], Data$`m/z`[i], Data$Mass[i],Entry,Composition,Neu5Gc_n,Neu5Ac_n,HexNAc_n,Hex_n,Fuc_n,"1",deviation) 
              Data_out<-rbind(Data_out, MS1)
            }  else {}
          }
        }
      }
    }
  }
}

colnames(Data_out)<-Data_out[1,]
Data_out<-Data_out[-1,]

write.table(Data_out,row.names=F,col.names=T,file="Data_out.csv",sep=",")