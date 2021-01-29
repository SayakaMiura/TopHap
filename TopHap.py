<<<<<<< HEAD
import Functions3
import os
import glob

CutFre=0.05 #vf
THCutFre=0.05 #hf
BooC=100 #bootstrap replicates
BooFillPer=0.02 #hf for bootstrap replicates
Exclude='Excluded.txt'
Rpath='Rscript'#'Rscript'#'\"C:\\Program Files\\R\\R-4.0.3\\bin\\Rscript\"'
Fol='Alignment\\'
OutGfas=Fol+'Outgroup.fasta'
PositionForce='MinorVariants.txt'
FasTaLs=[Fol+'Note.txt',Fol+'Note-Continent.txt',Fol+'Note-Country.txt','AdditionalData.txt']
SC='n'
Miss='0'
#Fol='C:\\Users\\tuf78332\\Desktop\\TopHap\\TallG5clone50\\'
#OutGfas=Fol+'Out.fasta'
#FasTaLs=[Fol+'Note.txt']


PositionForce=open(PositionForce,'r').readlines()[1:]
PositionForceLs=[]#from 1
PosExcludeLs=[]#from 1
Exclu=open(Exclude,'r').readlines()
for i in Exclu:
   PosExcludeLs.append(int(i)) 
for i in PositionForce:
    #if i.split('\t')[1].strip()=='Exclude': PosExcludeLs.append(int(i.split('\t')[0]))
    #else: 
       PositionForceLs.append(int(i.split('\t')[0]))
print ('forced positions',PositionForceLs)	
print ('excluded positions number',len(PosExcludeLs))	
#open('A','r').readlines()
PosDic={} #pos from 1
PosMajorNuc={}
TempIDls=[]
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:		
            Fas=Fol+Name+Fas0.replace(' ','_')+'.fasta'	
            TempIDls.append(Fol+Name+Fas0.replace(' ','_'))			
            print ('make PosInfo file',Fas)
            if os.path.exists(Fas[:-6]+'_PosInfo.txt')!=True: Functions3.GetMajorFre(Fas,CutFre) #pos from 1
            PosTa=open(Fas[:-6]+'_PosInfo.txt','r').readlines()[1:]
            for i in PosTa:
                i=i.strip().split('\t')
                if i[1]!='NA' and i[4]!='NA':
                    if float(i[3])>CutFre and (1-CutFre)>float(i[3]) and float(i[6])>CutFre and (1-CutFre)>float(i[6]) and PosExcludeLs.count(int(i[0]))==0:
                        PosDic[int(i[0])]=PosDic.get(int(i[0]),[])+[Fas0] 
                if i[1]!='NA': 
                 if (1-CutFre)<=float(i[3]) and PosExcludeLs.count(int(i[0]))==0:
                   PosMajorNuc[int(i[0])]=PosMajorNuc.get(int(i[0]),[])+[i[1]]

for Pos in PosMajorNuc:
    MnucLs=list(set(PosMajorNuc[Pos]))
    #if Pos==23063: print (MnucLs)    	
    #Find=PosDic.get(Pos,'NA')	
    if len(MnucLs)>1:
       PosDic[Pos]=PosDic.get(Pos,[])+['fixed']
       print ('more than one major',Pos,MnucLs)	   
#open('A','r').readlines()	
for i in PositionForceLs:
   if PosExcludeLs.count(i)==0:
    PosDic[i]=PosDic.get(i,[])+['force']
   else:
     print ('forced position is excluded, because no information in the 68K database',i)   
print ('common variant Count',len(PosDic))						
PosLs=list(PosDic.keys())
PosLs.sort()
out='Pos from 1\tData Ls\n'
for Pos in PosLs:
    out+=str(Pos)+'\t'+';'.join(PosDic[Pos])+'\n'
Functions3.GetOut('CommonVarPosLs.txt',out)	

TopHap2Data={}
#SubData2TopHapSeq2HapID={}
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
      #  TopHapSeq2HapID={}		
        if Fas0!='NA' and Fas0.find('-not used')==-1:		
            Fas=Fol+Name+Fas0.replace(' ','_')+'.fasta'	
            print ('make Hap fas',Fas)
            if os.path.exists(Fas[:-6]+'_Hap.fasta')!=True: Functions3.MakeHapFas(Fas,PosLs)
            if os.path.exists(Fas[:-6]+'_Hap_hapID.fasta')!=True: 
                  if SC!='y': Functions3.TopHap_haplotype_filtering(THCutFre,Fas[:-6]+'_Hap.fasta')	
                 # else: os.system(CallPy+' TopHap_haplotype_filtering_SC.py '+str(THCutFre)+' '+Fas[:-6]+'_Hap.fasta '+Miss)					  
            TopHapFasLs=glob.glob(Fas[:-6]+'_Hap_hapID_Filter*.fasta')
            if len(TopHapFasLs)!=1:
               print (TopHapFasLs)
               open('A','r').readlines()
            else:
               HapLs,Hap2Seq=Functions3.ReadFasSeq(TopHapFasLs[0])
               for Hap in HapLs:			   
                  Seq=Hap2Seq[Hap]			   
                  TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]
                 # TopHapSeq2HapID[Seq]=TopHapSeq2HapID.get(Seq,[])+[Hap]
             #  SubData2TopHapSeq2HapID[Fas0]=TopHapSeq2HapID				  
print ('tophap c',len(TopHap2Data))
TopID=1
ID2TopHap={}
outF=''
ID2fasID={}
out='HapID\tData list\n'
TopHap2ID={}	
for Hap in TopHap2Data:
    ID='All'+str(TopID)
    ID2fasID[ID]=ID+'_'+TopHap2Data[Hap][0].replace(' ','_')+'_'+str(len(TopHap2Data[Hap]))	
    outF+='>'+ID+'_'+TopHap2Data[Hap][0].replace(' ','_')+'_'+str(len(TopHap2Data[Hap]))+'\n'+Hap+'\n'
    out+=ID+'\t'+';'.join(TopHap2Data[Hap])+'\n'
    TopHap2ID[Hap]=ID	
    ID2TopHap[ID]=Hap	
    TopID+=1
Ols,O2Seq=Functions3.ReadFasSeq(OutGfas)
outG=''
for O in Ols:
    Seq=O2Seq[O]
    Hap=''	
    for Pos in PosLs:
       Hap+=Seq[Pos-1]
    outG +=O+'\n'+Hap+'\n'
outF=outG+outF	
Functions3.GetOut(Fol+'OutG.fasta',outG)
Functions3.GetOut('TopHap.fasta',outF)
Functions3.GetOut('TopHap.txt',out)	
         #   open('A','r').readlines()			
print ('annotation')
out='Group\tSeq\tTopHapID\n'	
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:				
            Fas=Fol+Name+Fas0.replace(' ','_')+'_Hap.fasta'	
            SeqLs,Seq2Hap=Functions3.ReadFasSeq(Fas)
            for Seq in SeqLs:
                 Hap=Seq2Hap[Seq]
                 Tid=TopHap2ID.get(Hap,'NA')
                 out+=Fas0+'\t'+Seq.replace('>','')+'\t'+Tid+'\n'
Functions3.GetOut('TopHap_anno.txt',out)

Data2Count={}
Data2Tot={}	
TopHapOr=list(ID2TopHap.keys())
print ('count')
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:
            Fas=Fol+Name+Fas0.replace(' ','_')+'_Hap.fasta'	
            StLs,St2Hap=Functions3.ReadFasSeq(Fas)
            TopH2Count={}			
            for St in StLs:
                Hap=St2Hap[St]			
                ID=TopHap2ID.get(Hap,'NA')
                TopH2Count[ID]=TopH2Count.get(ID,0)+1
            Data2Count[Fas0]=TopH2Count
            Data2Tot[Fas0]=len(StLs)			
out='TopHap\t'
out1='TopHap\t'
Datals=list(Data2Count.keys())
out+='\t'.join(Datals)+'\n'
out1+='\t'.join(Datals)+'\n'
for ID in TopHapOr:
    out+=ID2fasID[ID]
    out1+=ID2fasID[ID]	
    for Data in Datals:
        CountDic=Data2Count[Data]
        C=CountDic.get(ID,0)
        out+='\t'+str(C)
        Frac=1.0*C/Data2Tot[Data]		
        out1+='\t'+str(Frac)		
    out+='\n'
    out1+='\n'	
out+='Total'	
for Data in Datals:
    out+='\t'+str(Data2Tot[Data])	
Functions3.GetOut('TopHap_count.txt',out)
Functions3.GetOut('TopHap_count_freq.txt',out1)	
        	
print ('infer MP')
os.system('megacc -a infer_MP_nucleotide.mao -d ' + 'TopHap.fasta' +' -o ' + 'TopHap_allMP.nwk')
AncFileLs=glob.glob('TopHap_allMP_ancestral_states_*.txt')
for AncF in AncFileLs:
               os.remove(AncF)	  
os.remove('TopHap_allMP_summary.txt')			   
print ('do bootstrap')
Functions3.BootstrapSeq(BooC,BooFillPer,'TopHap_allMP.nwk','TopHap.fasta',Rpath)
os.remove('AllResampTrees.nwk')	
os.remove('Rplots.pdf')	
os.remove('TopHap_count_freq.txt')
os.remove('TopHap.fasta')
Ls=glob.glob('TopHap_allMP_*_prune.nwk')
for i in Ls:
    os.remove(i)
print (TempIDls)
for i in TempIDls:
    Ls=glob.glob(i+'_Hap_hapID*')
    os.remove(i+'_Hap.fasta')
    for i in Ls:
      os.remove(i)	
=======
import Functions3
import os
import glob

CutFre=0.05 #vf
THCutFre=0.05 #hf
BooC=100 #bootstrap replicates
BooFillPer=0.02 #hf for bootstrap replicates
Exclude='Excluded.txt'
Rpath='Rscript'#'Rscript'#'\"C:\\Program Files\\R\\R-4.0.3\\bin\\Rscript\"'
Fol='Alignment\\'
OutGfas=Fol+'Outgroup.fasta'
PositionForce='MinorVariants.txt'
FasTaLs=[Fol+'Note.txt',Fol+'Note-Continent.txt',Fol+'Note-Country.txt','AdditionalData.txt']
SC='n'
Miss='0'
#Fol='C:\\Users\\tuf78332\\Desktop\\TopHap\\TallG5clone50\\'
#OutGfas=Fol+'Out.fasta'
#FasTaLs=[Fol+'Note.txt']


PositionForce=open(PositionForce,'r').readlines()[1:]
PositionForceLs=[]#from 1
PosExcludeLs=[]#from 1
Exclu=open(Exclude,'r').readlines()
for i in Exclu:
   PosExcludeLs.append(int(i)) 
for i in PositionForce:
    #if i.split('\t')[1].strip()=='Exclude': PosExcludeLs.append(int(i.split('\t')[0]))
    #else: 
       PositionForceLs.append(int(i.split('\t')[0]))
print ('forced positions',PositionForceLs)	
print ('excluded positions number',len(PosExcludeLs))	
#open('A','r').readlines()
PosDic={} #pos from 1
PosMajorNuc={}
TempIDls=[]
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:		
            Fas=Fol+Name+Fas0.replace(' ','_')+'.fasta'	
            TempIDls.append(Fol+Name+Fas0.replace(' ','_'))			
            print ('make PosInfo file',Fas)
            if os.path.exists(Fas[:-6]+'_PosInfo.txt')!=True: Functions3.GetMajorFre(Fas,CutFre) #pos from 1
            PosTa=open(Fas[:-6]+'_PosInfo.txt','r').readlines()[1:]
            for i in PosTa:
                i=i.strip().split('\t')
                if i[1]!='NA' and i[4]!='NA':
                    if float(i[3])>CutFre and (1-CutFre)>float(i[3]) and float(i[6])>CutFre and (1-CutFre)>float(i[6]) and PosExcludeLs.count(int(i[0]))==0:
                        PosDic[int(i[0])]=PosDic.get(int(i[0]),[])+[Fas0] 
                if i[1]!='NA': 
                 if (1-CutFre)<=float(i[3]) and PosExcludeLs.count(int(i[0]))==0:
                   PosMajorNuc[int(i[0])]=PosMajorNuc.get(int(i[0]),[])+[i[1]]

for Pos in PosMajorNuc:
    MnucLs=list(set(PosMajorNuc[Pos]))
    #if Pos==23063: print (MnucLs)    	
    #Find=PosDic.get(Pos,'NA')	
    if len(MnucLs)>1:
       PosDic[Pos]=PosDic.get(Pos,[])+['fixed']
       print ('more than one major',Pos,MnucLs)	   
#open('A','r').readlines()	
for i in PositionForceLs:
   if PosExcludeLs.count(i)==0:
    PosDic[i]=PosDic.get(i,[])+['force']
   else:
     print ('forced position is excluded, because no information in the 68K database',i)   
print ('common variant Count',len(PosDic))						
PosLs=list(PosDic.keys())
PosLs.sort()
out='Pos from 1\tData Ls\n'
for Pos in PosLs:
    out+=str(Pos)+'\t'+';'.join(PosDic[Pos])+'\n'
Functions3.GetOut('CommonVarPosLs.txt',out)	

TopHap2Data={}
#SubData2TopHapSeq2HapID={}
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
      #  TopHapSeq2HapID={}		
        if Fas0!='NA' and Fas0.find('-not used')==-1:		
            Fas=Fol+Name+Fas0.replace(' ','_')+'.fasta'	
            print ('make Hap fas',Fas)
            if os.path.exists(Fas[:-6]+'_Hap.fasta')!=True: Functions3.MakeHapFas(Fas,PosLs)
            if os.path.exists(Fas[:-6]+'_Hap_hapID.fasta')!=True: 
                  if SC!='y': Functions3.TopHap_haplotype_filtering(THCutFre,Fas[:-6]+'_Hap.fasta')	
                 # else: os.system(CallPy+' TopHap_haplotype_filtering_SC.py '+str(THCutFre)+' '+Fas[:-6]+'_Hap.fasta '+Miss)					  
            TopHapFasLs=glob.glob(Fas[:-6]+'_Hap_hapID_Filter*.fasta')
            if len(TopHapFasLs)!=1:
               print (TopHapFasLs)
               open('A','r').readlines()
            else:
               HapLs,Hap2Seq=Functions3.ReadFasSeq(TopHapFasLs[0])
               for Hap in HapLs:			   
                  Seq=Hap2Seq[Hap]			   
                  TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]
                 # TopHapSeq2HapID[Seq]=TopHapSeq2HapID.get(Seq,[])+[Hap]
             #  SubData2TopHapSeq2HapID[Fas0]=TopHapSeq2HapID				  
print ('tophap c',len(TopHap2Data))
TopID=1
ID2TopHap={}
outF=''
ID2fasID={}
out='HapID\tData list\n'
TopHap2ID={}	
for Hap in TopHap2Data:
    ID='All'+str(TopID)
    ID2fasID[ID]=ID+'_'+TopHap2Data[Hap][0].replace(' ','_')+'_'+str(len(TopHap2Data[Hap]))	
    outF+='>'+ID+'_'+TopHap2Data[Hap][0].replace(' ','_')+'_'+str(len(TopHap2Data[Hap]))+'\n'+Hap+'\n'
    out+=ID+'\t'+';'.join(TopHap2Data[Hap])+'\n'
    TopHap2ID[Hap]=ID	
    ID2TopHap[ID]=Hap	
    TopID+=1
Ols,O2Seq=Functions3.ReadFasSeq(OutGfas)
outG=''
for O in Ols:
    Seq=O2Seq[O]
    Hap=''	
    for Pos in PosLs:
       Hap+=Seq[Pos-1]
    outG +=O+'\n'+Hap+'\n'
outF=outG+outF	
Functions3.GetOut(Fol+'OutG.fasta',outG)
Functions3.GetOut('TopHap.fasta',outF)
Functions3.GetOut('TopHap.txt',out)	
         #   open('A','r').readlines()			
print ('annotation')
out='Group\tSeq\tTopHapID\n'	
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:				
            Fas=Fol+Name+Fas0.replace(' ','_')+'_Hap.fasta'	
            SeqLs,Seq2Hap=Functions3.ReadFasSeq(Fas)
            for Seq in SeqLs:
                 Hap=Seq2Hap[Seq]
                 Tid=TopHap2ID.get(Hap,'NA')
                 out+=Fas0+'\t'+Seq.replace('>','')+'\t'+Tid+'\n'
Functions3.GetOut('TopHap_anno.txt',out)

Data2Count={}
Data2Tot={}	
TopHapOr=list(ID2TopHap.keys())
print ('count')
for FasLs in FasTaLs:
    print (FasLs)
    if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
    else: Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
        Fas0=i.split('\t')[0]
        if Fas0!='NA' and Fas0.find('-not used')==-1:
            Fas=Fol+Name+Fas0.replace(' ','_')+'_Hap.fasta'	
            StLs,St2Hap=Functions3.ReadFasSeq(Fas)
            TopH2Count={}			
            for St in StLs:
                Hap=St2Hap[St]			
                ID=TopHap2ID.get(Hap,'NA')
                TopH2Count[ID]=TopH2Count.get(ID,0)+1
            Data2Count[Fas0]=TopH2Count
            Data2Tot[Fas0]=len(StLs)			
out='TopHap\t'
out1='TopHap\t'
Datals=list(Data2Count.keys())
out+='\t'.join(Datals)+'\n'
out1+='\t'.join(Datals)+'\n'
for ID in TopHapOr:
    out+=ID2fasID[ID]
    out1+=ID2fasID[ID]	
    for Data in Datals:
        CountDic=Data2Count[Data]
        C=CountDic.get(ID,0)
        out+='\t'+str(C)
        Frac=1.0*C/Data2Tot[Data]		
        out1+='\t'+str(Frac)		
    out+='\n'
    out1+='\n'	
out+='Total'	
for Data in Datals:
    out+='\t'+str(Data2Tot[Data])	
Functions3.GetOut('TopHap_count.txt',out)
Functions3.GetOut('TopHap_count_freq.txt',out1)	
        	
print ('infer MP')
os.system('megacc -a infer_MP_nucleotide.mao -d ' + 'TopHap.fasta' +' -o ' + 'TopHap_allMP.nwk')
AncFileLs=glob.glob('TopHap_allMP_ancestral_states_*.txt')
for AncF in AncFileLs:
               os.remove(AncF)	  
os.remove('TopHap_allMP_summary.txt')			   
print ('do bootstrap')
Functions3.BootstrapSeq(BooC,BooFillPer,'TopHap_allMP.nwk','TopHap.fasta',Rpath)
os.remove('AllResampTrees.nwk')	
os.remove('Rplots.pdf')	
os.remove('TopHap_count_freq.txt')
os.remove('TopHap.fasta')
Ls=glob.glob('TopHap_allMP_*_prune.nwk')
for i in Ls:
    os.remove(i)
print (TempIDls)
for i in TempIDls:
    Ls=glob.glob(i+'_Hap_hapID*')
    os.remove(i+'_Hap.fasta')
    for i in Ls:
      os.remove(i)	
>>>>>>> f3cf616 (Initial commit)
	