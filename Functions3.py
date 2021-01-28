import os
from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import glob
def ReadFasSeq(pos72):
 StLs=[]
 St2Seq={}
 pos72=open(pos72,'r').readlines()
 for i in pos72:
    if i[0]=='>': 
       ID=i.strip()
       StLs.append(ID)
       St2Seq[ID]=''
    else: St2Seq[ID]+=i.strip()
 return StLs,St2Seq 	
def get_refHap(AllFas,AllAnnoTa,ExtPosLsTa):
    print ('get ref align')
    AllAnnoTa=open(AllAnnoTa,'r').readlines()[1:]
    Ref2Seq={}
    for i in AllAnnoTa:
      if i.strip()!='':	
        i=i.split()	
        if i[2].strip()=='0' and i[3].strip()=='0':
           Ref2Seq[i[1]]=i[0]
    #print (Ref2Seq)
    Seq2Ref=InvertDic(Ref2Seq)	
    ExtPosLsTa=open(ExtPosLsTa,'r').readlines()[1:]
    PosDic={} #from 0
    PosLs=[]	
    C=0
    for i in ExtPosLsTa:
        if i.strip()!='':
            PosDic[C]=int(i.split('\t')[0])-1
            PosLs.append(int(i.split('\t')[0])-1)			
            C+=1
    print (PosDic,PosLs)			
    ExtractPos(PosLs,AllFas,list(Seq2Ref.keys()),Seq2Ref,AllFas[:-6]+'_refHap.fasta')	#count from 0
    	

    return PosDic
def GetHap(Seq,PosLs): #from 0
    Hap=''
    for Pos in PosLs:
       if len(Seq)<=Pos: 
          Hap+='?'
          open('A','r').readlines()		  
       else: Hap+=Seq[Pos]
    return Hap
def GetMostSim(Hap,RefHap):
    DifC2HapLs={}
    Len=len(Hap)	
    for Ref in RefHap:
        RefSeq=RefHap[Ref]	
        if len(RefSeq)!=Len: 
            print (Hap,RefSeq,Len,len(RefSeq))
            open('a','r').readlines()
        else:
           DC=CountDifNum(Hap,RefSeq)
           DifC2HapLs[DC]=DifC2HapLs.get(DC,[])+[Ref.replace('>','')]
    DifCLs=list(DifC2HapLs.keys())
    DifCLs.sort()
    DifC=DifCLs[0]	
    Gid=';'.join(DifC2HapLs[DifC])	

    return Gid,DifC	
def Map_toRefHap(ID2Seq,PosLs,RefHap,ID2Info,OutFile): #from 0
    out='Strain\tLocation\tSampTime\tGroup\tDifC\n'
    St2Gid={}	
    for ID in ID2Seq:
       Seq=ID2Seq[ID]
       Hap=GetHap(Seq,PosLs) #from 0	
       Gid,DifC=GetMostSim(Hap,RefHap)	
       out+=ID+'\t'+ID2Info[ID]+'\t'+Gid+'\t'+str(DifC)+'\n'
       St2Gid[ID]=Gid+'\t'+str(DifC)	   
    GetOut(OutFile,out)	 
    return St2Gid	
	
def ExtractPos(PosLs,OriFas,SeqLs_toExtract,Seq2Ref_newID,OutFile): #from 0
    print ('reading input fas')
    StLs,St2Seq=ReadFasSeq(OriFas)
    out=''	 
    print ('making hap')		
    for St in SeqLs_toExtract:

        Seq=St2Seq['>'+St]
        Hap=''
        for Pos in PosLs:
           Hap+=Seq[Pos]
        if 	Seq2Ref_newID!={}: ID=Seq2Ref_newID[St]
        else: ID=St	
        out+='>'+ID+'\n'+Hap+'\n'		
    GetOut(OutFile,out) 
def GetMajorFre(Fas,CutFre): #pos from 1
   StLs,St2Seq=ReadFasSeq(Fas)
   SeqC=len(StLs)   
   Len=len((St2Seq[StLs[0]]))
   print ('seq count',SeqC,'nuc count',Len)
   c=0
   ATGC=['A','T','G','C'] 
   out='Pos from1\tMajorNuc\tMajorCou\tMajorFre\tMinorNuc\tMinorCou\tMinorFre\n'   
   while c<Len:
      Nuc2Cou={}
      for St in StLs:
         Nuc=St2Seq[St][c]
        # print (c,St,Nuc)		 
         if ATGC.count(Nuc)==1:		 
            Nuc2Cou[Nuc]=Nuc2Cou.get(Nuc,0)+1
   #   print (c,Nuc2Cou)	
      	
      if len(Nuc2Cou)==0: out+=str(c+1)+'\tNA\tNA\tNA\tNA\tNA\tNA\n'
      elif len(Nuc2Cou)==1:	  
        out+=str(c+1)
        for Nuc in Nuc2Cou:
            Cou=Nuc2Cou[Nuc]		
            Freq=1.0*Cou/SeqC	
            out+='\t'+Nuc+'\t'+str(Cou)+'\t'+str(Freq)
        out+='\tNA\tNA\tNA\n'		
      else: 	  
       Cou2Nuc={}
       for Nuc in Nuc2Cou:
          Cou=Nuc2Cou[Nuc]
          Cou2Nuc[Cou]=Cou2Nuc.get(Cou,[])+[Nuc]
       CouLs=list(Cou2Nuc.keys())
       CouLs.sort(reverse=True)
     #  print (CouLs)	  
       In=0
       i=0	  
       while In<2:
         Cou=CouLs[i]	   
         NucLs=Cou2Nuc[Cou]
        # print (c,Nuc2Cou,Nuc2Cou,NucLs,Cou)
       #  open('A','r').readlines()		 
         #NucLs=Cou2Nuc[Cou]
         for Nuc in NucLs:
            Freq=1.0*Cou/SeqC		 
            if In>1: pass		 
            elif In==0:
			
               out+=str(c+1)+'\t'+Nuc+'\t'+str(Cou)+'\t'+str(Freq)+'\t'
            else: out+=Nuc+'\t'+str(Cou)+'\t'+str(Freq)+'\n'
            In+=1
         i+=1
      c+=1		 
   OutF=Fas[:-6]+'_PosInfo.txt'
   GetOut(OutF,out)
def MakeHapFas(Fas,PosLs): #Pos from 1
    StLs,St2Seq=ReadFasSeq(Fas)
    out=''	
    for St in StLs:
        out+=St+'\n'
        Seq=St2Seq[St]
        Hap=''
        NAPos=[]		
        for Pos in PosLs:
            if len(Seq)<(Pos-1):
                 NAPos.append(Pos)			
                 #print (Pos,' not available',Pos)
                 Hap+='?'				 
            else: Hap+=Seq[Pos-1]
        out+=Hap+'\n'
    print ('Pos not available',set(NAPos))		
    GetOut(Fas[:-6]+'_Hap.fasta',out)   
def AddOutGroupSeq(OriFas,OutGFas):
   StLs,St2Seq=ReadFasSeq(OriFas)
   OutLs,Out2Seq=ReadFasSeq(OutGFas)
   out=''
   for i in OutLs:
       out+=i+'\n'+Out2Seq[i]+'\n'
   for i in StLs:
       out+=i+'\n'+St2Seq[i]+'\n'
   GetOut(OriFas[:-6]+'_without.fasta',out)	   
def InvertDic(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=St
 return Hap2ID   
 
def CountDifNum_RmMiss(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif+=1
                c+=1
            return Dif
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif

def prune_tree(OriNwk,ExtraLs):
   TreeLs=open(OriNwk,'r').readlines()
   Cou=1 
   for Tst in TreeLs:
    tree=Phylo.read(StringIO(Tst.strip()), "newick")   
  #  tree = Phylo.read(OriNwk, 'newick')
    for tip in ExtraLs:	
         tree.prune(tip)
    Phylo.write(tree, OriNwk[:-4]+'_'+str(Cou)+'_prune.nwk','newick')
    Cou+=1	
def Count_item(File,Col):
    File=open(File,'r').readlines()[1:]
    Item2C={}
    for i in File:
       i=i.strip().split('\t')
       Item2C[i[Col]]=Item2C.get(i[Col],0)+1
    return Item2C
def TopHap_haplotype_filtering(MinSeqCFre,pos72): #pos72=fasta
    HapOut=pos72[:-6]+'_hapID.fasta'
    StLs=[]
    St2Seq={}
    pos72=open(pos72,'r').readlines()
    for i in pos72:
        if i[0]=='>': 
           ID=i.strip()
           StLs.append(ID)
           St2Seq[ID]=''
        else: St2Seq[ID]+=i.strip()
    print (len(StLs))
    MinSeqC=MinSeqCFre*len(StLs)
    print ('min seqC',MinSeqC)
    Seq2C={}  
    for St in StLs:
        Seq=St2Seq[St]
        Seq2C[Seq]=Seq2C.get(Seq,0)+1
    HapID2Seq={}
    Seq2HapID={}
    Hap=1
    out=''
    outFil=''
    for Seq in Seq2C:
        C=Seq2C[Seq]
        HapID='Hap'+str(Hap)+'_'+str(C)             
        out+='>'+HapID+'\n'+Seq+'\n'
        if C>MinSeqC: outFil+='>'+HapID+'\n'+Seq+'\n' 
        Seq2HapID[Seq]=HapID
        Hap+=1
    OutF=open(HapOut,'w')
    OutF.write(out)
    OutF.close()
    OutF=open(HapOut[:-6]+'_Filter'+str(MinSeqC)+'.fasta','w')
    OutF.write(outFil)
    OutF.close()        
    out='ID\tHapID\n'
    for St in StLs:
        Seq=St2Seq[St]
        HapID=Seq2HapID[Seq]
        out+=St.replace('>','')+'\t'+HapID+'\n'
    OutF=open(HapOut[:-6]+'.txt','w')
    OutF.write(out)
    OutF.close()    
def InvertDic(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=St
 return Hap2ID        
def BootstrapSeq(BooC,FillPer,OriTreeFile,OriFas,Rpath):
    #BooC=100 #100
    #FillPer=0.03#0.05#24#40 for 68k 24 for 29k 
    #Cut=90 #90
    
    Fol='Alignment\\'
    OutFolMain='Bootstrap\\'#'C:\\Users\\kumarlab\\Desktop\\TopHap\\68k_SNV72\\Boo\\'
    
    #OriTree=Fol+'TopHap.nwk'
    #OriFas=Fol+'TopHap.fasta'
    OutGFas=Fol+'OutG.fasta'#'OutGHap.fasta'
    FasTaLs=[Fol+'Note.txt',Fol+'Note-Continent.txt',Fol+'Note-Country.txt','AdditionalData.txt']
    PosDic={} #pos from 1
    Boo2FasLs={}
    go='n'
    BooFolLs=[]    
    for FasLs in FasTaLs:
        print (FasLs)
        if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
        else: Name=''	
        FasLs=open(FasLs,'r').readlines()[1:]
        for i in FasLs:
            Fas0=i.split('\t')[0]
            if Fas0!='NA' and Fas0.find('-not used')==-1:	
                OutFol=OutFolMain+Fas0.replace(' ','_')+'\\'
                if os.path.exists(OutFol)!=True: os.mkdir(OutFol)
                BooFolLs.append(OutFol)            			
                Fas=Fol+Name+Fas0.replace(' ','_')+'_Hap.fasta'	
                StLs,St2Seq=ReadFasSeq(Fas)
                Len=len(StLs)
                print (Len)
                Fill=FillPer*Len
                print ('seq C to keep TopHap',Fill,FillPer)
                Boo=1
                Hap2RepC={}
                while Boo<=BooC:
                  if os.path.exists(OutFol+'Rep'+str(Boo)+'.fasta')!=True:			
                    Samp=list(np.random.choice(StLs,replace=True,size=Len))
                    print ('sampling seq',Boo,len(Samp),len(set(Samp)))
                    out=''
                    SeqID=1	
                    for i in Samp:
                        out+=i+'_'+str(SeqID)+'\n'+St2Seq[i]+'\n'
                        SeqID+=1	
                    GetOut(OutFol+'Rep'+str(Boo)+'.fasta',out)
                    print ('getting top haplotype')	
                  if os.path.exists(OutFol+'Rep'+str(Boo)+'_hapID_Filter'+str(Fill)+'.fasta')!=True:  TopHap_haplotype_filtering(FillPer,OutFol+'Rep'+str(Boo)+'.fasta')	#out=Rep1_hap_Filter40.fasta
                  HapLs,Hap2Seq=ReadFasSeq(OutFol+'Rep'+str(Boo)+'_hapID_Filter'+str(Fill)+'.fasta')
                  Boo2FasLs[Boo]=Boo2FasLs.get(Boo,[])+[OutFol+'Rep'+str(Boo)+'_hapID_Filter'+str(Fill)+'.fasta']				
                  Boo+=1
    print ('make fasta for each boo')
    OriFas_pru=OriFas[:-6]+'_prune.fasta'
    AllRepTrees=''    
    Boo=1
    RepFasLs=[]
    out=''
    while Boo<=BooC:
      FasLs=Boo2FasLs[Boo]
      Hap2C={}
      for Fas in FasLs:
          IDLs,ID2Seq=ReadFasSeq(Fas)
          for ID in IDLs:
              Hap=ID2Seq[ID]	  
              Hap2C[Hap]=Hap2C.get(Hap,0)+1
      ID=1
      for Hap in Hap2C:
          out+='>A'+str(ID)+'_'+str(Hap2C[Hap])+'\n'+Hap+'\n'
          ID+=1
      if os.path.exists(OutFolMain+'Res'+str(Boo)+'.fasta') !=True: GetOut(OutFolMain+'Res'+str(Boo)+'.fasta',out)
      RepFasLs.append(OutFolMain+'Res'+str(Boo)+'.fasta') 
      Boo+=1  
    print ('boo C',len(RepFasLs))
    StLs,St2Seq=ReadFasSeq(OriFas)
    Hap2ID=InvertDic(St2Seq)	
    Hap2Count={}
    for St in StLs:
    
        Hap2Count[St]=0
    OutLs,Out2Seq=ReadFasSeq(OutGFas)
    OutG=''
    for i in OutLs:
        OutG+=i+'\n'+Out2Seq[i]+'\n'
    
    for RepFas in RepFasLs:
        RStLs,RSt2Seq=ReadFasSeq(RepFas)
        Rhap2ID=InvertDic(RSt2Seq)	
    
        for Rhap in Rhap2ID:	
            if Rhap in Hap2ID: Hap2Count[Hap2ID[Rhap]]+=1
    ExtraTip=[]
    out=OutG
    for Hap in Hap2Count:
        if Hap2Count[Hap]!=len(RepFasLs) and OutLs.count(Hap)==0: ExtraTip.append(Hap.replace('>',''))
        elif OutLs.count(Hap)==0: out+=Hap+'\n'+St2Seq[Hap]+'\n'
    GetOut(OriFas_pru,out)
    #open('A','r').readlines()	
    print ('extra tip',ExtraTip)
    #for OriTree in OriTreeLs:	
    prune_tree(OriTreeFile,ExtraTip)	
    	
    for RepFas in RepFasLs:
        #print RepFas
       # open('A','r').readlines()	
        RStLs,RSt2Seq=ReadFasSeq(RepFas)
        Rhap2ID=InvertDic(RSt2Seq)	
        out=OutG	
        ID=RepFas.split('\\')[-1].split('_')[0]
        IDC=1	
     #   print (Rhap2ID)	
        ExtraLs=[]
        for Rhap in Rhap2ID:	
            if Rhap in Hap2ID: 
               NewID=Hap2ID[Rhap].replace('>','')
              # print (NewID)		   
               if ExtraTip.count(NewID)!=0: ExtraLs.append(NewID)		   
            else:
                NewID=ID+'_'+str(IDC)
                ExtraLs.append(NewID)			
                IDC+=1
    			
            out+='>'+NewID+'\n'+Rhap+'\n'
        if os.path.exists(RepFas[:-6]+'_Rename.fasta')!=True: GetOut(RepFas[:-6]+'_Rename.fasta',out)
        if os.path.exists(RepFas[:-6]+'_Rename.nwk')!=True: 
          os.system('megacc -a infer_MP_nucleotide.mao -d ' + RepFas[:-6]+'_Rename.fasta' +' -o ' + RepFas[:-6]+'_Rename.nwk')	
          AncFileLs=glob.glob(RepFas[:-6]+'_Rename_ancestral_states_*.txt')
          for AncF in AncFileLs:
               os.remove(AncF)	
        #  os.remove(RepFas[:-6]+'_summary.txt')			   
         # print (RepFas[:-6]+'_Rename.nwk')	
    	  
        #  open('A','r').readlines()
        #else:
        print (ID,len(ExtraLs),len(RStLs),len(RStLs)-len(ExtraLs))	
        if os.path.exists(RepFas[:-6]+'_Rename.nwk')==True: prune_tree(RepFas[:-6]+'_Rename.nwk',ExtraLs)
        PruTreeLs=glob.glob(RepFas[:-6]+'_Rename_*_prune.nwk')	
        for PruTre in PruTreeLs:
            TrLines=open(PruTre,'r').readlines()	
            for Tr in TrLines:			
                AllRepTrees+=Tr	
          #  os.remove(PruTre)
			
   # print ('Run sumTree.R and ReNameTreeFas.py')	
    GetOut('AllResampTrees.nwk',AllRepTrees)
    OriPruneLs=glob.glob(OriTreeFile[:-4]+'_*_prune.nwk')
    print ('target tres',OriPruneLs)	
    trees = list(Phylo.parse('AllResampTrees.nwk', "newick"))	
    c=1	
    for OriTreFile in OriPruneLs:	
        print ('score tree',OriTreFile,ID)	 	
        BooSum(OriTreFile,str(c),Rpath)	
  
     #   target_tree =Phylo.read(StringIO(open(OriTreFile,'r').readlines()[0].strip()), "newick")  
      #  support_tree = get_support(target_tree, trees)
       # Phylo.write(support_tree, 'TopHap_bootstrap'+str(c)+'.nwk','newick')
        c+=1
        os.remove(OriTreFile)
    print (BooFolLs)		
    for i in BooFolLs:
       Files=glob.glob(i+'Rep*')
       for i in Files:
          os.remove(i)      
    Files=glob.glob('Bootstrap\\Res*')
    for i in Files:
       os.remove(i)
def BooSum(OriTreFile,ID,Rpath):
    cwd = os.getcwd()
    cwd=cwd.replace('\\','/')
    Rout=''
    Rout+='library(ape)\n'
    Rout+='true_tree <- ape::read.tree(\"'+cwd+'/'+OriTreFile+'\")\n'
    Rout+='lf <- list.files(\"'+cwd+'/Bootstrap/\", pattern = \"_prune.nwk\", full.names = TRUE)\n'
    Rout+='x<- ape::rmtree(length(lf), true_tree$Nnode)  ## length(lf) = num of replicate trees\n'
    Rout+='for(j in 1:length(lf)){\n'
    Rout+='  x[j] <- list(ape::read.tree(lf[j]))\n'
    Rout+='}\n'
    Rout+='b <- phangorn::plotBS(true_tree, x, p =10,  \"phylogram\")\n'
    Rout+='true_boot_support <- as.numeric(b$node.label)\n'
    Rout+='ape::write.tree(b, file = \"'+cwd+'/TopHap_bootstrap'+ID+'.nwk\")\n'
    GetOut('RunR.r',Rout)
    os.system(Rpath+' RunR.r')
   # open('A','r').readlines()	
    os.remove('RunR.r')	
def GetOut(OutFile,outIn):
 OutF=open(OutFile,'w')
 OutF.write(outIn)
 OutF.close() 	