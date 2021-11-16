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

def prune_tree(OriNwk,ExtraLs,OutSeq):
   if os.path.exists('rooted.nwk')==True: os.remove('rooted.nwk')	
   trees = list(Phylo.parse(OriNwk, 'newick'))
   for tree in trees:  
       tree = tree.root_with_outgroup({'name': OutSeq})
   Phylo.write(trees, 'rooted.nwk', "newick")     
   TreeLs=open('rooted.nwk','r').readlines()
   Cou=1 
   for Tst in TreeLs:
    tree=Phylo.read(StringIO(Tst.strip()), "newick") 
    for tip in ExtraLs:	
         tree.prune(tip)
    Phylo.write(tree, OriNwk[:-4]+'_'+str(Cou)+'_prune.nwk','newick')
    Cou+=1

def root_tree(OriNwk,Root):
    Out=OriNwk[:-4]+'_rooted.nwk'
    trees = list(Phylo.parse(OriNwk, 'newick'))
    for tree in trees:
       tree = tree.root_with_outgroup({'name': Root})
    Phylo.write(trees, Out, "newick")


def InvertDic(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=St
 return Hap2ID  
def InvertDic1(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
 return Hap2ID   

def BootstrapSeq_afterBoo(BooC,FillPer,OriTreeFile,OriFas,Rpath,Fol,FasTaLs):
    OutFolMain='Bootstrap\\'#'C:\\Users\\kumarlab\\Desktop\\TopHap\\68k_SNV72\\Boo\\'
    OutGFas=Fol+'OutG.fasta'#'OutGHap.fasta'
    PosDic={} #pos from 1
    Boo2FasLs={}
    go='n'
    BooFolLs=[]  
    print (FasTaLs)
    for FasLs in FasTaLs:
        if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
        else: Name=''	
        FasLs=open(FasLs,'r').readlines()[1:]
        for i in FasLs:
                print ('making bootstrap TopHap seqs',i)		   
                Fas0=i.strip().split('\t')[0]
                OutFol=OutFolMain+Fas0.replace(' ','_')+'\\'
                if os.path.exists(OutFol)!=True: os.mkdir(OutFol)          			
                Fas=Fol+Name+Fas0.replace(' ','_')#+'_Hap.fasta'				
                StLs,St2Seq=ReadFasSeq(Fas)
                Len=len(StLs)
                Fill=FillPer*Len
                Boo=1
                 
                while Boo<=BooC:			  			  
                    Samp=list(np.random.choice(StLs,replace=True,size=Len))
                    Seq2Cou={}					
                    for i in Samp:
                       Seq2Cou[St2Seq[i]]=Seq2Cou.get(St2Seq[i],0)+1	
                    out=''
                    IDTop=1					
                    for Seq in Seq2Cou:
					 
                         if Seq2Cou[Seq]>=Fill:
                              out+='>'+str(IDTop)+'\n'+Seq+'\n'							  
                              IDTop+=1							   
                    GetOut(OutFol+'Rep'+str(Boo)+'TopHap.fasta',out)							  			
                    Boo+=1
    print ('make fasta for each boo')
    OriFas_pru=OriFas[:-6]+'_prune.fasta'
    AllRepTrees=''    
    Boo=1
    RepFasLs=[]

    while Boo<=BooC:
      print (Boo,'make fasta for each boo')	
      TSeq2C={}	
      for FasLs in FasTaLs:	
        FasLs=open(FasLs,'r').readlines()[1:]	  
        for i in FasLs:
              #  print (i)		   
                Fas0=i.strip().split('\t')[0]
                TopHapSeq=OutFolMain+Fas0.replace(' ','_')+'\\'	+'Rep'+str(Boo)+'TopHap.fasta'
                Tls,T2Seq=ReadFasSeq(TopHapSeq)
                for T in T2Seq:
                     Seq=T2Seq[T]
                     TSeq2C[Seq]=TSeq2C.get(Seq,0)+1	
      
      ID=1
      out=''	  
      for Hap in TSeq2C:
          out+='>A'+str(ID)+'_'+str(TSeq2C[Hap])+'\n'+Hap+'\n'
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
        if Hap2Count[Hap]!=len(RepFasLs) and OutLs.count(Hap)==0: ExtraTip.append(Hap.replace('>','').replace('__','_'))
        elif OutLs.count(Hap)==0: out+=Hap+'\n'+St2Seq[Hap]+'\n'
    GetOut(OriFas_pru,out)
    print ('extra tip',ExtraTip,len(ExtraTip))
    prune_tree(OriTreeFile,ExtraTip,OutLs[0].replace('>',''))		
    for RepFas in RepFasLs:
        RStLs,RSt2Seq=ReadFasSeq(RepFas)
        Rhap2ID=InvertDic(RSt2Seq)	
        out=OutG	
        ID=RepFas.split('\\')[-1].split('_')[0]
        IDC=1	
        ExtraLs=[]
        for Rhap in Rhap2ID:	
            if Rhap in Hap2ID: 
               NewID=Hap2ID[Rhap].replace('>','').replace('__','_')   
               if ExtraTip.count(NewID)!=0: ExtraLs.append(NewID)		   
            else:
                NewID=ID+'_'+str(IDC)
                ExtraLs.append(NewID)			
                IDC+=1
    			
            out+='>'+NewID+'\n'+Rhap+'\n'
        if os.path.exists(RepFas[:-6]+'_Rename.fasta')!=True: 
              GetOut(RepFas[:-6]+'_Rename.fasta',out)
		  
        if os.path.exists(RepFas[:-6]+'_Rename.nwk')!=True: 
          os.system('megacc -a infer_MP_nucleotide.mao -d ' + RepFas[:-6]+'_Rename.fasta' +' -o ' + RepFas[:-6]+'_Rename.nwk')	
          AncFileLs=glob.glob(RepFas[:-6]+'_Rename_ancestral_states_*.txt')
          for AncF in AncFileLs:
               os.remove(AncF)	

        print (ID,len(ExtraLs),len(RStLs),len(RStLs)-len(ExtraLs))	
        if os.path.exists(RepFas[:-6]+'_Rename.nwk')==True: prune_tree(RepFas[:-6]+'_Rename.nwk',ExtraLs,OutLs[0].replace('>',''))
        PruTreeLs=glob.glob(RepFas[:-6]+'_Rename_*_prune.nwk')	
        for PruTre in PruTreeLs:
            TrLines=open(PruTre,'r').readlines()	
            for Tr in TrLines:			
                AllRepTrees+=Tr	

    GetOut('AllResampTrees.nwk',AllRepTrees)
    OriPruneLs=glob.glob(OriTreeFile[:-4]+'_*_prune.nwk')
    print ('target tres',OriPruneLs)	
    trees = list(Phylo.parse('AllResampTrees.nwk', "newick"))	
    c=1	
    for OriTreFile in OriPruneLs:	
        print ('score tree',OriTreFile,ID)	 	
        BooSum(OriTreFile,str(c),Rpath)	

        c+=1
        os.remove(OriTreFile)	   
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

    os.remove('RunR.r')	
def SumSub(AnnoSubFileDic): #TopHap_anno_Sub*.txt
    VarPosLs=[]	
def GetOut(OutFile,outIn):
 OutF=open(OutFile,'w')
 OutF.write(outIn)
 OutF.close() 	