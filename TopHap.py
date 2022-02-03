import Functions3
import os
import glob
import sys
import time

StartTime='Start time: '+time.strftime("%Y-%m-%d %H:%M")+'\n'
THCutFre=float(sys.argv[1])#0.05 #hf
BooC=int(sys.argv[2])#100 #bootstrap replicates
BooFillPer=THCutFre#0.05 #hf for bootstrap replicates
Rpath='Rscript'
OutGfas='Outgroup.fasta'
PosStart=0
path = os.path.join(os.path.expanduser('~'), os.getcwd())

if os.path.exists('Bootstrap')==True:
     print ('Please delete Bootstrap directory and try again')
      

else:
 os.mkdir('Bootstrap')
 if sys.argv.count('-Hap')==1:
    Fol=sys.argv[sys.argv.index('-Hap')+1]+os.sep
    FasLs0=glob.glob(Fol+'*.fasta')
    out='Input\n'
    for Fas in FasLs0:
       InLs,InSeq=Functions3.ReadFasSeq(Fas)	
       if Fas!=Fol+'OutG.fasta' and len(InLs)>=500:	
          out+=Fas.split(os.sep)[-1]+'\n'
    Functions3.GetOut('InputFasList.txt',out)
    FasTaLs=['InputFasList.txt']
    	
 else: print ('please provide the directory of haplotype alignments')	

				  
 TopHap2Data={}
 for FasLs in FasTaLs:
    print (FasLs)
    Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
            print ('making TopHap sequences for the slice',i)	
            Fas0=i.strip().split('\t')[0] ####Fix
            Fas=Fol+Fas0	 ####Fix		
            Seq2Cou={}
            StLs,St2Seq=Functions3.ReadFasSeq(Fas)
            Fill=1.0*len(StLs)*THCutFre			
            for i in StLs:
                       Seq2Cou[St2Seq[i]]=Seq2Cou.get(St2Seq[i],0)+1	
				
            for Seq in Seq2Cou:
					 
                  if Seq2Cou[Seq]>=Fill:			
                       TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]	####Fix	

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
 PosLs0=open(Fol+'Haplotypes.txt','r').readlines()
 PosLs=[]
 for i in PosLs0: 
       if PosStart==0: PosLs.append(int(i)+1)
       else: PosLs.append(int(i))	
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
 Ols,O2Seq=Functions3.ReadFasSeq(Fol+'OutG.fasta')
	
 Functions3.GetOut('TopHap.fasta',outF)

        	
 print ('infer MP')
 os.system('megacc -a infer_MP_nucleotide.mao -d ' + 'TopHap.fasta' +' -o ' + 'TopHap_allMP.nwk')
 AncFileLs=glob.glob('TopHap_allMP_ancestral_states_*.txt')
 for AncF in AncFileLs:
               os.remove(AncF)	  
 if os.path.exists('TopHap_allMP_summary.txt')==True: os.remove('TopHap_allMP_summary.txt')			   
 print ('do bootstrap')

 Functions3.BootstrapSeq_afterBoo(BooC,BooFillPer,'TopHap_allMP.nwk','TopHap.fasta',Rpath,Fol,FasTaLs)


 os.remove('AllResampTrees.nwk')	
 os.remove('Rplots.pdf')	

 Ls=glob.glob('TopHap_allMP_*_prune.nwk')
 for i in Ls:
    os.remove(i)
 Ls=glob.glob('TopHap_allMP_changes_list_*.txt')
 for i in Ls:
    os.remove(i)	

 os.remove('rooted.nwk') 
 os.remove('TopHap.fasta') 
 os.remove('TopHap_allMP.nwk') 
EndTime='End time: '+time.strftime("%Y-%m-%d %H:%M")+'\n'
Functions3.GetOut('Time.txt',StartTime+EndTime)	
