import Functions3
import os
import glob
import sys
import time

StartTime='Start time: '+time.strftime("%Y-%m-%d %H:%M")+'\n'
THCutFre=0.05 #hf
BooC=100 #bootstrap replicates
BooFillPer=0.05 #hf for bootstrap replicates
Rpath='Rscript'

if os.path.exists('Bootstrap')==True:
     print ('Please delete Bootstrap directory and try again')
      

else:
 os.mkdir('Bootstrap')
 if sys.argv.count('-Hap')==1:
    Fol=sys.argv[sys.argv.index('-Hap')+1]+'\\'
    FasLs0=glob.glob(Fol+'*_Hap.fasta')
    out='Input\n'
    for Fas in FasLs0:
       out+=Fas.split('\\')[-1]+'\n'
    Functions3.GetOut('InputFasList.txt',out)
    FasTaLs=['InputFasList.txt']
    HapGive='y'	

 TopHap2Data={}
 for FasLs in FasTaLs:
    print (FasLs)
    Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
            print ('making TopHap sequences for the slice',i)	
            Fas0=i.strip().split('\t')[0]
            Fas=Fol+Fas0[:-10]+'.fasta'			
            Seq2Cou={}
            StLs,St2Seq=Functions3.ReadFasSeq(Fas[:-6]+'_Hap.fasta')
            Fill=1.0*len(StLs)*THCutFre			
            for i in StLs:
                       Seq2Cou[St2Seq[i]]=Seq2Cou.get(St2Seq[i],0)+1	
				
            for Seq in Seq2Cou:
					 
                  if Seq2Cou[Seq]>=Fill:			
                       TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]			
					  
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
 Ols,O2Seq=Functions3.ReadFasSeq(Fol+'OutG.fasta')
 outG=''
 for O in Ols:
    Seq=O2Seq[O]
    outG +=O+'\n'+Seq+'\n'
 outF=outG+outF	
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
