import Functions3
import os
import sys
import glob
import shutil


OriFas=sys.argv[1]
OriTree=sys.argv[2]
QueryFas=sys.argv[3]


RaxMLdir=sys.argv[4]
CommandRaxML=RaxMLdir+' -f v -s REFQUEFAS -t REFREE -m GTRGAMMA -n OUT -T 2'
 

OutLbLs=['Outgroup1','Outgroup2','Outgroup3']
OutGSeqID=OutLbLs[0] 
wdir= os.path.join(os.path.expanduser('~'), os.getcwd())
TMPdir=wdir+ os.sep+ 'RaxMLoutputsAll'
if os.path.exists(TMPdir)!=True:
    os.mkdir(TMPdir)
TMPdir+= os.sep	
OUTdir=wdir+os.sep+'RaxMLoutputs'	
if os.path.exists(OUTdir)!=True:
    os.mkdir(OUTdir)
OUTdir+= os.sep	
OriTree=open(OriTree,'r').readlines()[0]
Qls,Q2seq=Functions3.ReadFasSeq(QueryFas)
Ols,O2seq=Functions3.ReadFasSeq(OriFas)
out='Query\tMostSimRef\tDifCount\tDiffPos(from 1)\tExtMutLs\tRecMutLs\tMultiMutLs\tCoinfection\tAttached\n'
out1=out#'Query\tMostSimRef\tDifCount\tMethod\tFirstSplit\tCategory\n'
	   

TopHapSeqLs=[]			   
for Q in Qls:
  Go='y'

  if Go=='y':
    Qseq=Q2seq[Q]
    Osimid,DifC=Functions3.GetMostSim(Qseq,O2seq)  
    print (Q,Osimid,DifC)
    if DifC==0:
       TopHapSeqLs.append(Q)	

    else:

           OsimidLs=Osimid.split(';')
           Dif=[]		   
           for Oid in OsimidLs:		   
              Oseq=O2seq['>'+Oid]
              DifPosLs=Functions3.CountDifPos_from1(Qseq,Oseq)
              Dif.append(';'.join(map(str,DifPosLs)))			  
           BaseInfo=Q+'\t'+Osimid+'\t'+str(DifC)+'\t'+';;'.join(Dif)+'\t'
           InTree=	TMPdir+Q.replace('>','')+'.nwk'	   
           Functions3.GetOut(InTree,OriTree)	
           NewFas=TMPdir+Q.replace('>','')+'.fas'		   
           Functions3.makeRaxMLin(O2seq,{Q:Qseq},NewFas,'')			   
           CL=CommandRaxML.replace('-n OUT','-n '+Q.replace('>','')).replace('REFQUEFAS',NewFas).replace('REFREE',InTree)
           print (CL)

           os.system(CL)	
           OutLs=glob.glob('RAxML_*.'+Q.replace('>',''))
           for Out in OutLs:
               shutil.copyfile(Out, TMPdir+Out)		   
               os.remove(Out)	
           os.remove('RAxML_portableTree.'+Q.replace('>','')+'.jplace')			   
           NewTree_str=open(TMPdir+'RAxML_labelledTree.'+Q.replace('>','')).readlines()[0].replace('QUERY___','')	
           Functions3.cleanRaxMLtree(NewTree_str,OUTdir+'RAxML_labelledTree.'+Q.replace('>','')+'.nwk',OutGSeqID)
			   

os.remove('Unroot_rooted.nwk')  
print ('the following quesries were not attached because they were identical to TopHap sequencies:')
print ('\n'.join(TopHapSeqLs))	
		  
	 
