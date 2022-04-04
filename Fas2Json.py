import sys

#count from 0
Fa=sys.argv[1]
OutFa=sys.argv[2]
OutSeq=open(OutFa,'r').readlines()[0].strip()
Len=len(OutSeq)

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
def GetOut(OutFile,outIn):
 OutF=open(OutFile,'w')
 OutF.write(outIn)
 OutF.close() 	 
CellLs,Cell2Seq=ReadFasSeq(Fa) 


In=[]
for Cell in CellLs:
      Seq=Cell2Seq[Cell]
      if len(Seq)!=len(OutSeq):		
        print ('skiped because length is different: ',Cell,len(Seq))
      else:		
        VarLs=[]
        C=0
        while C<Len:
            if Seq[C]!=OutSeq[C]: VarLs+=[str(C),'\"'+Seq[C]+'\"']
            C+=1
        In.append('{'+'\"V\":['+','.join(VarLs)+'],\"I\":[\"20220101\",\"'+Cell.replace('>','')+'\",\"NA\",\"'+Cell.replace('>','')+'\",\"'+Cell.replace('>','')+'\"]}')
out='{\"sequences\":['+','.join(In)+']}'
GetOut(Fa[:-3]+'.json',out)	
		
        		
    	