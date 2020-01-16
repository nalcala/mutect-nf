import pysam
import re
import sys

regex = re.compile('[0-9]+N[0-9]+D[0-9]+N') 
sj = ''
infile = pysam.AlignmentFile(sys.argv[1], "rb")
outfile = pysam.AlignmentFile(sys.argv[2], "wb", template=infile)
for s in infile:
    s2 = s
    whMC = [i for i in range(len(s2.tags)) if s2.tags[i][0]=="MC" ]
    if len(whMC)>0:
        ctmp = s2.get_tag('MC')
        rctmp = regex.search(ctmp)
        while rctmp!=None:
            ctmpshort = ctmp[rctmp.start():rctmp.end()]
            indsshort  = [i for i in range(len(ctmpshort)) if ctmpshort[i].isalpha()]
            ndn = int( ctmpshort[0:indsshort[0]] ) + int( ctmpshort[(indsshort[0]+1):indsshort[1]]) + int( ctmpshort[(indsshort[1]+1):indsshort[2]] )
            ctmp = ctmp[0:rctmp.start()]+str(ndn)+"N"+ctmp[rctmp.end():len(ctmp)]
            #print(ctmp)
            #inds  = [(i,ctmp[i]) for i in range(len(ctmp)) if ctmp[i].isalpha()]
            #inds2 = [(ctmp[0:inds[0][0]],inds[0][1]) ] + [(ctmp[(inds[i][0]+1):inds[i+1][0]],inds[i+1][1]) for i in range(len(inds)-1)]
            #ndn   = [(i,int(inds2[i][0]) + int(inds2[i+1][0]) + int(inds2[i+2][0])) for i in range(len(inds2)-2) if (inds2[i][1]=="N" and inds2[i+1][1]=="D") and inds2[i+2][1]=="N" ]
            #inds3 = inds2[0:(ndn[0][0])] + [(str(ndn[0][1]),"N")] + inds2[(ndn[0][0])+3 : len(inds2) ]
            #ctmp2 = sj.join([sj.join(x) for x in inds3 ])
            rctmp = regex.search(ctmp)
        s2.set_tag( 'MC', ctmp )
    outfile.write(s2)
infile.close()
outfile.close()