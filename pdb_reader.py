import numpy as np
from collections import Counter
from collections import OrderedDict

def aminofrequency(line):
    frequency = dict(Counter(line.split()))
    frequency = OrderedDict(sorted(frequency.items(), key=lambda t: t[0]))
    return frequency


def new_dihedral(p):
    p0 = np.array(p[0])
    p1 = np.array(p[1])
    p2 = np.array(p[2])
    p3 = np.array(p[3])
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    b1 /= np.linalg.norm(b1)
    # v = projection of b0 onto plane perpendicular to b1
    # w = projection of b2 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the dihedral angle
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
    

#LISTS USED
l3  = []
l2  = []
helix = []
sheet = []
lig = []
chain = ['#']
amino = ''
TITLE = ''
title = []

#READING PDB FILE
for line in open('2wsc.pdb'):
    l1 = list(line.split())
    l3.append(l1)

target = open("2wsc_output.txt",'w')

sum = 0

for i in l3:
    if (i[0] == "SEQRES"):
        if (chain[-1] != i[2]):
            chain.append(i[2])
        for j in range (0,4):
            i.pop(0)
        sum += len(i)
        for k in i:
            amino = amino + " " + k
    elif (i[0] == "FORMUL"):   
        lig.append(i[2])
    elif (i[0] == "HELIX"):
        helix.append(str("HELIX" + "\t" + i[5] + "," + i[8] + "    " + "\t" + i[-4] + "\t" + i[-1]))
    elif (i[0] == "SHEET"):
        sheet.append(str("SHEET" + "\t" + i[5] + "," + i[8] + "    " + "\t\t" + i[10]))
    elif (i[0] == "TITLE"):
            l2 = []
            l2.append(i)
            l2[0].pop(0)
            if (TITLE == ''):
                TITLE += ' '.join(l2[0])
            else:
                l2[0].pop(0)
                TITLE += ' ' + ' '.join(l2[0])
    else:  
        continue    
TITLE += '\n'
lig.sort()
ligandcons = ' '.join(lig)

aminodic = aminofrequency(amino)
if 'UNK' in aminodic:
    a = int(aminodic['UNK'])
else:
    a = 0;
chain.pop(0)
chain.sort()
chainconst = ','.join(chain)

target.write (TITLE)

str1 = 'LENGTH'+'\t'+str(sum)+'\n'
target.write (str1)

str1 = 'CHAINS'+'\t'+str(len(chain))+'\t'+str(chainconst)+'\n'
target.write (str1)

for key in aminodic:
        str1 = str(key)+'\t'+str(round((float(aminodic[key])/sum),13))+'\n'
        target.write (str1)

str1 = 'UNKNOWN '+str(a)+'\n'
target.write (str1)

for i in helix: 
    target.write (i)
    target.write ('\n')

for i in sheet:
    target.write (i)
    target.write ('\n')

target.write ('LIGANDS ')
target.write (ligandcons)
target.write ('\n')
C = []
N = []
CA = []
ch = []

for i in l3:
    if (i[0] == "ATOM"):
        #if (ch[-1] == '#'):
         #   ch.append(i[4])
        #elif (ch[-1] != i[4]):
          #  ch.append(i[4])
           # angleprint(C,N,CA)
        cd = []
        cd.append(float(i[6]))

        cd.append(float(i[7]))
        cd.append(float(i[8]))

        if (i[2] == "C"):
            C.append(cd)
        elif (i[2] == "N"):
            ch.append(i[4])
            N.append(cd)
        elif (i[2] == "CA"):
            CA.append(cd)

strl = "CHAIN-"+ch[1]
phi = [strl]
phi.append('NA')
psi = ['']
omega = ['']

for i in range(0,len(C)-1):
    ph = [];
    ph.append(C[i])
    ph.append(N[i+1])
    ph.append(CA[i+1])
    ph.append(C[i+1])
    phi.append(new_dihedral(ph))
    ps = [];
    ps.append(N[i])
    ps.append(CA[i])
    ps.append(C[i])
    ps.append(N[i+1])
    psi.append(new_dihedral(ps))
    o = [];
    o.append(CA[i])
    o.append(C[i])
    o.append(N[i+1])
    o.append(CA[i+1])
    omega.append(new_dihedral(o))
    #if (i > 0 and ch[i+1] != ch [i]):
    if ( ch[i] != ch [i+1]):
        strl = "CHAIN-" + ch[i+1]
        phi[-1] = (strl)
        phi.append('NA')
        omega.pop(-1)
        psi.pop(-1)
        omega.append('NA')
        psi.append('NA')
        omega.append("")
        psi.append("")

psi.append('NA')
omega.append('NA')
        

for i in (zip(phi,psi,omega)):
    str1 =str(i[0])+'        \t'+str(i[1]) +'      \t'+str(i[2])+'\n'
    target.write (str1)
