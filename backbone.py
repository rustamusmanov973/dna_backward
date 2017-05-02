#This program is created on the basis of "backward.py" script freely available at "https://github.com/rustamusmanov973/dna_backward/blob/master/backward.py".

import sys, random, math, re, os, itertools
import Mapping

NucleicAcids  = "C G A T U DC DG DA DT DCYT DGUA DADE DTHY CYT GUA ADE THY URA".split()
nucleic_stuff = list(NucleicAcids)

def inter(a1,a2,m,n):

    return([float(a1[0]) +(m*(float(a2[0])-float(a1[0]))/(m+n)), float(a1[1])+(m*(float(a2[1])-float(a1[1]))/(m+n)), float(a1[2])+(m*(float(a2[2])-float(a1[2]))/(m+n)) ])




def kick(x,u):
    return x+(random.random()-0.5)*u

def norm2(a):
    return sum([i*i for i in a])

def norm(a):
    return math.sqrt(norm2(a))

def normalize(a):
    f = norm(a)
    if f < 1e-8:
        return (0,0,0)
    else:
        return [i/f for i in a]

def iprod(a,b):
    return sum([i*j for i,j in zip(a,b)])

def mvmul(A,b):
    return [iprod(a,b) for a in A]

def det(A):
    (a,d,g),(b,e,h),(c,f,k) = A
    return a*(e*k-f*h)-b*(k*d-f*g)+c*(d*h-e*g)

def m_inv(A):
    u,v,w = A
    d = 1.0/det(A)
    I = zip(*(crossprod(v,w),crossprod(w,u),crossprod(u,v)))
    return [[d*i for i in j] for j in I]

def vr(a):
    return [i-round(i) for i in a]

def dist(a,b):
    return math.sqrt(norm2([i-j for i,j in zip(a,b)]))
def vsub(a,b):
    return [i-j for i,j in zip(a,b)]

def vadd(a,b):
    return [i+j for i,j in zip(a,b)]

def svmul(s,a):
    return [s*i for i in a]

def crossprod(a,b):
    return a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]

def unbreak(R,box,inv=None):
    if not inv:
        inv = m_inv(box)
    # Subtract coordinates of first atom
    # Convert to box coordinates by multiplying with inverse box
    # Truncate vector to remove box shifts
    # Add back coordinates of first atom
    xyz = [vadd(R[0][4:7],mvmul(box,vr(mvmul(inv,vsub(i[4:7],R[0][4:7]))))) for i in R]
    return [i[:4]+tuple(j) for i,j in zip(R,xyz)]

d2r = 3.14159265358979323846264338327950288/180
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"

def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")

def pdbAtom(a,strict=False):
    # With strict format, the residue field is three characters long
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>       atom name,   res name,     res id, chain,       x,            y,             z
    if strict:
        return (str(a[12:16]),str(a[17:20]),int(a[22:26]),a[21],float(a[30:38])/10,float(a[38:46])/10,float(a[46:54])/10)
    else:
        return (str(a[12:16]),str(a[17:21]),int(a[22:26]),a[21],float(a[30:38]),float(a[38:46]),float(a[46:54]))


def pdbBoxRead(a):
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac) , math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [[0.1*fa, 0, 0], [0.1*fb*cg, 0.1*fb*sg, 0], [wx, wy, wz]]

def cos_angle(a,b):
    p = sum([i*j for i,j in zip(a,b)])
    q = math.sqrt(sum([i*i for i in a])*sum([j*j for j in b]))
    return min(max(-1,p/q),1)


def pdbBoxString(b):
    # Box vectors
    u, v, w  = (b[0],b[3],b[4]), (b[5],b[1],b[6]), (b[7],b[8],b[2])

    # Box vector lengths
    nu,nv,nw = [math.sqrt(norm2(i)) for i in (u,v,w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(cos_angle(v,w))/d2r
    beta  = nu*nw == 0 and 90 or math.acos(cos_angle(u,w))/d2r
    gamma = nu*nv == 0 and 90 or math.acos(cos_angle(u,v))/d2r

    return pdbBoxLine % (10*norm(u),10*norm(v),10*norm(w),alpha,beta,gamma)


def pdbOut(atom,i=1):
    insc    = atom[2]>>20
    resi    = atom[2]-(insc<<20)
    x,y,z   = atom[4:7]
    pdbline = "ATOM  %5d %4s %4s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
    return pdbline%((i,atom[0][:4],atom[1][:4],atom[3],atom[2],10*x,10*y,10*z,1,0))


def groAtom(a):
    # 012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return (str(a[10:15]), str(a[5:10]),   int(a[:5]), " ", float(a[20:28]),float(a[28:36]),float(a[36:44]))


def get_pcc_xyz(r):
    pcc = {}
    for i in r:
        if i[0].strip() == "PX":
            pcc["PX"] = i[4:7]
        if i[0].strip() == "C5X":
            pcc["C5X"] = i[4:7]
        if i[0].strip() == "C1A":
            pcc["C1A"] = i[4:7]
        if i[0].strip() == "C1T":
            pcc["C1T"] = i[4:7]
        if i[0].strip() == "C1G":
            pcc["C1G"] = i[4:7]
        if i[0].strip() == "C1C":
            pcc["C1C"] = i[4:7]


    return pcc

def is_terminal(a, b, box, invbox):
    return not a or not b or dist(a,b,box,invbox) > 0.7


class Structure:
    def __init__(self, other, strict=False):

        lines = open(other).readlines()

        # Try extracting PDB atom/hetatm definitions and set the box
        self.box = None
        rest = []
        self.atoms = [pdbAtom(i, strict) for i in lines if isPDBAtom(i) or rest.append(i)]
        if not self.atoms:
            # This should be a GRO file - get the atom count
            n = int(lines[1]) + 2
            self.atoms = [groAtom(i) for i in lines[2:n]]
            b = [float(i) for i in lines[n].split()] + 6 * [0]  # Padding for rectangular boxes

        else:
            # Make sure there is a box definition
            b = [i for i in rest if i.startswith("CRYST1")]
        print("self.atom s")


        print(self.atoms)
        # Build a residue list
        self.residues = [[self.atoms[0]]]
        for i in self.atoms[1:]:
            if i[1:4] != self.residues[-1][-1][1:4]:
                self.residues.append([])
            self.residues[-1].append(i)

        # Extract the sequence
        self.sequence = [i[0][1].strip() for i in self.residues]
        self.nterm = [False]*len(self.sequence)
        self.cterm = [False]*len(self.sequence)
        print("self.residues")
        print(self.residues)
        print("self.sequence")
        print(len(self.sequence))
        print(self.sequence)
        # PBC handling
        # To 'unbreak' residues, subtract the coordinates of the first atom
        # convert to box coordinates and truncate. Convert back to Cartesian
        # coordinates and add to the coordinates of the firsxt atom.
        A, B = None, None

        # Check for protein chains and breaks
        # List the coordinates for amino acid backbone


        # Set chain backbones based on termini. Begin with a 'chain' unless the first residue is protein.
        backbone = [[]]
        for k in self.residues:
            backbone[-1].append(get_pcc_xyz(k))
        # For each protein chain, positions are estimated for backbone atoms and set as a dictionary:
        # {"N": (x,y,z), "CA": (x,y,z), ...}
        self.backbone = []
        print("backbone")
        print(backbone)
        for chain in backbone:

            # Set a dictionary for each residue. The dictionary will contain entries
            # N, H, CA, HA, C, O
            bb = [dict() for i in chain]


            for i in range(len(bb)):
                if i+1 < len(bb):
                    if self.sequence[i] == 'DA':
                        if 'PX' in chain[i] and 'PX' in chain[i+1]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1A'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1A']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1A'], 2, 1))
                        elif 'PX' in chain[i+1]:
                            bb[i]["O5'"] = chain[i]['C5X']
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1A'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1A']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1A'], 2, 1))
                        elif 'PX' in chain[i]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = chain[i]['C5X']
                            bb[i]["C3'"] = chain[i]['C5X']
                            bb[i]["O3'"] = chain[i]['C5X']
                            bb[i]["C2'"] = chain[i]['C1A']
                            bb[i]["C1'"] = chain[i]['C1A']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1A'], 2, 1))

                    elif self.sequence[i] == 'DT':
                        if 'PX' in chain[i] and 'PX' in chain[i+1]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1T'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1T']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1T'], 2, 1))
                        elif 'PX' in chain[i+1]:
                            bb[i]["O5'"] = chain[i]['C5X']
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1T'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1T']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1T'], 2, 1))
                        elif 'PX' in chain[i]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = chain[i]['C5X']
                            bb[i]["C3'"] = chain[i]['C5X']
                            bb[i]["O3'"] = chain[i]['C5X']
                            bb[i]["C2'"] = chain[i]['C1T']
                            bb[i]["C1'"] = chain[i]['C1T']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1T'], 2, 1))


                    elif self.sequence[i] == 'DG':
                        if 'PX' in chain[i] and 'PX' in chain[i+1]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1G'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1G']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1G'], 2, 1))
                        elif 'PX' in chain[i+1]:
                            bb[i]["O5'"] = chain[i]['C5X']
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1G'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1G']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1G'], 2, 1))
                        elif 'PX' in chain[i]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = chain[i]['C5X']
                            bb[i]["C3'"] = chain[i]['C5X']
                            bb[i]["O3'"] = chain[i]['C5X']
                            bb[i]["C2'"] = chain[i]['C1G']
                            bb[i]["C1'"] = chain[i]['C1G']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1G'], 2, 1))

                    elif self.sequence[i] == 'DC':
                        if 'PX' in chain[i] and 'PX' in chain[i+1]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1C'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1C']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1C'], 2, 1))
                        elif 'PX' in chain[i+1]:
                            bb[i]["O5'"] = chain[i]['C5X']
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 3))
                            bb[i]["C3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 1, 1))
                            bb[i]["O3'"] = tuple(inter(chain[i]['C5X'], chain[i+1]['PX'], 3, 1))
                            bb[i]["C2'"] = tuple(inter(chain[i]['C1C'], chain[i+1]['PX'], 1, 2))
                            bb[i]["C1'"] = chain[i]['C1C']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1C'], 2, 1))
                        elif 'PX' in chain[i]:
                            bb[i]["P"] = chain[i]['PX']
                            bb[i]["O1P"] = chain[i]['PX']
                            bb[i]["O2P"] = chain[i]['PX']
                            bb[i]["O5'"] = tuple(inter(chain[i]['C5X'],chain[i]['PX'],1,1))
                            bb[i]["C5'"] = chain[i]['C5X']
                            bb[i]["C4'"] = chain[i]['C5X']
                            bb[i]["C3'"] = chain[i]['C5X']
                            bb[i]["O3'"] = chain[i]['C5X']
                            bb[i]["C2'"] = chain[i]['C1C']
                            bb[i]["C1'"] = chain[i]['C1C']
                            bb[i]["O4'"] = tuple(inter(chain[i]['C5X'], chain[i]['C1C'], 2, 1))
            print("self.backbone")

            self.backbone.extend(bb)
        print(self.backbone)



# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        if self.func == bool:
            return self.value != None
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

# Description
desc = ""

# Option list
options = [
#   option           type         number     default   description
    ("-f",        Option(str,           1,         None, "Input  GRO/PDB structure")),
    ("-o",        Option(str,           1,         None, "Output GRO/PDB structure")),
    ("-raw",      Option(str,           1,         None, "Projected structure before geometric modifications")),
#    ("-c",        Option(str,           1,         None, "Output GRO/PDB structure of expanded CG beads for position restraints")),
    ("-n",        Option(str,           1,         None, "Output NDX index file with default groups")),
    ("-p",        Option(str,           1,         None, "Atomistic target topology")),
    ("-po",       Option(str,           1,         None, "Output target topology with matching molecules list")),
    ("-pp",       Option(str,           1,         None, "Processed target topology, with resolved #includes")),
    ("-atomlist", Option(str,           1,         None, "Atomlist according to target topology")),
    ("-fc",       Option(float,         1,          200, "Position restraint force constant")),
    ("-to",       Option(str,           1,         None, "Output force field")),
    ("-from",     Option(str,           1,         None, "Input force field")),
    ("-strict",   Option(bool,          0,         None, "Use strict format for PDB files")),
    ("-nt",       Option(bool,          0,         None, "Use neutral termini for proteins")),
    ("-sol",      Option(bool,          0,         None, "Write water")),
    ("-kick",     Option(float,         1,            0, "Random kick added to output atom positions")),
    ]


# Parsing arguments
args = sys.argv[1:]
if '-h' in args or '--help' in args:
    print "\n",__file__
    print desc or "\nSomeone ought to write a description for this script...\n"
    for thing in options:
        print type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing
    print
    sys.exit()


# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])



# Process the command line - list the options that were given
opts = []
while args:
    opts.append(args.pop(0))
    options[opts[-1]].setvalue([args.pop(0) for i in range(options[opts[-1]].num)])

mapping = Mapping.get()

out = open("inter.pdb", "w")

predet = open(options["-f"].value, "r")
for line in predet:
    p1 = re.sub("DTX","DT ", line)
    p2 = re.sub("DCX","DC ", p1)
    p3 = re.sub("DGX","DG ", p2)
    p4 = re.sub("DAX","DA ", p3)
    out.write(p4)
out.close()
predet.close()
struc = Structure("inter.pdb")



counter  =  0
out      = []
cg       = []
raw      = []
sol      = []
ions     = []

for residue,bb,nterm,cterm in zip(struc.residues,struc.backbone,struc.nterm,struc.cterm):


    counter += 1

    # Unpack first atom
    first, resn, resi, chain, x, y, z = residue[0]

    # Extract residue name and atom list
    resn  = resn.strip()
    atoms = [i[0].strip() for i in residue]

    # we have no topology provided, so :
    target = None


    o, r = mapping[resn].do(residue,target,bb,nterm,cterm,options["-nt"])
    out.extend(o)
    raw.extend(r)
new_out = []
for line in out:
    new_line = tuple([line[0], line[1], line[2], line[3], line[4]/10 , line[5]/10, line[6] /10 ])
    new_out.append(new_line)
out = new_out
dev = open(options["-o"].value, "w")
dev.write("Backmapped structure from Sirah to ParmBSC1\n")
dev.write("%5d\n"%len(out))

# Atoms
idx = 1
for atom in out:
    # Regular atom
    nam,res,id,chn,x,y,z = atom
    dev.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(id%1e5,res,nam,idx%1e5,x,y,z))

    idx += 1

# Box
dev.write("  13.03970  13.03970  13.03970   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000")

# Close if we were writing to file
dev.close()















































































































































































