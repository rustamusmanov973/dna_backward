




import glob, os, re, math, sys, random

# Should version this... 130502-11 TAW
# More extensive support for geometric operations


# In the context of this module, a residue is a list of atoms, where each atom/item
# is a list or tuple with at least 7 values:
#
# (atom name (str), residue (str), residue id (int), chain (char), x (float), y (float), z (float))

# Crude mass for weighted averages. No consideration of united atoms.
# This will probably give only minor deviations, while also giving less headache
# We add B with a mass of 32, so
# BB and SC* will have equal weights
def norm2(a):
    return sum([i*i for i in a])

def norm(a):
    return math.sqrt(norm2(a))
def dist(a,b):
    return math.sqrt(norm2([i-j for i,j in zip(a,b)]))

_mass = {'YYYY': 243424324}

# Normalization factor for geometric modifiers (nanometer)
_normfac = 0.125

# Listing of aminoacids
_aminoacids = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "ACE", "NH2",
]
DA_attributes = "C1A C2A N1A N6A".split()
DT_attributes = "C1T O2T N3T O4T".split()

def mirror_DA(x, y, z, residue):
    for i in residue:
        if i[0].strip() == "C1A":
            a1 = tuple([i[4], i[5], i[6]])
        elif i[0].strip() == "N6A":
            a2 = tuple([i[4], i[5], i[6]])
    a0 = tuple([x, y, z])
    r01 = dist(a0, a1)
    r02 = dist(a0, a2)
    r12 = dist(a1, a2)
    p = 0.5*(r01+r02+r12)
    h = (2*math.sqrt(p*(p-r01)*(p-r02)*(p-r12))/r12)
    h1 = math.sqrt(r01*r01-h*h)
    h2 = math.sqrt(r02*r02-h*h)
    ah = tuple([(h1*a2[0] + h2*a1[0])/(h1+h2), (h1*a2[1] + h2*a1[1])/(h1+h2), (h1*a2[2] + h2*a1[2])/(h1+h2)])
    am = [2*ah[0]-a0[0], 2*ah[1]-a0[1], 2*ah[2]-a0[2]]
    return (am)


def mirror_DT(x, y, z, residue):

    for i in residue:
        if i[0].strip() == "C1T":
            a1 = tuple([i[4], i[5], i[6]])
        elif i[0].strip() == "O4T":
            a2 = tuple([i[4], i[5], i[6]])

    a0 = tuple([x, y, z])
    r01 = dist(a0, a1)
    r02 = dist(a0, a2)
    r12 = dist(a1, a2)
    p = 0.5*(r01+r02+r12)
    h = (2*math.sqrt(p*(p-r01)*(p-r02)*(p-r12))/r12)
    h1 = math.sqrt(r01*r01-h*h)
    h2 = math.sqrt(r02*r02-h*h)
    ah = tuple([(h1*a2[0] + h2*a1[0])/(h1+h2), (h1*a2[1] + h2*a1[1])/(h1+h2), (h1*a2[2] + h2*a1[2])/(h1+h2) ])
    am = [2*ah[0]-a0[0], 2*ah[1]-a0[1], 2*ah[2]-a0[2]]
    return(am)

def mirror_DG(x, y, z, residue):
    for i in residue:
        if i[0].strip() == "O6G":
            a1 = tuple([i[4], i[5], i[6]])
        elif i[0].strip() == "C1G":
            a2 = tuple([i[4], i[5], i[6]])
    a0 = tuple([x, y, z])
    r01 = dist(a0, a1)
    r02 = dist(a0, a2)
    r12 = dist(a1, a2)
    p = 0.5*(r01+r02+r12)
    h = (2*math.sqrt(p*(p-r01)*(p-r02)*(p-r12))/r12)
    h1 = math.sqrt(r01*r01-h*h)
    h2 = math.sqrt(r02*r02-h*h)
    ah = tuple([(h1*a2[0] + h2*a1[0])/(h1+h2), (h1*a2[1] + h2*a1[1])/(h1+h2), (h1*a2[2] + h2*a1[2])/(h1+h2) ])
    am = [2*ah[0]-a0[0], 2*ah[1]-a0[1], 2*ah[2]-a0[2]]

    return(am)

def mirror_DC(x, y, z, residue):
    for i in residue:
        if i[0].strip() == "C1C":
            a1 = tuple([i[4], i[5], i[6]])
        elif i[0].strip() == "N4C":
            a2 = tuple([i[4], i[5], i[6]])
    a0 = tuple([x, y, z])
    r01 = dist(a0, a1)
    r02 = dist(a0, a2)
    r12 = dist(a1, a2)
    p = 0.5*(r01+r02+r12)
    h = (2*math.sqrt(p*(p-r01)*(p-r02)*(p-r12))/r12)
    h1 = math.sqrt(r01*r01-h*h)
    h2 = math.sqrt(r02*r02-h*h)
    ah = tuple([(h1*a2[0] + h2*a1[0])/(h1+h2), (h1*a2[1] + h2*a1[1])/(h1+h2), (h1*a2[2] + h2*a1[2])/(h1+h2) ])
    am = [2*ah[0]-a0[0], 2*ah[1]-a0[1], 2*ah[2]-a0[2]]
    return(am)

u = 10000
# Determine average position for a set of atoms
def _average(mas, residue):
    a = []
    for j in mas:

        if j[0] == 1:
            a.append(j[1])
        else:
            if residue[0][1].strip() == 'DA':
                p = tuple([j[1][0], j[1][1], j[1][2], j[1][3]] + mirror_DA(j[1][4], j[1][5], j[1][6], residue))
                a.append(p)
            elif residue[0][1].strip() == 'DT':
                p = tuple([j[1][0], j[1][1], j[1][2], j[1][3]] + mirror_DT(j[1][4], j[1][5], j[1][6], residue))
                a.append(p)
            elif residue[0][1].strip() == 'DG':
                p = tuple([j[1][0], j[1][1], j[1][2], j[1][3]] + mirror_DG(j[1][4], j[1][5], j[1][6], residue))
                a.append(p)
            elif residue[0][1].strip() == 'DC':
                p = tuple([j[1][0], j[1][1], j[1][2], j[1][3]] + mirror_DC(j[1][4], j[1][5], j[1][6], residue))
                a.append(p)
    print("a")
    print(a)
    mxyz = [(_mass.get(i[0][0], 1), i[4], i[5], i[6]) for i in a if i]  # Masses and coordinates
    mw = [sum(i) for i in
        zip(*[(m * x, m * y, m * z, m) for m, x, y, z in mxyz])]
    return [i / mw[3] for i in mw]  # Centre of mass


def kick(x,u):
    return x+(random.random()-0.5)*u

def _vsub(a, b):
    return [i - j for i, j in zip(a, b)]


def _vadd(a, b):
    return [i + j for i, j in zip(a, b)]


def _normalize(a):
    l = math.sqrt(sum([i * i for i in a]))
    return [i / l for i in a]


def _crossprod(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _r(a, kick):
    return a + random.random() * kick - kick / 2


class ResidueMap:
    def __init__(self, target=None, source=None, atoms=None, mod=[], name=""):
        print(atoms)
        if atoms:
            # Setting mapping from an atom list
            # Format is:
            # number, aa, cg beads
            x = [i[1] for i in atoms]
            y = [i[2:] for i in atoms]
            # For those atoms for which no source list is given
            # set the source equal to that of the previous one
            for i in range(len(y)):
                if not y[i]:
                    y[i] = y[i - 1]

        if source:
            y = source

        if target:

            if not atoms:
                x = target

            assert len(x) == len(y)

            # Case of forward mapping: atomistic to martini
            # Initialize dictionary
            d = dict(zip(target, [[] for i in target]))

            # Fill entries
            # The mapping is specified in full both ways, which
            # means that, e.g. for Martini, the result may differ
            # from the original mapping definition. This should
            # add stability, and is required to allow the double
            # mappings used in Martini for some residues.
            for u, v in zip(y, x):
                for j in u:
                    d[j].append(v)

            self.atoms = target
            self.map = d
        else:
            self.atoms = x
            print(x)
            self.map = dict(zip(x, y))
        # Special cases
        self.mod = mod

    def do(self, residue, target=None, coords=False, nterm=False, cterm=False, nt=False, kick=0.05, ):
        # Given a set of source atoms with coordinates
        # return the corresponding list of mapped atoms
        # with suitable starting coordinates.

        # If a target list is given, match every atom against
        # the atoms in the ResidueMap. If an atom is not in
        # the definition, then it is returned with the
        # coordinates of the last atom that was in the list.

        # For amino acids, nterm/cterm will cause extra hydrogens/
        # oxygen to be added at the start/end of the residue.
        # If nt (neutral termini) is true, two hydrogens will be
        # added in stead of three if nterm is true and an additional
        # hydrogen will be added at the end if cterm is true.

        # Unpack first atom
        first, resn, resi, chain, x, y, z = residue[0]
        resn = resn.strip()

        # Check whether a target was supplied
        set_termini = not target

        # Target atoms list
        if target:
            print("sdfsdf"*90)
            # A target atom list can be given as a list of names
            # or as a list of atoms. In the latter case, each element
            # will be a list or tuple and we extract the names.
            if type(target[0]) in (list, tuple):
                target = [i[0].strip() for i in target]
            elif type(target) == str:
                target = target.split()
            else:
                # Make a copy to leave the original list untouched
                target = list(target)

        else:
            target = list(self.atoms)
        # Atoms we have; the source dictionary
        atoms = [i[0].strip() for i in residue]
        have = dict(zip(atoms, residue))
        #for i in have:
         #   print(i)
          #  modified_tuple = tuple([have[i][0], have[i][1], have[i][2], have[i][3], -have[i][4], -have[i][5], -have[i][6]])
           # have["-"+str(i)] = modified_tuple

        # Set array for output; residue to return
        out = []
        five_prime = False
        if "PX" not in have:
            target.remove("P")
            target.remove("O1P")
            target.remove("O2P")
            five_prime = True
        # The target list is leading. Make sure that the atom list matches.
        # So, the actual atom list is built from the target list:
        atomlist = [i for i in target if i in self.atoms]
        # Go over the target particles; the particles we want
        # If we have a target topology, then there may be
        # additional particles to want, especially hydrogens.
        # These will be assigned to the previous heavy atom
        # from the want list.

        for want in atomlist:
            print("have")
            print(have)
            mas = []
            for i in self.map[want]:
                if i[0] == "-":
                    mas.append([-1, have.get(i[1:])])
                else:
                    mas.append([1, have.get(i)])

            if coords:
                got = coords.get(want, _average(mas, residue))


            else:
                #got = _average([have.get(i) for i in self.map[want]])
                got = _average(mas, residue)

            if not got:
                print "Problem determining mapping coordinates for atom %s of residue %s." % (target[0], resn)
                print "atomlist:", atomlist
                print "want:", want, self.map[want]
                print "have:", have.keys()
                print "Bailing out..."
                print target
                sys.exit(1)

            # This logic reads the atom we want together with all atoms
            # that are in the target list, but not in the residue
            # definition in the mapping dictionary, up to the next atom
            # that is in that definition.
            while target and (target[0] == want or target[0] not in self.atoms):
                name = target.pop(0)
                if five_prime == False:
                    out.append((name, resn, resi, chain, got[0], got[1], got[2]))
                else:
                    out.append((name, resn.strip()+"5 ", resi, chain, got[0], got[1], got[2]))





                            # Now add a random kick whenever needed to ensure that no atoms overlap
        for i in range(len(out)):
            for j in range(i):
                if out[i][4] == out[j][4] and out[i][5] == out[j][5] and out[i][6] == out[j][6]:
                    # Equal coordinates: need fix
                    x, y, z = out[i][4:7]
                    kick = 0.90
                    out[i] = out[i][:4] + (_r(x, kick), _r(y, kick), _r(z, kick))

        # Create a lookup dictionary
        atoms = dict(zip([i[0] for i in out], range(len(out))))

        # Create a coordinate dictionary
        # This allows adding dummy particles for increasingly complex operations
        coord = dict([(i[0], (_r(i[4], 1e-5), _r(i[5], 1e-5), _r(i[6], 1e-5))) for i in out])
        # Correct terminal amino acid N/C positions, to increase stability of
        # subsequent corrections.
        if resn in _aminoacids and nterm:
            t = atoms.get("N")

            # Set N at 120 degree angle with respect to CA-C
            b = coord.get("CA")
            c = coord.get("C")

            # Only act if we really have backbone atoms
            if t != None and b != None and c != None:
                u = b[0] - c[0], b[1] - c[1], b[2] - c[2]
                l = math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
                try:
                    v = _normalize(_crossprod(u, (u[0] + 1, u[1], u[2])))
                except ZeroDivisionError:
                    # Oh, so that vector was parallel. Then this must work.
                    v = _normalize(_crossprod(u, (u[0], u[1] + 1, u[2])))

                coord["N"] = (
                b[0] + u[0] / 2 + .866 * l * v[0], b[1] + u[1] / 2 + .866 * l * v[1], b[2] + u[2] / 2 + .866 * l * v[2])
                out[t] = out[t][:4] + coord["N"]

        if resn in _aminoacids and cterm:
            t = atoms.get("C")

            # Set N at 120 degree angle with respect to CA-N
            b = coord.get("CA")
            c = coord.get("N")

            # Only act if we really have backbone atoms
            if t != None and b != None and c != None:
                u = b[0] - c[0], b[1] - c[1], b[2] - c[2]
                l = math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
                try:
                    v = _normalize(_crossprod(u, (u[0] + 1, u[1], u[2])))
                except ZeroDivisionError:
                    # Oh, so that vector was parallel. Then this must work.
                    v = _normalize(_crossprod(u, (u[0], u[1] + 1, u[2])))

                coord["C"] = (
                b[0] + u[0] / 2 + .866 * l * v[0], b[1] + u[1] / 2 + .866 * l * v[1], b[2] + u[2] / 2 + .866 * l * v[2])
                out[t] = out[t][:4] + coord["C"]

        # Before making modifications, save the raw projection results
        raw = [i for i in out]

        # Now add a random kick again whenever needed to ensure that no atoms overlap
        for i in range(len(out)):
            for j in range(i):
                if out[i][4] == out[j][4] and out[i][5] == out[j][5] and out[i][6] == out[j][6]:
                    # Equal coordinates: need fix
                    x, y, z = out[i][4:7]
                    out[i] = out[i][:4] + (_r(x, kick), _r(y, kick), _r(z, kick))

        return out, raw

# These are the modifier tags. They signify a specific
# operation on the atoms listed.
_mods = ("chiral","trans","cis","out")


# These are the default tags. Other tags should be
# coarse grained model names.
_tags = _mods + ("molecule","mapping","atoms")

def _init():
    molecules = []
    mapping   = {}
    cg        = []
    aa        = []
    ff        = []
    mod       = []
    cur       = []
    mol       = []
    tag       = re.compile('^ *\[ *(.*) *\]')

    # Read all .map residue definitions in the module directory
    for filename in glob.glob(os.path.dirname(__file__)+"/*.map"):
        # Default CG model is martini.
        cg_ff = "martini"

        for line in open(filename):

            # Strip leading and trailing spaces
            s = line.strip()

            # Check for directive
            if s.startswith("["):

                # Extract the directive name
                cur = re.findall(tag,s)[0].strip().lower()

                if not cur in _tags: # cur == "martini":
                    cg_ff = cur
                    cg    = []

                # The tag molecule starts a new molecule block
                # The tag mapping starts a new mapping block for a specific force field
                # In both cases, we need to purge the stuff that we have so far and
                # empty the variables, except 'mol'
                if cur in ("molecule","mapping"):
                    # Check whether we have stuff
                    # If so, we purge
                    if aa:
                        for ffi in ff:
                            for m in mol:
                                try:
                                    mapping[(m, 'sirah' ,  ffi  )] = ResidueMap(atoms=aa,mod=mod,name=m)
                                except:
                                    print "Error reading %s to %s mapping for %s (file: %s)."%(cg_ff,ffi,m,filename)

                    # Reset lists
                    aa,ff,mod = [],[],[]

                    # Reset molecule name list if we have a molecule tag
                    if cur == "molecule":
                        mol = []

                continue

            # Remove comments
            s = s.split(';')[0].strip()

            if not s:
                # Skip empty lines
                continue

            elif cur == "molecule":
                mol.extend(s.split())

            elif not cur in _tags:  # cur == "martini":
                # Martini coarse grained beads in topology order
                cg.extend(s.split())

            elif cur == "mapping":
                # Multiple force fields can be specified
                ff.extend(s.split())

            elif cur == "atoms":
                # Atom list for current force field molecule definition
                aa.append(s.split())

            elif cur in _mods:
                # Definitions of modifying operations
                mod.append((cur,s.split()))


    # At the end we may still have rubbish left:
    if aa:
        for ffi in ff:
            for m in mol:
                try:
                    mapping[(m, 'sirah' ,  ffi  )] = ResidueMap(atoms=aa,mod=mod,name=m)
                except:
                    print "Error reading %s to %s mapping for %s (file: %s)."%(cg_ff,ffi,m,filename)


    return mapping


mapping = _init()


def get(target="ParmBSC1",source="sirah"):
    D = dict([(i[0],mapping[i]) for i in mapping.keys() if i[1] == source and i[2] == target])
    print "Residues defined for transformation from %s to %s:"%(source,target)
    print D.keys()
    return D

mapping_mimic = get()




