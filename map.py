import numpy as np
indata = open("coord.txt")
mas = indata.readlines()
cor = [1,2,3,4]
all_cor = []
for line in mas:
    cor = [1,2,3,4]
    line = line.strip()
    inte = line.split(" ")
    cor[0] = inte[0].split(".")[1]
    cor[1], cor[2], cor[3] = float(inte[2]), float(inte[3]), float(inte[4])
    all_cor.append(cor)
    res = line[1:3]
cor_A = {"P":"PX", "C5'":"C5X", "C1'":"C1A", "C2":"C2A", "N1":"N1A", "N6":"N6A"}
cor_T = {"P":"PX", "C5'":"C5X", "C1'":"C1T", "O2":"O2T", "N3":"N3T", "O4":"O4T"}
cor_G = {"P":"PX", "C5'":"C5X", "C1'":"C1G", "N2":"N2G", "N1":"N1G", "O6": "O6G"}
cor_C = {"P":"PX", "C5'":"C5X", "C1'":"C1C", "O2":"O2C", "N3":"N3C", "N4":"N4C"}

a1 = tuple(all_cor[0][1:4])
a2 = tuple(all_cor[1][1:4])
a3 = tuple(all_cor[2][1:4])
a0 = tuple(all_cor[3][1:4])
x1, x2, x3 = zip(a1,a2,a3)[0]
y1, y2, y3  = zip(a1,a2,a3)[1]
z1, z2, z3 = zip(a1,a2,a3)[2]
x, y, z = a0

a = np.array([[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]])
b = np.array([x, y, z])
pre = np.linalg.solve(a,b)
minim = min(abs(pre))
reduced = [round(float(i)/float(minim), 1) for i in pre]
quot = [int(i*10) for i in reduced]

out = open("lal_out.txt", "w")
ou = ""
for i in range(3):
    if res == "DA":

        if quot[i]>0:
            ou += ((cor_A[all_cor[i][0]]+" ")*quot[i])
        else:
            ou +=(("-"+cor_A[all_cor[i][0]]+" ")*(-quot[i]))
    if res == "DT":

        if quot[i]>0:
            ou += ((cor_T[all_cor[i][0]]+" ")*quot[i])
        else:
            ou +=(("-"+cor_T[all_cor[i][0]]+" ")*(-quot[i]))
    if res == "DG":
        if quot[i]>0:
            ou += ((cor_G[all_cor[i][0]]+" ")*quot[i])
        else:
            ou +=(("-"+cor_G[all_cor[i][0]]+" ")*(-quot[i]))
    if res == "DC":

        if quot[i]>0:
            ou += ((cor_C[all_cor[i][0]]+" ")*quot[i])
        else:
            ou +=(("-"+cor_C[all_cor[i][0]]+" ")*(-quot[i]))

print(ou)

pro = ou.split()
mir = pro[-1]
co = []
for i in pro:
    if i == mir:
        co.append("-"+mir)
    else:
        co.append(i)
ou2 = ""
for i in co:
    ou2+= i+" "
print(ou2)
