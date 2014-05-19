'''
Created on May 9, 2014

@author: bertalan
'''
import numpy as np
import kevrekidis as kv
import tomSims
def convert(string):
    out = {}
    lis = string.strip().split('\n')
    for i in range(len(lis)):
        if len(lis[i].split()) == 2:
            key, val = [x.strip() for x in lis[i].split()][:2]  # split on whitespace
            out[key.lower()] = float(val)
    return out

lel = convert('''Kd_pcn     0.3
kdn_pcn     0.02
ksP     1.2
Kdp_pcc     0.1
kdn_pcc     0.02
Kdp_pcn     0.1
Kd_pcc     0.3
CCP     0
V1C     1.2
V4PC     0.2
V3B     1.4
kdn_bnp     0.02
V2PC     0.2
k4     0.4
ksC     3.2
kdn_rc     0.02
vdPCN     1.4
PCCP     0
BCP     0
kdn_in     0.02
Kp_c     1.006
Kdp_bn     0.1
PC     0
Kdp_bc     0.1
Kp_bc     1.006
BC     0
vmB     1.3
vmC     2
Kp_p     1.006
vdPCC     1.4
RN     0
vmR     1.6
BN     0
vmP     2.2
k6     0.4
kdmb     0.02
V2P     0.6
CC     0
vdBC     3
k2     0.4
V1P     9.6
V2C     0.2
V2B     0.2
vdPC     3.4
Kd_bn     0.3
Kd_rc     0.3
kdmc     0.02
h     2
Kd_rn     0.3
k9     0.8
kdmr     0.02
KAR     0.6
kdmp     0.02
KAP     0.6
vdCC     1.4
PCP     0
KIB     2.2
Kd_c     0.3
k1     0.8
Kdp_c     0.1
Kd_p     0.3
vdRC     4.4
vdRN     0.8
ksB     0.32
n_p     2
V3PC     2.4
kdn_cp     0.02
PCN     0
PCC     0
V1B     1.4
KAC     0.6
kdn_pp     0.02
vdIN     1.6
V1PC     2.4
kdn_rn     0.02
n_c     2
vsP     2.4
Kd_in     0.3
vsR     1.6
KmC     0.4
KmB     0.4
k3     0.8
kdn_pcnp     0.02
vsB     1.8
vsC     2.2
k7     1
KmR     0.4
k5     0.8
KmP     0.3
ksR     1.7
vdBN     3
k8     0.2
Kp_bn     1.006
BNP     0
kdn_p     0.02
RC     0.5
kdn_bn     0.02
PCNP     0
kdn_pccp     0.02
k10     0.4
Kp_pcn     1.006
IN     0
Kp_pcc     1.006
kdnc     0.02
Kd_bc     0.3
V4B     0.4
MC     2
MB     8
kdn_Bc     0.02
m     2
Kdp_p     0.1
kdn_bcp     0.02
MP     4
MR     1''')

to = convert('''vsp0    1.02
a    10
b    4
KD    2.5
RT    100
N    100
kD    10
v0    0.5
v1    5
v2    5
vP    1
Ka    2.5
VMK    8
K1    0.01
K2    0.01
CBT    1
CT    1.1
KC    0.3''')

vas = convert('''Ek    -97
THETA    37
gK0    9.7
Vgk    10
Kgk    10
gNa    36
Ena    45
vkk    3.3
Kkk    0.02
nkk    0.1
vvo    0.09
Kvo    4.5
nvo    4.5
v    2
v1    0.0003
betaIP3    0.5
VM2    149.5
K2    5
n    2.2
VM3    400
KR    3
m    6
KA    0.67
p    4.2
kf    0.001
Caex    5
vCa    12.3
Kca    22
nCa    2.2
vKCa    3
KKCa    0.16
nKCa    -1
EL    -29
GABA0    0.2
vGABA    19
KGABA    3
gGABA    12.3
Cl0    1
    
vCl1    15.5
vCl2    19
KCl1    4
KCl2    1
nCl    -0.2
Clex    114.5
vex1    105
Kex1    574.05
nex1    2.5
vex2    4.4
Kex2    1
nex2    -1
Eex    0
Pca    0.05
Pna    0.036
Pcl    0.3
Kex    1
Naex    145
vPK    1.9
KPK    1
npk    -2
VR    0.41
KR    34
vVIP    0.5
KVIP    15
nVIP    1.9
kdVIP    0.5
ndVIP    0.2
VMK    5
KMK    2.9
Vbeta    2
Kbeta    2
CT    1.6
Kc    0.15
KD    0.08
k1    0.45
KAP    0.6
vsp0    1
Cm    5''')

kv.pp((lel, to, vas))
print "lel    to    vas"
for k in lel.keys():
    l = lel[k]
    if k in to:
        t = to[l]
    else:
        t = None
        
    if k in vas:
        v = vas[k]
    else:
        v = None
    print "%s    %s    %s" % (l, t, v)