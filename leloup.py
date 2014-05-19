import numpy as np

import kevrekidis as kv

from parameterHandling import Parameters
class LeloupParameters(Parameters):
    def __init__(self, forcingShape=None, *args, **kwargs):
        super(LeloupParameters, self).__init__(*args, **kwargs)
        self.forcingShape = forcingShape
    def loadDefaults(self):
        from parameters.leloupDefaults import paramsString
        self.addFromList(paramsString.split('\n'))
        

    
def leloupX0params():
    p = LeloupParameters()
    p.loadDefaults()
    initial = {  # From:
                 #  Vasalou, C. & Henson, M. A. A multiscale model to
                 #  investigate circadian rhythmicity of pacemaker neurons in
                 #  the suprachiasmatic nucleus. PLoS Comput. Biol. 6, e1000706
                 #  (2010).
                'BC': 2.41,
                'BCP': 0.48,
                'BN': 1.94,
                'BNP': 0.32,
                'CC': 2.231,
                'CCP': 9.00,
                'IN': 0.05,
#                 'MB': 7.94,
                'MB': 5,
#                 'MC': 2.00,
                'MC': 3.00,
#                 'MP': 2.80,
                'MP': 1.0,
                'PC': 0.40,
                'PCC': 1.26,
                'PCCP': 0.20,
                'PCN': 0.16,
                'PCNP': 0.091,
                'PCP': 0.13
              } 
    p.addFromDict(initial)
    out = (
             p.MP, p.MC, p.MB,  # mRNA
             p.PC, p.CC, p.PCP, p.CCP,  # PER & CRY
             p.PCC, p.PCN, p.PCCP, p.PCNP, # PER-CRY
             p.BC, p.BCP, p.BN, p.BNP,  # BMAL1
             p.IN,  # inactive supercomplex
#              p.MR, p.RC, p.RN  # REV-ERB-\alpha
         )
    return out

    
def leloupX0array():
    p = LeloupParameters()
    p.loadDefaults()
    return np.array(leloupX0params())



def LofT(t, L, Lbase, shape=None):
    if shape is None:
        # Steady light
        out = L - Lbase
        if hasattr(t, '__getitem__'):
            out = [out] * len(t)
        return out
    elif shape.lower() == "sin":
        # sinusoidal forcing
        return (np.sin(t * np.pi / 12.0) + 1)/2 * (L - Lbase) + Lbase
    else:
            # square-wave forcing
        if isinstance(t, np.ndarray):
            return Lbase + (t % 24 < 12).astype(int) * L
        else:
            return Lbase + int(t % 24 < 12) * L



def mRNA(vsP, p, MP, MC, MB, BN, Lt):
    dMPdt = +(Lt + 1) * vsP * BN ** p.n_p / (p.KAP ** p.n_p + BN ** p.n_p) - p.vmP * MP / (p.KmP + MP) - p.kdmp * MP
    dMCdt = +p.vsC * BN ** p.n_c / (p.KAC ** p.n_c + BN ** p.n_c) - p.vmC * MC / (p.KmC + MC) - p.kdmc * MC
    dMBdt = +p.vsB * p.KIB ** p.m / (p.KIB ** p.m + BN ** p.m) - p.vmB * MB / (p.KmB + MB) - p.kdmb * MB # eqn (3)
    return dMPdt, dMCdt, dMBdt


def protein(p, MP, MC, PC, CC, PCP, CCP, PCC):
    dPCdt = +p.ksP * MP - p.V1P * PC / (p.Kp_p + PC) + p.V2P * PCP / (p.Kdp_p + PCP) + p.k4 * PCC - p.k3 * PC * CC - p.kdn_p * PC
    dCCdt = +p.ksC * MC - p.V1C * CC / (p.Kp_c + CC) + p.V2C * CCP / (p.Kdp_c + CCP) + p.k4 * PCC - p.k3 * PC * CC - p.kdnc * CC
    dPCPdt = +p.V1P * PC / (p.Kp_p + PC) - p.V2P * PCP / (p.Kdp_p + PCP) - p.vdPC * PCP / (p.Kd_p + PCP) - p.kdn_pp * PCP
    dCCPdt = +p.V1C * CC / (p.Kp_c + CC) - p.V2C * CCP / (p.Kdp_c + CCP) - p.vdCC * CCP / (p.Kd_c + CCP) - p.kdn_cp * CCP
    return dPCdt, dCCdt, dPCPdt, dCCPdt


def perCry(p, PC, CC, PCC, PCN, PCCP, PCNP, BN, IN):
    dPCCdt = -p.V1PC * PCC / (p.Kp_pcc + PCC) + p.V2PC * PCCP / (p.Kdp_pcc + PCCP) - p.k4 * PCC + p.k3 * PC * CC + p.k2 * PCN - p.k1 * PCC - p.kdn_pcc * PCC
    dPCNdt = -p.V3PC * PCN / (p.Kp_pcn + PCN) + p.V4PC * PCNP / (p.Kdp_pcn + PCNP) - p.k2 * PCN + p.k1 * PCC - p.k7 * BN * PCN + p.k8 * IN - p.kdn_pcn * PCN
    dPCCPdt = +p.V1PC * PCC / (p.Kp_pcc + PCC) - p.V2PC * PCCP / (p.Kdp_pcc + PCCP) - p.vdPCC * PCCP / (p.Kd_pcc + PCCP) - p.kdn_pccp * PCCP
    dPCNPdt = +p.V3PC * PCN / (p.Kp_pcn + PCN) - p.V4PC * PCNP / (p.Kdp_pcn + PCNP) - p.vdPCN * PCNP / (p.Kd_pcn + PCNP) - p.kdn_pcnp * PCNP
    return dPCCdt, dPCNdt, dPCCPdt, dPCNPdt


def bmal1(p, MB, PCN, BC, BCP, BN, BNP, IN):
    dBCdt = +p.ksB * MB - p.V1B * BC / (p.Kp_bc + BC) + p.V2B * BCP / (p.Kdp_bc + BCP) - p.k5 * BC + p.k6 * BN - p.kdn_Bc * BC # kdn_Bc is an irregularity in the parameters file conventions
    dBCPdt = +p.V1B * BC / (p.Kp_bc + BC) - p.V2B * BCP / (p.Kdp_bc + BCP) - p.vdBC * BCP / (p.Kd_bc + BCP) - p.kdn_bcp * BCP # weird that it's not p.kdn_bc. TODO ask about this. Other kdn_*p above are similar.
    dBNdt = -p.V3B * BN / (p.Kp_bn + BN) + p.V4B * BNP / (p.Kdp_bn + BNP) + p.k5 * BC - p.k6 * BN - p.k7 * BN * PCN + p.k8 * IN - p.kdn_bn * BN
    dBNPdt = +p.V3B * BN / (p.Kp_bn + BN) - p.V4B * BNP / (p.Kdp_bn + BNP) - p.vdBN * BNP / (p.Kd_bn + BNP) - p.kdn_bnp * BNP
    return dBCdt, dBCPdt, dBNdt, dBNPdt


def inactiveComplex(p, PCN, BN, IN):
    dINdt = -p.k8 * IN + p.k7 * BN * PCN - p.vdIN * IN / (p.Kd_in + IN) - p.kdn_in * IN
    return dINdt

def diff(X, t, parameters=None, L=.2, Lbase=0.0, vsP=None, forcingShape=None):
    '''X.size = 16
      L is light level'''
                
    if parameters is None:
        p = LeloupParameters()
        p.loadDefaults()
    else:
        p = parameters
    if vsP is None:
        vsP = p.vsP
    
    # genetic model from Vasalou, C., Herzog, E. D. & Henson, M. A.
    #  A multicellular model for intercellular synchronization in circadian
    #  neural networks: supplementary material. Biophys. J. 101, (2011).
    #
    #  SBML code and parameters files are available at:
    #  millar.bio.ed.ac.uk/PEBrown/CircadianModelling/NewModels/NewModels.htm
    
    # Variables
    (
     MP, MC, MB,  # mRNA
     PC, CC, PCP, CCP,  # PER & CRY
     PCC, PCN, PCCP, PCNP, # PER-CRY
     BC, BCP, BN, BNP,  # BMAL1
     IN,  # inactive supercomplex
#      MR, RC, RN  # REV-ERB-\alpha
     ) = X[:16]
    
    
    Lt = LofT(t, L, Lbase, shape=p.forcingShape)
    #mRNA for Per, Cry, and Bmal1
    dMPdt, dMCdt, dMBdt = mRNA(vsP, p, MP, MC, MB, BN, Lt)
#   dMBdt = p.vsB*p.KIB ** p.m/(p.KIB**p.m + RN**p.m) - p.vmB*MB/(p.KmB + MB) - p.kdmb*MB  # eqn (3')
            
    
    # Protein, phosphorylated and unphosphorylated, for PER and CRY in cytosol
    dPCdt, dCCdt, dPCPdt, dCCPdt = protein(p, MP, MC, PC, CC, PCP, CCP, PCC)
    
    # PER-CRY complex in cytosol and nucleus
    dPCCdt, dPCNdt, dPCCPdt, dPCNPdt = perCry(p, PC, CC, PCC, PCN, PCCP, PCNP, BN, IN)
    
    # BMAL1
    dBCdt, dBCPdt, dBNdt, dBNPdt = bmal1(p, MB, PCN, BC, BCP, BN, BNP, IN)
    
    # Inactive PER-CRY/CLOCK-BMAL1 complex in nucleus
    dINdt = inactiveComplex(p, PCN, BN, IN)
    
#     # Addition of REV-ERBalpha; mRNA and protein
#     dMRdt = p.vsR*BN**p.h/(p.KAR**p.h + BN**p.h)\
#             -p.vmR*MR/(p.KmR + MR)\
#             -p.kdmr*MR
#     dRCdt = p.ksR*MR - p.k9*RC + p.k10*RN - p.vdRC*RC/(p.Kd_rc + RC) - p.kdn_rc*RC
#     dRNdt = p.k9*RC - p.k10*RN - p.vdRN*RN/(p.Kd_rn + RN) - p.kdn_rn*RN
    
    return np.array(
                    (
                     dMPdt, dMCdt, dMBdt,  # mRNA
                     dPCdt, dCCdt, dPCPdt, dCCPdt,  # PER & CRY
                     dPCCdt, dPCNdt, dPCCPdt, dPCNPdt, # PER-CRY
                     dBCdt, dBCPdt, dBNdt, dBNPdt,  # BMAL1
                     dINdt,  # inactive supercomplex
#                      dMRdt, dRCdt, dRNdt  # REV-ERB-\alpha
                     )
                    )
    
    
def main(pdict=None, icdict=None, L=.2, Lbase=0.0, tmax=72, filename="leloup_main_result"):
    parameters = LeloupParameters()
    parameters.loadDefaults()
    if pdict is not None:
        parameters.addFromDict(pdict)
        
    def diff_(X, t):
        return diff(X, t, parameters=parameters, L=L, Lbase=Lbase)
    
    IC = leloupX0params()
    if icdict is not None:
        indices = {x.glabel():i for (i,x) in enumerate(IC)}
        IC = list(IC)
        from parameterHandling import FloatParameter
        for key, val in icdict.items():
            IC[indices[key]] = FloatParameter(val)
            IC[indices[key]].slabel(key)
    
    X, T = kv.dynSys.integration.integrate(IC, diff_,
                                           tmin=0, tmax=tmax, giveTime=True, nstep=1e4)
    fig, ax = kv.fa(figsize=(8,5))
    fig.subplots_adjust(bottom=.14)
    X0p = leloupX0params()
    for i in range(3):  # plot just the mRNA
        x = X0p[i]
        ax.plot(T, X[:,i], label=x.glabel())
        ax.axhline(X[0,i])
    twinx = ax.twinx()
    Lt = LofT(T, L, Lbase, parameters.forcingShape)
    twinx.plot(T, Lt, label="light level", color="magenta")
    ax.set_zorder(twinx.get_zorder()+1)  # put ax in front of twinx
    ax.patch.set_visible(False) # hide the 'canvas' 
    ax.legend(loc='upper right')
    twinx.legend(loc='lower right')
    ax.set_ylabel('mRNA levels [nM]')
    twinx.set_ylabel('light level $L$')
    ax.set_xlabel('time $t$ [hr]')
    kv.plotting.saveFig(fig, filename)


def showSpecialICs(X0):
    defaultsList = leloupX0params()
    defaults = {x.glabel(): x.gval() for x in defaultsList}
    this = {x.glabel(): x0 for (x, x0) in zip(defaultsList, X0)}
    return showSpecialValues(defaults, this)


def showSpecialParameters(parameters):
    defaults = LeloupParameters()
    defaults.loadDefaults()
    return showSpecialValues(defaults.__dict__, parameters.__dict__)


def showSpecialValues(normal, this):
    out = []
    for k in normal:
        if normal[k] != this[k]:
            out.append(k + str(this[k]))
    return '-'.join(out)


def varyParameters():
    X0p = leloupX0params()

    def mkplot(X, T, parameters, X0):
        fig, ax = kv.fa()
        for i in range(3):  # plot just the mRNA
            x = X0p[i]
            ax.plot(T, X[:,i], label=x.glabel())
        ax.legend(loc='best')
        filename = []
        for st in showSpecialICs(X0), showSpecialParameters(parameters):
            if st != '':
                filename.append(st)
        filename = '__'.join(filename)
        kv.plotting.saveFig(fig, kv.utils.sanitize(filename), exts=('png',))
    
                
    
    parameters = LeloupParameters()
    parameters.loadDefaults()
    pnames = parameters.__dict__.keys()
    pnames.sort()
    for i in range(len(X0p) + len(pnames)):
        from tomSims import logging
        IC = np.array(X0p)
        parameters = LeloupParameters()
        parameters.loadDefaults()
        if i < len(X0p):
            logging.debug('setting IC %d: %s' % (i, X0p[i].glabel()))
            if IC[i] != 0:
                IC[i] *= 2.231
            else:
                IC[i] += 24.0
        else:
            i -= len(X0p)
            logging.debug('setting parameter %d: %s' % (i, pnames[i]))
            parameters = LeloupParameters()
            parameters.loadDefaults()
            if parameters.__dict__[pnames[i]] != 0:
                parameters.__dict__[pnames[i]] *= 2.231
            else:            
                parameters.__dict__[pnames[i]] += 24.0
            
        X, T = kv.dynSys.integration.integrate(IC, diff, tmin=0, tmax=72, giveTime=True, nstep=3e3)
        mkplot(X, T, parameters, IC)


def showIC():
    X0 = leloupX0params()
    kv.pp({x.glabel():x.gval() for x in X0})


def showParams():
    p = LeloupParameters()
    p.loadDefaults()
    kv.pp(p.__dict__)
#     for x in p.__dict__.values():
#         print x.glabel(), ",", x.gval()
    
    
if __name__=="__main__":
#     main(L=0, Lbase=0, tmax=72)
#     main(L=0, Lbase=0, tmax=72, filename="unforced")
#     main(L=.2, Lbase=0, tmax=72, filename="periodically forced")
    icdict = {x:64 for x in (
#                              'PC', 'CC',
#                              'PCP', 'CCP',
#                              'PCC',
#                              'PCN',
#                              'PCCP', 'PCNP',
#                              'BC',
#                              'BCP',
                             'BN',
#                              'BNP',
#                              'IN'
                             )}
    pdict = {'vmC': 2.0, 'vmP': 1.9, 'vmB': 1.0}
    main(icdict=icdict, pdict=pdict, L=0, Lbase=0, tmax=420, filename="unforced_%s" % str(icdict))
#     main(icdict=icdict, L=0.2, Lbase=0, tmax=72, filename="periodically forced_%s" % str(icdict))
#     main(icdict={"MP":1.6, "MC":1.4, "MB":3.4}, L=.2, Lbase=0, tmax=420, filename="MP1.6_MC1.5_MB3.4")
#     kv.plotting.show()
    
