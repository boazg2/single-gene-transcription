import numpy as np, os, sys
import math
from param_var import *

"""
    Functions used by ../bin/twin.py    
"""

# supercoiling densities
def _sigma(rnap: RNAP, loc='up'):
    return (rnap.Lk[loc] - rnap.Lk0[loc]) / rnap.Lk0[loc]

# Simulation runs
def generate_run(modelP: ModelParam, simuP: SimuParam):

    traj = Trajectory()
    RNAP_list = []
    # DNA-bound RNAPs

    nextevent2iter = {tag: -1 for tag in ["b", "oc", "esc", "gyr", "top"]}
    # next trial for binding, OC formation, promoter escape, and gyrase action
    # -1 to avoid initial True
    # added gyrase to function
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = int(Z / modelP.coarse_g.tau_0)
    
    if modelP.coarse_g.k_top > 0:
        Ztop = np.random.exponential(scale=1/modelP.coarse_g.k_top, size=None)
        nextevent2iter["top"] = int(Ztop / modelP.coarse_g.tau_0)
    else:
        nextevent2iter["top"] = simuP.Niterations
    
    if modelP.coarse_g.kmax_gyr > 0:
        Zgyr = np.random.exponential(scale=1/modelP.coarse_g.kmax_gyr, size=None)
        nextevent2iter["gyr"] = int(Zgyr / modelP.coarse_g.tau_0)
    else:
        nextevent2iter["gyr"] = simuP.Niterations

    write_transcripts_on_the_fly(traj, simuP, header=True)
    write_follow_promoter(traj, None, modelP, simuP, header=True)
    while traj.niter < simuP.Niterations and traj.Ntranscripts < simuP.Ntranscripts_max:
        verbosing(simuP, traj, RNAP_list)
        write_follow_promoter(traj, RNAP_list, modelP, simuP)

        # Binding
        if traj.niter == nextevent2iter["b"]:
            binding_stage(modelP, RNAP_list, traj, nextevent2iter)

        # OC formation
        if traj.niter == nextevent2iter["oc"]:
            oc_formation_stage_sigmoidal(modelP, RNAP_list, traj, nextevent2iter)

        # Promoter escape
        if traj.niter == nextevent2iter["esc"]:
            escape_stage(modelP, RNAP_list, traj)

        # Topoisomerases
        if traj.niter == nextevent2iter["top"]:
            
            
            topoI_sigmoidal(modelP, RNAP_list)
            
            Ztop = np.random.exponential(scale=1/modelP.coarse_g.k_top, size=None)
            nextevent2iter["top"] = traj.niter + max(1, int(Ztop / modelP.coarse_g.tau_0))

        if traj.niter == nextevent2iter["gyr"]:
            
            #call gyrase procession
            gyrase_processive(modelP, RNAP_list)
            
            #calculat time for next binding attempt
            Zgyr = np.random.exponential(scale=1/modelP.coarse_g.kmax_gyr, size=None)
            nextevent2iter["gyr"] = traj.niter + max(1, int(Zgyr / modelP.coarse_g.tau_0))
            

        # Elongation
        elongation_stage(modelP, RNAP_list)

        # Termination
        termination_stage(modelP, simuP, RNAP_list, traj)

        traj.niter += 1
        traj.time = traj.niter * modelP.coarse_g.tau_0

    return


# Transcription stages
def binding_stage(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter):
    """Binding of the RNAP at the promoter"""

    if (
        not RNAP_list
        or RNAP_list[-1].X
        > modelP.gene.rnap_xi + modelP.rnap.excluded_length_elongation
    ):
        # promoter is free => binding occurs!

        if traj.ref_binding_time > 0:
            traj.b2b_times["mean"] = (
                traj.b2b_times["n"] * traj.b2b_times["mean"]
                + traj.time
                - traj.ref_binding_time
            ) / (traj.b2b_times["n"] + 1)
            traj.b2b_times["n"] += 1

        traj.ref_binding_time = traj.time

        # SETTING UP NEW RNAP
        rnap = RNAP()

        rnap.tb = traj.time
        rnap.t_elongating = False
        rnap.X = modelP.gene.rnap_xi
        rnap.Lk0['up'] = modelP.gene.Lk0_rnap_xi

        # SIGMA PROPERTIES
        if not RNAP_list:
            # a single RNAP => topological properties dictated by the entire topological domain
            rnap.sigma['up'] = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
            ) / modelP.gene.Lk0_domain
            rnap.Lk0['down'] = modelP.gene.Lk0_domain - modelP.gene.Lk0_rnap_xi
        else:
            # multiple RNAP => topological properties dictated by the immediate downstream RNAP
            rnap.sigma['up'] = RNAP_list[-1].sigma['up']
            rnap.Lk0['down'] = RNAP_list[-1].Lk0['up'] - rnap.Lk0['up']

        rnap.Lk['up'] = (1 + rnap.sigma['up']) * rnap.Lk0['up']

        rnap.sigma['down'] = rnap.sigma['up']
        rnap.Lk['down'] = (1 + rnap.sigma['down']) * rnap.Lk0['down']

        # ADDING THE NEW RNAP
        RNAP_list.append(rnap)

        # NEXT STAGE: OC FORMATION
        Z = np.random.exponential(scale=modelP.promoter.ko_s, size=None)
        nextevent2iter["oc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1

    # NEXT BINDING TRIAL (INDEPENDENT WHETHER BINDING HAS OCCURED OR NOT)
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = np.max(
        (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
    )  # at least niter + 1

    return

def oc_formation_stage_sigmoidal(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter):
    """OC formation if sigma <= threshold"""
    
    if(np.random.binomial(1, 1 / (1 + np.exp(-(RNAP_list[-1].sigma['up'] - modelP.promoter.sigma_o) / modelP.promoter.o_width)))):
        
        # NEXT OC FORMATION TRIAL
        Z = np.random.exponential(scale=modelP.promoter.ko_s, size=None)
        nextevent2iter["oc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1
        
    else:
        # sigma is below threshold: OC formation occurs!

        RNAP_list[-1].tocf = traj.time

        traj.ocf_times["mean"] = (
            traj.ocf_times["n"] * traj.ocf_times["mean"]
            + RNAP_list[-1].tocf
            - RNAP_list[-1].tb
        ) / (traj.ocf_times["n"] + 1)
        traj.ocf_times["n"] += 1

        # NEXT STAGE: ESCAPE EVENT
        Z = np.random.exponential(scale=modelP.promoter.ke_s, size=None)
        nextevent2iter["esc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1

    return


def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
    """promoter escape => RNAP is now in elongating mode"""

    RNAP_list[-1].t_elongating = True
    RNAP_list[-1].tesc = traj.time

    traj.esc_times["mean"] = (
        traj.esc_times["n"] * traj.esc_times["mean"]
        + RNAP_list[-1].tesc
        - RNAP_list[-1].tocf
    ) / (traj.esc_times["n"] + 1)
    traj.esc_times["n"] += 1

    traj.initiation_times["mean"] = (
        traj.initiation_times["n"] * traj.initiation_times["mean"]
        + RNAP_list[-1].tesc
        - RNAP_list[-1].tb
    ) / (traj.initiation_times["n"] + 1)
    traj.initiation_times["n"] += 1

    # UPDATING DOWNSTREAM RNAP
    # change Lk because the new elongating RNA becomes a barrier
    # one can actually check the conservation of the Lk
    if len(RNAP_list) > 1:
        RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
        RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']

    return
    
def topoI_sigmoidal(modelP: ModelParam, RNAP_list):
    """topoI activity under processive model"""

    if RNAP_list:
        topoI_sigmoidal_RNAPpresent(modelP, RNAP_list)
    else:
        topoI_sigmoidal_RNAPabsent(modelP)
        
    return

def topoI_sigmoidal_RNAPpresent(modelP: ModelParam, RNAP_list):
    #draw upstream vs downstream (1 is downstream, 0 is upstream)
    location = np.random.binomial(1, modelP.gene.Ldown / (modelP.gene.tss + modelP.gene.Ldown))
    
    #set sigma
    if(location):
        sigma = RNAP_list[0].sigma['down']
    else:
        sigma = RNAP_list[-1].sigma['up']
    #check for successful binding
    if(np.random.binomial(1, 1 / (1 + np.exp(-(sigma - modelP.topoI.sigma_active) / modelP.topoI.topoI_width)))):
       return
    
    #split into five cases
    if(location):
        if(not RNAP_list[0].t_elongating):
            #case 1: 1 RNAP, not elongating (if downstreammost RNAP not elongating, it must be at promotor, so it is also the upstreammost RNAP)
            
            #calculate variables
            domain_length_topo = modelP.gene.L_domain
            sigma = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
            
            #update gene
            modelP.gene.Lk_domain += 1
            
            #update RNAP
            RNAP_list[0].sigma['up'] = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
            RNAP_list[0].Lk['up'] = (1 + RNAP_list[0].sigma['up']) * RNAP_list[0].Lk0['up']
            RNAP_list[0].sigma['down'] = RNAP_list[
                0
            ].sigma['up']  # because RNAP is not a barrier
            RNAP_list[0].Lk['down'] = (1 + RNAP_list[0].sigma['down']) * RNAP_list[
                0
            ].Lk0['down']
            
        else:
            #case 2: acts downstream of elonating RNAP
            
            #calculate variables
            domain_length_topo = RNAP_list[0].Lk0['down'] * modelP.dna.n
            sigma = RNAP_list[0].sigma['down']
            
            #update domain
            modelP.gene.Lk_domain += 1

            #update RNAP 0
            RNAP_list[0].Lk['down'] += 1
            RNAP_list[0].sigma['down'] = _sigma(RNAP_list[0], 'down')
            
    else:
        if(not RNAP_list[-1].t_elongating):
            if(len(RNAP_list) < 2):
                #case 4 (which is just case 1 again of 1 RNAP w/o elongation)
                
                #calculate variables
                domain_length_topo = modelP.gene.L_domain
                sigma = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                    ) / modelP.gene.Lk0_domain
                
                #update gene
                modelP.gene.Lk_domain += 1
                
                #update RNAP
                RNAP_list[0].sigma['up'] = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                    ) / modelP.gene.Lk0_domain
                RNAP_list[0].Lk['up'] = (1 + RNAP_list[0].sigma['up']) * RNAP_list[0].Lk0['up']
                RNAP_list[0].sigma['down'] = RNAP_list[
                    0
                ].sigma['up']  # because RNAP is not a barrier
                RNAP_list[0].Lk['down'] = (1 + RNAP_list[0].sigma['down']) * RNAP_list[
                    0
                ].Lk0['down']
            else:
                #case 5: acts upstream of several RNAPs, of which the RNAP -1 is not elongating
                
                #calculate variables
                domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
                sigma = RNAP_list[-2].sigma['up']
                
                #update gene
                modelP.gene.Lk_domain += 1
                
                #update RNAP -2
                RNAP_list[-2].Lk['up'] += 1
                RNAP_list[-2].sigma['up'] = _sigma(RNAP_list[-2], 'up')
                
                #update RNAP -1
                RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']
                RNAP_list[-1].sigma['down'] = RNAP_list[-2].sigma['up']
                RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']
                RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[-1].Lk0['down']
                
        else:
            #case 3: Acts upstream of elongating RNAP
            
            #calculate variables
            domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
            sigma = RNAP_list[-1].sigma['up']
            
            #update gene
            modelP.gene.Lk_domain += 1
            
            #update RNAP -1
            RNAP_list[-1].Lk['up'] += 1
            RNAP_list[-1].sigma['up'] = _sigma(RNAP_list[-1], 'up')      
            
    return

def topoI_sigmoidal_RNAPabsent(modelP: ModelParam):
    """For when topoI acts on whole domain"""
    
    sigma = (
        modelP.gene.Lk_domain - modelP.gene.Lk0_domain
    ) / modelP.gene.Lk0_domain
    
    #check if binding occurs
    if(np.random.binomial(1, 1 / (1 + np.exp(-(sigma - modelP.topoI.sigma_active) / modelP.topoI.topoI_width)))):
       return
    
    modelP.gene.Lk_domain += 1

    return
    

def gyrase_processive(modelP: ModelParam, RNAP_list):
    """gyrase activity under processive model"""

    if RNAP_list:
        gyrase_processive_RNAPpresent(modelP, RNAP_list)
    else:
        gyrase_processive_RNAPabsent(modelP)

    return

def gyrase_processive_RNAPpresent(modelP: ModelParam, RNAP_list):
    #draw upstream vs downstream (1 is downstream, 0 is upstream)
    location = np.random.binomial(1, modelP.gene.Ldown / (modelP.gene.tss + modelP.gene.Ldown))
    
    #set sigma
    if(location):
        sigma = RNAP_list[0].sigma['down']
    else:
        sigma = RNAP_list[-1].sigma['up']
    #check for successful binding
    if(np.random.binomial(1, 1 / (1 + np.exp((sigma - modelP.gyrase.sigma_active) / modelP.gyrase.gyrase_width)))):
       return
    
    #split into six cases
    if(location):
        if(not RNAP_list[0].t_elongating):
            #case 1: 1 RNAP, not elongating (if downstreammost RNAP not elongating, it must be at promotor, so it is also the upstreammost RNAP)
            
            #calculate variables
            domain_length_topo = modelP.gene.L_domain
            sigma = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
            DLk = DLk_Gyrase_processive(domain_length_topo, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)
            
            #update gene
            modelP.gene.Lk_domain += DLk
            
            #update RNAP
            RNAP_list[0].sigma['up'] = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
            RNAP_list[0].Lk['up'] = (1 + RNAP_list[0].sigma['up']) * RNAP_list[0].Lk0['up']
            RNAP_list[0].sigma['down'] = RNAP_list[
                0
            ].sigma['up']  # because RNAP is not a barrier
            RNAP_list[0].Lk['down'] = (1 + RNAP_list[0].sigma['down']) * RNAP_list[
                0
            ].Lk0['down']
            
        else:
            #case 2: acts downstream of elonating RNAP
            
            #calculate variables
            domain_length_topo = RNAP_list[0].Lk0['down'] * modelP.dna.n
            sigma = RNAP_list[0].sigma['down']
            DLk = DLk_Gyrase_processive(domain_length_topo, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)
            
            #update domain
            modelP.gene.Lk_domain += DLk

            #update RNAP 0
            RNAP_list[0].Lk['down'] += DLk
            RNAP_list[0].sigma['down'] = _sigma(RNAP_list[0], 'down')
            
    else:
        if(not RNAP_list[-1].t_elongating):
            if(len(RNAP_list) < 2):
                #case 4 (which is just case 1 again of 1 RNAP w/o elongation)
                
                #calculate variables
                domain_length_topo = modelP.gene.L_domain
                sigma = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                    ) / modelP.gene.Lk0_domain
                DLk = DLk_Gyrase_processive(domain_length_topo, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)
                
                #update gene
                modelP.gene.Lk_domain += DLk
                
                #update RNAP
                RNAP_list[0].sigma['up'] = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                    ) / modelP.gene.Lk0_domain
                RNAP_list[0].Lk['up'] = (1 + RNAP_list[0].sigma['up']) * RNAP_list[0].Lk0['up']
                RNAP_list[0].sigma['down'] = RNAP_list[
                    0
                ].sigma['up']  # because RNAP is not a barrier
                RNAP_list[0].Lk['down'] = (1 + RNAP_list[0].sigma['down']) * RNAP_list[
                    0
                ].Lk0['down']
            else:
                #case 5: acts upstream of several RNAPs, of which the RNAP -1 is not elongating
                
                #calculate variables
                domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
                sigma = RNAP_list[-2].sigma['up']
                DLk = DLk_Gyrase_processive(domain_length_topo, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)
                
                #update gene
                modelP.gene.Lk_domain += DLk
                
                #update RNAP -2
                RNAP_list[-2].Lk['up'] += DLk
                RNAP_list[-2].sigma['up'] = _sigma(RNAP_list[-2], 'up')
                
                #update RNAP -1
                RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']
                RNAP_list[-1].sigma['down'] = RNAP_list[-2].sigma['up']
                RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']
                RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[-1].Lk0['down']
                
        else:
            #case 3: Acts upstream of elongating RNAP
            
            #calculate variables
            domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
            sigma = RNAP_list[-1].sigma['up']
            DLk = DLk_Gyrase_processive(domain_length_topo, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)
            
            #update gene
            modelP.gene.Lk_domain += DLk
            
            #update RNAP -1
            RNAP_list[-1].Lk['up'] += DLk
            RNAP_list[-1].sigma['up'] = _sigma(RNAP_list[-1], 'up')      
            
    return

def gyrase_processive_RNAPabsent(modelP: ModelParam):
    
    """For when processive gyrase acts on whole domain"""
    # gyrase
    sigma = (
        modelP.gene.Lk_domain - modelP.gene.Lk0_domain
    ) / modelP.gene.Lk0_domain
    
    #check if binding occurs
    if(np.random.binomial(1, 1 / (1 + np.exp((sigma - modelP.gyrase.sigma_active) / modelP.gyrase.gyrase_width)))):
       return
    
    modelP.gene.Lk_domain += DLk_Gyrase_processive(modelP.gene.L_domain, sigma, modelP.gyrase.gyrase_expected_actions, modelP.gyrase.gyrase_sigma_procession)

    return
    
def DLk_Gyrase_processive(domain_length_topo, sigma, expectation, sigma_bound):
    """Change in linking number from processive gyrase actions"""
    Lk0 = domain_length_topo / 10.5
    Lk = Lk0 * (1 + sigma)
    Lk_bound = Lk0 * (1 + sigma_bound)
    
    max_actions = max(0, 1 + math.floor((Lk - Lk_bound) / 2))
    
    drawn_actions = np.random.geometric(1 / (expectation + 1)) 
    #We need to add one to the expectation to account that the last cycle attempt doesn't result in an action
    
    DLk = -2 * min(max_actions, drawn_actions)
    
    return DLk
    
    
def elongation_stage(modelP: ModelParam, RNAP_list):
    """RNAPs translocation"""

    if RNAP_list:
        ix_ = np.arange(len(RNAP_list))
        np.random.shuffle(ix_)

        for ix_rnap in ix_:
            RNAP_translocation(ix_rnap, RNAP_list, modelP)

    return


# Translocations and their topological consequences
def RNAP_translocation(ix_rnap, RNAP_list, modelP: ModelParam):

    if (
        not RNAP_list[ix_rnap].t_elongating
        or not RNAP_list[ix_rnap].sigma['up'] >= modelP.rnap.sigma_stall
        or not RNAP_list[ix_rnap].sigma['down'] <= np.abs(modelP.rnap.sigma_stall)
    ):
        # not elongating or beyond sigma_stall = no translocation
        return

    # Rem1: linking numbers do not change, supercoiling density does
    # Rem2: no constrain on one RNAP overtaking another one (supercoiling constraints do the job)

    RNAP_list[ix_rnap].X += modelP.coarse_g.dx
    RNAP_list[ix_rnap].Lk0['up'] += modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma['up'] = _sigma(RNAP_list[ix_rnap], 'up')

    RNAP_list[ix_rnap].Lk0['down'] -= modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma['down'] = _sigma(RNAP_list[ix_rnap], 'down')

    # UPDATING DOWNSTREAM RNAP
    if ix_rnap > 0:
        RNAP_list[ix_rnap - 1].Lk0['up'] -= modelP.coarse_g.dLk
        RNAP_list[ix_rnap - 1].sigma['up'] = _sigma(RNAP_list[ix_rnap - 1], 'up')

    # UPDATING UPSTREAM RNAP
    if ix_rnap < len(RNAP_list) - 1:
        if RNAP_list[ix_rnap + 1].t_elongating:
            # the RNAP is elongating
            RNAP_list[ix_rnap + 1].Lk0['down'] += modelP.coarse_g.dLk
            RNAP_list[ix_rnap + 1].sigma['down'] = _sigma(
                RNAP_list[ix_rnap + 1], 'down'
            )
        elif (
            ix_rnap == len(RNAP_list) - 2
            and not RNAP_list[len(RNAP_list) - 1].t_elongating
        ):
            # the RNAP is NON-elongating
            RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']
            RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']

            RNAP_list[-1].Lk0['down'] += modelP.coarse_g.dLk
            RNAP_list[-1].sigma['down'] = RNAP_list[-2].sigma['up']
            RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[
                -1
            ].Lk0['down']

    return


def termination_stage(modelP: ModelParam, simuP: SimuParam, RNAP_list, traj: Trajectory):
    """Termination stage: transcript production by the most downstream RNAP"""

    if RNAP_list and RNAP_list[0].X >= modelP.gene.term:
        traj.Ntranscripts += 1

        # ELONGATION TIME
        traj.elongation_times["mean"] = (
            traj.elongation_times["n"] * traj.elongation_times["mean"]
            + traj.time
            - RNAP_list[0].tesc
        ) / (traj.elongation_times["n"] + 1)
        traj.elongation_times["n"] += 1

        # PRODUCTION TIME
        if traj.time_last_prod > 0:
            traj.prod_times["mean"] = (
                traj.prod_times["n"] * traj.prod_times["mean"]
                + traj.time
                - traj.time_last_prod
            ) / (traj.prod_times["n"] + 1)
            traj.prod_times["n"] += 1
        traj.time_last_prod = traj.time

        # TRANSCRIPT TERMINATION
        # 1. We update the upstream RNAP, if it exists        
        if len(RNAP_list) > 1: # at least 2 RNAPs: upating of the upstream RNAP (index = 1)
            RNAP_list[1].Lk0['down'] += RNAP_list[0].Lk0['down']
            if RNAP_list[1].t_elongating:
                # only downstream properties are updated                
                RNAP_list[1].Lk['down'] += RNAP_list[0].Lk['down']
                RNAP_list[1].sigma['down'] = _sigma(RNAP_list[1], 'down')
            else:
                # both upstrean and dowsntream properties are updtaed from sigma of the domain
                RNAP_list[1].sigma['up'] = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
                RNAP_list[1].Lk['up'] = (1 + RNAP_list[1].sigma['up']) * RNAP_list[1].Lk0['up']

                RNAP_list[1].sigma['down'] = RNAP_list[1].sigma['up']
                RNAP_list[1].Lk['down'] = (1 + RNAP_list[1].sigma['down']) * RNAP_list[1].Lk0['down']
        # 2: We remove the RNAP
        del RNAP_list[0]

        if not traj.Ntranscripts % simuP.Nevery_transcripts:
            write_transcripts_on_the_fly(traj, simuP)

    return

# I/O
def verbosing(simuP: SimuParam, traj: Trajectory, RNAP_list):
    """some info to STDOUT (verbosing mode, e.g., for debugging)"""

    if simuP.verbose and not traj.niter % simuP.verbose_every:
        if not len(RNAP_list):
            print(traj.niter, traj.Ntranscripts, len(RNAP_list), end="\r")
        else:
            if len(RNAP_list) == 1:
                print(
                    traj.niter,
                    traj.Ntranscripts,
                    len(RNAP_list),
                    RNAP_list[-1].sigma['up'],
                    RNAP_list[-1].sigma['down'],
                    end="\r",
                )
            else:
                print(
                    traj.niter,
                    traj.Ntranscripts,
                    len(RNAP_list),
                    RNAP_list[-1].sigma['up'],
                    RNAP_list[-1].sigma['down'],
                    RNAP_list[-1].X,
                    RNAP_list[-1].t_elongating,
                    RNAP_list[-2].sigma['up'],
                    RNAP_list[-2].sigma['down'],
                    RNAP_list[-2].X,
                    RNAP_list[-2].t_elongating,
                    RNAP_list[0].sigma['down'],
                    RNAP_list[0].X,
                    end="\r",
                )


def write_follow_promoter(traj: Trajectory, RNAP_list, modelP: ModelParam, simuP: SimuParam, header=False):
    """promoter properties"""
    #Modified: Added Lk_domain output
    fi = simuP.fo_out + "/traj_promoter.txt"

    if header:
        with open(fi, "w") as out:
            out.write("time\tsigma_avg\tsigma_up\tsigma_gene\tsigma_down\tt_upRNAPelong\n")
        return

    with open(fi, "a") as out:
        avg_sigma = (
            modelP.gene.Lk_domain - modelP.gene.Lk0_domain
        ) / modelP.gene.Lk0_domain
        if len(RNAP_list) > 0:
            sigma = RNAP_list[-1].sigma['up']
            down_sigma = RNAP_list[0].sigma['down']
        else:
            sigma = avg_sigma
            down_sigma = avg_sigma
        gene_sigma = (avg_sigma * modelP.gene.L_domain - modelP.gene.tss * sigma - modelP.gene.Ldown * down_sigma) / modelP.gene.L
        out.write(
            "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%d\n"
            % (
                str(traj.time),
                avg_sigma,
                sigma,
                gene_sigma,
                down_sigma,
                (len(RNAP_list) > 0 and RNAP_list[-1].t_elongating),
            )
        )
    return


def write_transcripts_on_the_fly(traj: Trajectory, simuP: SimuParam, header=False):
    """mean properties for every stage of the transcription process"""

    fi = simuP.fo_out + "/mean_properties.txt"

    if header:
        with open(fi, "w") as out:
            out.write(
                "transcripts_nb\ttime\tprod_rate\tmean_prod_time\tmean_bind_time\tmean_ocf_time\tmean_esc_time\tmean_init_time\tmean_elong_time\n"
            )
        return

    with open(fi, "a") as out:
        prod_rate = traj.Ntranscripts / traj.time
        out.write(
            "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
            % (
                traj.Ntranscripts,
                traj.time,
                prod_rate,
                traj.prod_times["mean"],
                traj.b2b_times["mean"],
                traj.ocf_times["mean"],
                traj.esc_times["mean"],
                traj.initiation_times["mean"],
                traj.elongation_times["mean"],
            )
        )

    return

# I/O
def output_variables(cmd, modelP: ModelParam, simuP: SimuParam):
    """writing out parameters and variables"""

    with open(simuP.fo_out + "/param_var.txt", "w") as out:
        out.write(cmd + "\n")

        out.write("\n")
        out.write("###########\n")
        out.write("# General #\n")
        out.write("###########\n")
        dicto = modelP.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("########\n")
        out.write("# Gene #\n")
        out.write("########\n")
        dicto = modelP.gene.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("############\n")
        out.write("# Promoter #\n")
        out.write("############\n")
        dicto = modelP.promoter.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("########\n")
        out.write("# RNAP #\n")
        out.write("########\n")
        dicto = modelP.rnap.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("#########\n")
        out.write("# TopoI #\n")
        out.write("#########\n")
        dicto = modelP.topoI.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("##########\n")
        out.write("# Gyrase #\n")
        out.write("##########\n")
        dicto = modelP.gyrase.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("#######\n")
        out.write("# DNA #\n")
        out.write("#######\n")
        dicto = modelP.dna.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("###################\n")
        out.write("# Coarse graining #\n")
        out.write("###################\n")
        dicto = modelP.coarse_g.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("##############\n")
        out.write("# Statistics #\n")
        out.write("##############\n")
        dicto = simuP.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")
