import argparse
import os, sys
import numpy as np
import shutil

"""
    Parameters and variables used by ../bin/twin.py    
"""


def parsing_cmd():

    parser = argparse.ArgumentParser(
        description="Run the twin transcriptional-loop model in presence of topoisomerases",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # POSITIONAL ARGUMENTS
    parser.add_argument(
        "output_folder", help="results output folder (-f to force overwriting)"
    )

    # OPTIONAL ARGUMENTS
    parser.add_argument(
        "-f", action="store_true", help="overwrite results output folder\n"
    )
    parser.add_argument(
        "-promfollow", action="store_true", help="following promoter status\n"
    )

    # GENE
    parser.add_argument(
        "-tss",
        metavar="X",
        dest="gene_tss",
        type=int,
        help="distance between TSS and upstream barrier (in bp)",
        default=55,
    )
    parser.add_argument(
        "-L",
        metavar="X",
        dest="gene_L",
        type=int,
        help="gene length (in bp)",
        default=900,
    )
    parser.add_argument(
        "-Ld",
        metavar="X",
        dest="gene_Ldown",
        type=int,
        help="distance seperating gene termination STOP and downstream barrier (in bp)",
        default=320,
    )

    # PROMOTER
    parser.add_argument(
        "-kb",
        metavar="X",
        dest="promoter_kb",
        type=float,
        help="promoter binding rate (per s)",
        default=1,
    )
    parser.add_argument(
        "-ko",
        metavar="X",
        dest="promoter_ko",
        type=float,
        help="OC formation rate (per s)",
        default=1,
    )
    parser.add_argument(
        "-so",
        metavar="X",
        dest="promoter_sigma_o",
        type=float,
        help="supercoiling density threshold for OC formation (no unit)",
        default=-0.05,
    )
    parser.add_argument(
        "-ob",
        metavar="X",
        dest="promoter_sigma_width",
        type=float,
        help="width of OC formation curve (no unit)",
        default=0.01,
    )
    parser.add_argument(
        "-ke",
        metavar="X",
        dest="promoter_ke",
        type=float,
        help="OC escape rate (per s)",
        default=1,
    )

    # RNAP
    parser.add_argument(
        "-vm",
        metavar="X",
        dest="rnap_vm",
        type=float,
        help="RNAP maximum translocation speed (in bp per s)",
        default=25,
    )
    parser.add_argument(
        "-G",
        metavar="X",
        dest="rnap_torque_stall",
        type=float,
        help="RNAP stalling torque (in pN.nm)",
        default=18.5,
    )
    parser.add_argument(
        "-rel",
        metavar="X",
        dest="rnap_excluded_length",
        type=float,
        help="RNAP excluded length (in bp)",
        default=30,
    )

    # TOPOI
    parser.add_argument(
        "-lansT",
        metavar="X",
        dest="topoI_lambda_ns",
        type=float,
        help="non-specific rate of upstream TopoI activity (per bp per s)",
        default=1e-4,
    )
    parser.add_argument(
        "-tsa",
        metavar="X",
        dest="topoI_sigma_active",
        type=float,
        help="Midpoint of topoI SC dependency (no unit)",
        default=-0.05,
    )
    parser.add_argument( #added v1.2 on 2024-01-11
        "-tb",
        metavar="X",
        dest="topoI_width",
        type=float,
        help="width of topoI SC dependency (no units)",
        default=0.01,
    )
    # GYRASE
    parser.add_argument(
        "-lansG",
        metavar="X",
        dest="gyrase_lambda_ns",
        type=float,
        help="non-specific rate of downstream gyrase activity (per bp per s)",
        default=1e-4,
    )
    parser.add_argument( #added v1.2 on 2024-01-11
        "-gsa",
        metavar="X",
        dest="gyrase_sigma_active",
        type=float,
        help="Midpoint of gyrase SC dependency (no unit)",
        default=None,
    )
    parser.add_argument( #added v1.2 on 2024-01-11
        "-gb",
        metavar="X",
        dest="gyrase_width",
        type=float,
        help="width of gyrase SC dependency (no units)",
        default=0.01,
    )
    parser.add_argument( #added v1.2 on 2024-01-11
        "-gsp",
        metavar="X",
        dest="gyrase_sigma_procession",
        type=float,
        help="SC above which gyrase catalytic cycle can proceed (no units)",
        default=None,
    )
    parser.add_argument( #added v1.2 on 2024-01-11
        "-mg",
        metavar="X",
        dest="gyrase_expected_actions",
        type=float,
        help="Expected number of actions from gyrase (no units)",
        default=4,
    )
    
    # DNA
    parser.add_argument(
        "-A",
        metavar="X",
        dest="dna_A",
        type=float,
        help="proportionality constant between torque and sigma for unstructured DNA (in pN.nm)",
        default=300,
    )

    # SIMULATION
    parser.add_argument(
        "-Ni",
        metavar="X",
        dest="simup_Niterations",
        type=int,
        help="maximum number of iterations",
        default=int(1e10),
    )
    parser.add_argument(
        "-Nt",
        metavar="X",
        dest="simup_Ntranscripts_max",
        type=int,
        help="maximum number of transcripts",
        default=int(1e4),
    )
    parser.add_argument(
        "-Net",
        metavar="X",
        dest="simup_Nevery_transcripts",
        type=int,
        help="writing out results every X",
        default=int(1e3),
    )

    # RESOLUTION
    parser.add_argument(
        "-dx",
        metavar="X",
        dest="coarse_g_dx",
        type=float,
        help="resolution, i.e., RNAP translocation step (in bp)",
        default=1,
    )

    # TESTING/DEBUGGING
    parser.add_argument(
        "-verbose",
        "-v",
        action="store_true",
        help="(for testing) provides some output of the numbers",
    )
    parser.add_argument(
        "-v_every",
        metavar="X",
        dest="simup_verbose_every",
        type=int,
        help="verbose every X iterations (-verbose is required)",
        default=int(1e3),
    )
    parser.add_argument(
        "-rseed",
        metavar="X",
        type=int,
        help="(for debugging) specifies the seed for random number generator",
    )
    #NEW PARAM
    parser.add_argument(
        "-si",
        metavar="X",
        type=float,
        help="Initial SC density",
        default=0,
    )

    args = parser.parse_args()

    if args.rseed:
        np.random.seed(seed=args.rseed)
    # fixing random seed

    # OUTPUT FOLDER
    if args.f and os.path.isdir(args.output_folder):
        shutil.rmtree(args.output_folder, ignore_errors=True)
    try:
        os.mkdir(args.output_folder)
    except OSError as err:
        print("\n" + str(err))
        sys.exit("(use -f to overwrite it)\n")

    return args


class SimuParam:
    """Simulation parameters"""

    def __init__(self, args):

        self.verbose = args.verbose
        self.verbose_every = args.simup_verbose_every
        self.promoter_to_follow = args.promfollow

        self.Niterations = args.simup_Niterations
        # total number of iterations
        self.Ntranscripts_max = args.simup_Ntranscripts_max
        # maximum number of transcripts
        self.Nevery_transcripts = args.simup_Nevery_transcripts
        # writing out every Nevery_transcripts

        self.fo_out = args.output_folder  # output folder


class ModelParam:
    """Modelling parameters"""

    def __init__(self, args):

        self.dna = self._DNA(args)
        self.gene = self._Gene(args, self.dna)
        self.promoter = self._Promoter(args)
        self.rnap = self._RNAP(args, self.dna)
        self.topoI = self._TopoI(args)
        self.gyrase = self._Gyrase(args)
        self.coarse_g = self._CoarseGraining(
            args, self.dna, self.gene, self.rnap, self.topoI, self.gyrase
        )

    class _DNA:
        """DNA microscopic parameters"""

        def __init__(self, args):

            self.n = 10.5  # bp
            # structural properties

            self.A = args.dna_A  # pN.nm
            # torque <-> sigma properties

    class _Gene:
        """Gene parameters, exlcuding promoter features"""

        def __init__(self, args, dna):

            self.tss = args.gene_tss
            # TSS position
            self.rnap_xi = self.tss - int(args.rnap_excluded_length / 2)
            # initial position of rnap at the promoter
            self.Lk0_rnap_xi = self.rnap_xi / dna.n
            # TSS distance to upstream barrier in units of relaxed linking humber

            self.L = args.gene_L
            # gene length (in bp)
            self.Ldown = args.gene_Ldown
            # downstream distance to barrier (in bp)

            self.term = self.tss + self.L
            # location of termination 
            self.L_domain = self.tss + self.L + self.Ldown
            # total domain length
            self.Lk0_domain = self.L_domain / dna.n
            # in unit of linking number in the relaxed state
            self.Lk_domain = self.Lk0_domain * (args.si + 1)
            #MODIFIED: multiplied by arg.si + 1 term to set proper sigma
            # current state

    class _Promoter:
        """Promoter features"""

        def __init__(self, args):

            # BINDING RATE
            self.kb = args.promoter_kb
            self.kb_s = 1 / self.kb
            # corresponding scale (s)

            # FORMATION OF THE OPEN COMPLEX
            self.ko = args.promoter_ko
            self.ko_s = 1 / self.ko
            # corresponding scale (s)
            self.sigma_o = args.promoter_sigma_o
            # supercoiling density threshold
            self.o_width = args.promoter_sigma_width
            
            # OC ESCAPE RATE
            self.ke = args.promoter_ke
            self.ke_s = 1 / self.ke
            # corresponding scale (s)

    class _RNAP:
        """RNAP microscopic parameters"""

        def __init__(self, args, dna):

            self.excluded_length_elongation = args.rnap_excluded_length
            # in bp

            self.vm = args.rnap_vm
            # max velocity of the RNAP (bp per s)
            self.torque_stall = args.rnap_torque_stall
            # 18.5 pN in the presence of GreB (in vivo situation) (pN.nm)

            self.sigma_stall = -self.torque_stall / dna.A
            # (sigma's are signed here)

    class _TopoI:
        """TopoI activity"""

        def __init__(self, args):

            self.sigma_active = args.topoI_sigma_active
            # sigma below which TopoI is active

            self.lambda_ns = args.topoI_lambda_ns
            # non-specific (per base pair and per second)
            
            self.topoI_width = args.topoI_width

    class _Gyrase:
        """Gyrase activity
        rem: it is active whenever sigma is larger than stalling sigma's (i.e., stalling torque of the RNAP)
        """

        def __init__(self, args):

            self.lambda_ns = args.gyrase_lambda_ns
            # non-specific (per base pair and per second)
            
            self.sigma_active = args.gyrase_sigma_active
            # sigma above which gyrase is active, modified on 2024-01-11 in v1.2
            
            self.gyrase_width = args.gyrase_width
            # width of logistic curve for gyrase binding
            
            self.gyrase_sigma_procession = args.gyrase_sigma_procession
            # sigma above which gyrase can act
            
            self.gyrase_expected_actions = args.gyrase_expected_actions
            # expected number of gyrase actions per binding event
            
            if (self.sigma_active == None):
                self.sigma_active = -args.rnap_torque_stall / args.dna_A
            # if no sigma is specified, set the active gyrase level to the stalling torque (calculated same as in RNAP class)
            
            if (self.gyrase_sigma_procession == None):
                self.gyrase_sigma_procession = -args.rnap_torque_stall / args.dna_A
            # same as above for gyrase processivity
            
            

    class _CoarseGraining:
        """Coarse-graining parameters (modulated by resolution, i.e. dx)"""

        def __init__(self, args, dna, gene, rnap, topoI, gyrase):

            self.dx = args.coarse_g_dx
            # translocation step (in bp)
            self.dLk = self.dx / dna.n
            # translocation step in Lk unit
            self.tau_0 = self.dx / rnap.vm
            # corresponding time unit

            # TOPOI
            self.p_topoI_ns_per_bp = topoI.lambda_ns * self.tau_0
            # mean number of non-specific +1 (Lk) events per bp => Poisson process

            # GYRASE
            self.p_gyrase_ns_per_bp = gyrase.lambda_ns * self.tau_0 / 2
            # mean number of non-specific -2 (Lk) events per bp => Poisson process
            
            
            self.kmax_gyr = gyrase.lambda_ns * (gene.tss + gene.Ldown)
            # idealized binding rate for gyrase
            
            self.k_top = topoI.lambda_ns * (gene.tss + gene.Ldown)

    def _test(self):
        """Some tests"""

        # RNAP
        if self.rnap.sigma_stall > 0:
            sys.exit("stalling sigma's must be none-positive")

        # TOPOI
        if self.topoI.sigma_active <= self.rnap.sigma_stall:
            sys.exit("sigma_active <= sigma_stall is not biologically relevant")


class RNAP:
    def __init__(self):

        self.tb = 0
        # time at binding event
        self.tocf = 0
        # time at OC formation event
        self.tesc = 0
        # time at OC escape

        self.t_elongating = False
        # = True if RNAP is elongating

        self.X = 0
        # position in base pairs

        # UPSTREAM and DOWNSTREAM TOPOLOGICAL QUANTITIES
        self.Lk0 = {'up':None, 'down':None}
        # contour length in units of Lk (=relaxed linking number)
        self.Lk = {'up':None, 'down':None}
        # linking number
        self.sigma = {'up':None, 'down':None}
        # corresponding supercoiling density
        

class Trajectory:
    def __init__(self):

        self.niter = 0
        # number of iteration
        self.time = 0
        # real time

        self.Ntranscripts = 0
        # number of produced transcripts

        # STATISTICAL ANALYSIS OF TIMES
        self.prod_times = {"n": 0, "mean": 0}
        # time between the production of 2 RNA transcripts
        self.time_last_prod = 0
        self.b2b_times = {"n": 0, "mean": 0}
        # time between two binding events
        self.ref_binding_time = 0
        self.ocf_times = {"n": 0, "mean": 0}
        # time for OC to form once bound
        self.esc_times = {"n": 0, "mean": 0}
        # time to escape promoter once OC is formed
        self.initiation_times = {"n": 0, "mean": 0}
        # total time for initiation once it is bound
        self.elongation_times = {"n": 0, "mean": 0}
        # time for elongation once it is escaped from promoter

        self.sigma_at_ocf = []
        # supercoiling density at of OC formation stage
