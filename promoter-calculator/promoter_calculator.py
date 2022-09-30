'''version 1.1'''

import random, sys, pickle, collections, operator, itertools, time, math, os
from .util import *
from collections import defaultdict
import numpy as np
import math
from copy import copy
import argparse
import os
from dataclasses import dataclass

# k and BETA for La Fleur dataset
LOGK   = -2.80271176
BETA    = 0.81632623

@dataclass
class PromoCalcResults:
    """Class to hold results of the promoter prediction"""
    promoter_sequence: str
    TSS: int
    UP: str
    hex35: str
    spacer: str
    hex10: str
    disc: str
    ITR: str
    dG_total: float
    dG_10: float
    dG_35: float
    dG_disc: float
    dG_ITR: float
    dG_ext10: float
    dG_spacer: float
    dG_bind: float
    dG_UP: float
    Tx_rate: float
    UP_position: int
    hex35_position: int
    spacer_position: int
    hex10_position: int
    disc_position: int
    fraction: int
    strand: str
    drop: bool
    length: int

    def __getitem__(self, key):
        if key in self.__annotations__:
            return getattr(self, key)
        else:
            raise KeyError(key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def keys(self):
        return self.__annotations__.keys()

def unpickler(infile):
    """ Unpickle a file """
    with open(infile, 'rb') as handle:
        obj= pickle.load(handle)
    return obj

def _revcomp(seq):
    """ Reverse complement a sequence """
    revcomp = {'U' : 'A', 'A' : 'T', 'G' : 'C', 'T' : 'A', 'C' : 'G'}
    return "".join([revcomp[letter] for letter in seq[::-1] ])

def get_matrices(two_mer_encoder, three_mer_encoder, spacer_encoder, coeffs):
    """ Get matrices for linear model """

    #Extract dG values from model coefficients
    ref10_0 = coeffs.tolist()[0:64]
    ref10_3 = coeffs.tolist()[64:128]
    ref35_0 = coeffs.tolist()[128:192]
    ref35_3 = coeffs.tolist()[192:256]
    discs   = coeffs.tolist()[256:256+64]
    x10     = coeffs.tolist()[256+64:256+64+16]
    spacs   = coeffs.tolist()[256+64+16:256+64+16+3]

    # make dG matrices for each feature
    dg10_0  = get_dg_matrices(ref10_0,three_mer_encoder)
    dg10_3  = get_dg_matrices(ref10_3,three_mer_encoder)
    dg35_0  = get_dg_matrices(ref35_0,three_mer_encoder)
    dg35_3  = get_dg_matrices(ref35_3,three_mer_encoder)
    dmers   = get_dg_matrices(discs,three_mer_encoder)
    x10mers = get_dg_matrices(x10,two_mer_encoder)
    spacers = get_dg_matrices(spacs, spacer_encoder)

    return dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers

def scan_arbitrary(inputSequence, model, inters, constraints, dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers):
    """Scan sequence left to right with no TSS information. Calc dG of all possible promoter configurations."""
    seq_query = {}
    upstream = constraints[0]
    downstream = constraints[1]
    sequence = upstream + inputSequence + downstream

    # first 20 nt will be initial UP candidate
    for i in range(0,len(sequence)):
        tempUP = sequence[i:i+24]
        temp35 = sequence[i+25:25+i+6]  # leaves 1 nt between h35 and UPs

        # bounds defined by what was present during parameterization
        for j in range(15,21):
            tempspacer = sequence[i+25+6:25+i+6+j]
            temp10     = sequence[25+i+j+6:25+i+j+12]
            for k in range(6,11):
                tempdisc  = sequence[25+i+j+12:25+i+j+12+k]
                tempITR   =sequence[25+i+j+12+k:45+i+j+12+k]
                if len(tempITR) < 20:
                    continue
                else:
                    model_results= linear_free_energy_model(tempUP,
                                                            temp35,
                                                            tempspacer,
                                                            temp10,
                                                            tempdisc,
                                                            tempITR,
                                                            dg10_0,
                                                            dg10_3,
                                                            dg35_0,
                                                            dg35_3,
                                                            dmers,
                                                            x10mers,
                                                            spacers,
                                                            model,
                                                            inters)
                    TSS_distance = i + len(tempUP) + len(temp35) + len(tempspacer) + len(temp10) + len(tempdisc)
                    # seq_query[(float(dG_bind), float(dG_total), TSS_distance)] = ((tempUP, temp35, tempspacer, temp10, tempdisc, tempITR),(dg10, dg35, dg_disc, dg_ITR, dg_ext10, dg_spacer, dg_UP))
                    seq_query[(float(model_results.dg_total), float(model_results.dg_apparent), TSS_distance)] = ((tempUP, temp35, tempspacer, temp10, tempdisc, tempITR),(model_results.dg_10, model_results.dg_35, model_results.dg_disc, model_results.dg_ITR, model_results.dg_ext10, model_results.dg_spacer, model_results.dg_UP))

    print("best: ", min(seq_query.items(), key=operator.itemgetter(0)))

    best = (collections.OrderedDict(sorted(seq_query.items())), min(seq_query.items(), key=operator.itemgetter(0)))
    return best, seq_query

def linear_free_energy_model(UP, h35, spacer, h10, disc, ITR, dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers, coeffs, inters):
    """ Calculate dG of a promoter configuration """

    prox_UP = UP[-int(len(UP)/2)::]
    dist_UP = UP[0:int(len(UP)/2)]

    # CATEGORICAL FEATURES
    ext10           = spacer[-3:-1] # TGN motif, contacts sigma
    hex10_0         = h10[0:3]
    hex10_3         = h10[3::]
    hex35_0         = h35[0:3]
    hex35_3         = h35[3::]
    disc_first_3mer = disc[0:3]
    spacer_length   = str(len(spacer))

    # NUMERICAL FEATURES
    _, _, dg_ITR     = calc_DNA_RNA_hybrid_energy(ITR) # calc R-loop strength
    rigidity                 = calc_rigidity(seq = UP + h35 + spacer[0:14])

    width_proxy_prox = calc_groove_width(prox_UP)
    width_proxy_dist = calc_groove_width(dist_UP)

    # NORMALIZE NUMERICAL FEATURES BY MAX IN TRAINING SET
    numericals         = np.array([width_proxy_dist, width_proxy_prox, dg_ITR, rigidity])
    normalizing_values = [256.0, 255.0, 4.300000000000002, 25.780434782608694]
    numerical_coefs    = np.array(coeffs.tolist()[-4::])
    normald            = np.divide(numericals,normalizing_values)
    dg_numerical       = np.multiply(normald, numerical_coefs)
    dg10      = dg10_0[hex10_0] + dg10_3[hex10_3]
    dg35      = dg35_0[hex35_0] + dg35_3[hex35_3]
    dg_disc   = dmers[disc_first_3mer]
    dg_ITR    = dg_numerical[-2]
    dg_ext10  = x10mers[ext10]

    x = float(spacer_length)
    dg_spacer = 0.1463*x**2 - 4.9113*x + 41.119

    dg_UP        = dg_numerical[0] + dg_numerical[1] + dg_numerical[-1]
    dG_apparent  = (dg10 + dg35 + dg_disc + dg_ITR + dg_ext10 + dg_spacer + dg_UP + inters[0] - LOGK)/BETA
    dG_total     = dg10 + dg35 + dg_disc + dg_ITR + dg_ext10 + dg_spacer + dg_UP + inters[0]

    @dataclass
    class PromoModelResults:
        """Class to hold results of the linear free energy model"""
        dg_total: float
        dg_apparent: float
        dg_10: float
        dg_35: float
        dg_disc: float
        dg_ITR: float
        dg_ext10: float
        dg_spacer: float
        dg_UP: float

    return PromoModelResults(dG_total,  # @TODO: Univy this with PromoPredictionResults
                             dG_apparent,
                             dg10,
                             dg35,
                             dg_disc,
                             dg_ITR,
                             dg_ext10,
                             dg_spacer,
                             dg_UP)

def predict(sequence, constraints):
    """ Predict the free energy of a promoter configuration """

    # Initialize model and matrices
    install_location = os.path.dirname(os.path.realpath(__file__))
    layer1 = np.load(install_location + '/free_energy_coeffs.npy')
    inters = np.load(install_location + '/model_intercept.npy')

    two_mer_encoder   = kmer_encoders(k = 2)
    three_mer_encoder = kmer_encoders(k = 3)
    spacer_encoder    = length_encoders(16, 18)
    dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers = get_matrices(two_mer_encoder = two_mer_encoder, three_mer_encoder = three_mer_encoder, spacer_encoder = spacer_encoder, coeffs = layer1)

    # Scan DNA and return predictions
    (_, result), query   = scan_arbitrary(inputSequence = sequence,
                                            model = layer1, inters = inters, constraints = constraints, dg10_0 = dg10_0, dg10_3 = dg10_3,
                                            dg35_0 =dg35_0, dg35_3 = dg35_3, dmers = dmers , x10mers = x10mers, spacers = spacers)

    @dataclass
    class PromoPredictResults:
        """Class to hold results of the promoter prediction"""
        dg_total: float
        query: dict
        UP: str
        h35: str
        spacer: str
        h10: str
        disc: str
        ITR: str

    dG_total, UP, h35, spacer, h10, disc, ITR = result[0][0], result[1][0][0],result[1][0][1],result[1][0][2],result[1][0][3], result[1][0][4], result[1][0][5]
    return PromoPredictResults(dG_total, query, UP, h35, spacer, h10, disc, ITR)

class Promoter_Calculator(object):
    """ Class to calculate the free energy of a promoter configuration """

    def __init__(self, organism = 'Escherichia coli str. K-12 substr. MG1655',
                       sigmaLevels = {'70' : 1.0, '19' : 0.0, '24' : 0.0, '28' : 0.0, '32' : 0.0, '38' : 0.0, '54' : 0.0}):

        # Initialize model and matrices
        path = os.path.dirname(os.path.abspath(__file__))
        self.layer1 = np.load(path + '/free_energy_coeffs.npy')
        self.inters = np.load(path + '/model_intercept.npy')

        self.two_mer_encoder   = kmer_encoders(k = 2)
        self.three_mer_encoder = kmer_encoders(k = 3)
        self.spacer_encoder    = length_encoders(16, 18)
        self.dg10_0, self.dg10_3, self.dg35_0, self.dg35_3, self.dmers, self.x10mers, self.spacers = get_matrices(two_mer_encoder = self.two_mer_encoder, three_mer_encoder = self.three_mer_encoder, spacer_encoder = self.spacer_encoder, coeffs = self.layer1)

        self.model = self.layer1
        self.organism = organism
        self.sigmaLevels = sigmaLevels

        if organism == 'in vitro':
            self.K = 42.00000
            self.BETA = 0.81632623
        elif organism == 'Escherichia coli str. K-12 substr. MG1655':
            self.K = 42.00000
            self.BETA = 1.636217004872062
        else:
            self.K = 42.00000
            self.BETA = 1.636217004872062

    # Identify promoter with minimum dG_total (across many possible promoter states) for each TSS position in an inputted sequence.
    def predict(self, sequence, TSS_range):
        """ Predict the free energy of a promoter configuration """

        UPS_length = 24
        HEX35_length = 6
        UPS_HEX35_SPACER = 1
        SPACER_length_range = [15, 21]
        HEX10_length = 6
        DISC_length_range = [6, 11]
        ITR_length = 20
        All_States = {}
        Min_States = {}

        #Specify fixed TSS range
        for TSS in range(TSS_range[0], TSS_range[1]):
            All_States[TSS] = {}
            for DISC_length in range(DISC_length_range[0],DISC_length_range[1]):

                if TSS - DISC_length >= 0 and TSS + ITR_length <= len(sequence):
                    tempdisc = sequence[ TSS - DISC_length : TSS  ]
                    tempITR  = sequence[ TSS : TSS + ITR_length]

                    for SPACER_length in range(SPACER_length_range[0], SPACER_length_range[1]):

                        if TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER >= 0:
                            temp10     = sequence[ TSS - DISC_length - HEX10_length : TSS - DISC_length]
                            tempspacer = sequence[ TSS - DISC_length - HEX10_length - SPACER_length : TSS - DISC_length - HEX10_length ]
                            temp35     = sequence[ TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length : TSS - DISC_length - HEX10_length - SPACER_length]
                            tempUP     = sequence[ TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER:  TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER]

                            model_results = linear_free_energy_model(tempUP,
                                                                     temp35,
                                                                     tempspacer,
                                                                     temp10,
                                                                     tempdisc,
                                                                     tempITR,
                                                                     self.dg10_0,
                                                                     self.dg10_3,
                                                                     self.dg35_0,
                                                                     self.dg35_3,
                                                                     self.dmers,
                                                                     self.x10mers,
                                                                     self.spacers,
                                                                     self.model,
                                                                     self.inters)
                            dG_bind = (model_results.dg_10 +
                                       model_results.dg_35 +
                                       model_results.dg_spacer +
                                       model_results.dg_ext10 +
                                       model_results.dg_UP)

                            Tx_rate = self.K * math.exp(- self.BETA * model_results.dg_total )


                            result = {'promoter_sequence' : sequence[TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER : TSS + ITR_length ],
                                      'TSS' : TSS, 'UP' : tempUP, 'hex35' : temp35, 'spacer' : tempspacer, 'hex10' : temp10, 'disc' : tempdisc, 'ITR' : tempITR,
                                      'dG_total' : model_results.dg_total, 'dG_10' : model_results.dg_10, 'dG_35' : model_results.dg_35, 'dG_disc' : model_results.dg_disc, 'dG_ITR' : model_results.dg_ITR, 'dG_ext10' : model_results.dg_ext10, 'dG_spacer' : model_results.dg_spacer, 'dG_UP' : model_results.dg_UP, 'dG_bind' : dG_bind,
                                      'Tx_rate' : Tx_rate,
                                      'UP_position' : [TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER - UPS_length,
                                                       TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER],
                                      'hex35_position' : [TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length, TSS - DISC_length - HEX10_length - SPACER_length],
                                      'spacer_position' : [TSS - DISC_length - HEX10_length - SPACER_length, TSS - DISC_length - HEX10_length],
                                      'hex10_position' : [TSS - DISC_length - HEX10_length, TSS - DISC_length],
                                      'disc_position' : [TSS - DISC_length, TSS]
                                      }

                            result = PromoCalcResults(**result, fraction=None, strand=None, drop=False, length = False)

                            All_States[TSS][ (DISC_length, SPACER_length) ] = result
                            if TSS in Min_States:
                                if result['dG_total'] < Min_States[TSS]['dG_total']:  Min_States[TSS] = result
                            else:
                                Min_States[TSS] = result

        return (Min_States, All_States)

    def run(self, sequence, TSS_range = None):
        """ Run the predictor """

        if TSS_range is None: TSS_range = [0, len(sequence)]

        self.sequence = sequence
        self.TSS_range = TSS_range
        self.TSS_range_rev = [len(sequence) - TSS_range[1], len(sequence) - TSS_range[0]]

        # print "self.TSS_range_rev: ", self.TSS_range_rev

        fwd_sequence = sequence
        rev_sequence = _revcomp(sequence)
        (Forward_Min_States, Forward_All_States) = self.predict(fwd_sequence, TSS_range = self.TSS_range)
        (Reverse_Min_States_Temp, Reverse_All_States_Temp) = self.predict(rev_sequence, TSS_range = self.TSS_range_rev)

        Reverse_Min_States = {}
        Reverse_All_States = {}
        # 0  <------>500 fwd 500 bp
        # 500<------>0   rev 500 bp
        #      275 TSS
        #      200-300
        #500-275 = 225
        #
        for TSS in Reverse_Min_States_Temp.keys():
            Reverse_Min_States[len(sequence) - TSS] = Reverse_Min_States_Temp[TSS]
            Reverse_All_States[len(sequence) - TSS] = Reverse_All_States_Temp[TSS]

        self.Forward_Predictions_per_TSS = Forward_Min_States
        self.Reverse_Predictions_per_TSS = Reverse_Min_States

    def output(self):
        """ Output the results as a dictionary """
        output = {'organism' : self.organism,
                  'sigmaLevels' : self.sigmaLevels,
                  'K' : self.K,
                  'beta' : self.BETA,
                  'sequence' : self.sequence,
                  'TSS_range' : self.TSS_range,
                  'Forward_Predictions_per_TSS' : self.Forward_Predictions_per_TSS,
                  'Reverse_Predictions_per_TSS' : self.Reverse_Predictions_per_TSS
                }
        return output.copy()

