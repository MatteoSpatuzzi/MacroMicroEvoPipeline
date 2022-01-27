#! \usr\bin\env python3

# File: setup_paml.py

# Sets up a standard PAML control file for sister pairs analysis and runs in console.

# Usage: 

# baseml : ./setup_paml.py pair1_noncoding.fasta baseml pair1.tre GTR
# codeml : ./setup_paml.py pair1_coding.fasta codeml pair1.tre --geneticCode standard

#Notes:

# baseml and codeml should be present on the command search path.

# For codeml, the alignment must be in frame and should not contain stop codons.
# If the process pauses, it is likely prompting the user what to do about stop codons - 
# best to remove these before analysis.

# The labels on the tips of the tree file should exactly match the sequence labels 
# in the alignment.
# For codeml, the tree file should be labelled with clade labels '$1', '$2' etc. to 
# specify which clades are subject to different omega (dN/dS) values.
# Alternatively, the "model" parameter could be set to 0 for non-varying rates.

### -------------------------DEPENDENCIES-----------------------###
import sys
import os
import re
import csv
import subprocess

from argparse import ArgumentParser

### ---------------------------GLOBAL-------------------------###

MODELS = "submodels.csv"

STD_GAMMA_CATS = 8

### -------------------------ARGUMENTS------------------------###

def make_parser():
    parser = ArgumentParser(description="Make a baseml or codeml control file given input alignments and trees.")
    parser.add_argument("alignment",
                        type = str,
                        help = "The alignment file to be analysed.")
    parser.add_argument("program",
                        type = str,
                        help = "The type of sequence - must be 'codeml' or 'baseml'",
                        choices = ["codeml", "baseml"])
    parser.add_argument("treefile",
                        type = str,
                        help = "The input tree file if required.")
    parser.add_argument("--model",
                        type = str,
                        help = "Only used for baseml runs. The model of DNA evolution.",
                        choices = ["JC", "F81", "K80", "HKY", "TrNef", 
                        "TrN", "TPM1", "TPM1uf", "TPM2", "TPM2uf", 
                        "TPM3", "TPM3uf", "TIM1ef", "TIM1", "TIM2ef", 
                        "TIM2", "TIM3ef", "TIM3", "TVMef", "TVM", "SYM", 
                        "GTR"],
                        default="GTR")
    parser.add_argument("--geneticCode",
                        type = str,
                        help = "Only required for codeml runs. Must be 'standard' or 'vert_mito'",
                        choices = ["standard", "vert_mito"])
    parser.add_argument("--doNotRun",
                        help = "Including this option just writes the control file without running PAML.",
                        action='store_true')
    return parser


### ------------------------MAIN PROCEDURE---------------------###

if __name__ == "__main__":

    parser = make_parser()
    args = parser.parse_args()
    
    if args.program == "codeml" and args.geneticCode is None:
        print ("Please specify the genetic code using --geneticCode=standard or --genetic_code=vert_mito")
    if args.program == "baseml" and args.geneticCode is not None:
        print("Warning: this is a baseml run. --geneticCode argument will be ignored.")
    
    
    ctl_outname = os.path.split(os.path.splitext(args.alignment)[0])[-1] + "_%s.ctl" %(args.program)

    modstr = ""
    if args.program == "baseml":
    
        model_dict = {}
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
        MODELS), 
        'r') as f:
            modreader = csv.DictReader(f)
            model_dict = {row["Model"] : row for row in modreader}
    
        modint = 9 if model_dict[args.model]["ueFreqs"] == 'u' else 10
        modstr = "%d [%s %s]" % (modint, 
                                    model_dict[args.model]["nRates"],
                                    model_dict[args.model]["rateGroups"])
    
    else:
        modstr = "7" # We just use GTR for codeml

    ctl_dict = {    "seqfile"     :    args.alignment, # The phylip or FASTA file with aligned sequence data
                    "treefile"    :    args.treefile,  # The reference tree/subtree. Codeml requires labels to be inserted to determine where to apply different nonsynonymous/synonymous substitution rates.
                    "outfile"     :    os.path.splitext(args.alignment)[0] + "_" + args.program + ".txt", # Location of PAML output

                    "noisy"       :    9, # Maximum output
                    "ndata"       :    1, # Number of data sets - almost always 1
                    "RateAncestor":    0, # Set to 1 only if doing ancestral sequence reconstruction
                    "getSE"       :    1, # Calculate standard errors for parameter values
                    
                    "runmode"     :    0, # Use the reference tree topology, do not perform phylogenetic reconstruction
                    
                    "clock"       :    0, # Estimate genetiv divergence in substitutions/site, do not estimate relative divergence times
                    
                    "fix_blength" :    0, # Ignore branch lengths in reference tree
                    
                    "Mgene"       :    0, # For alignments where different models are being fit to different data partitions. 0 = unpartitioned data

                    
                    "fix_alpha"   :    1, # Do not estimate the distribution of rates among DNA sites. Not needed for simple analyses.
                    "Malpha"      :    0, # Only doing one partition, so do not need different alpha shape parameters for different genes
                    "ncatG"       :    STD_GAMMA_CATS, # Number of categories into which evolutionary rates at sites fall - Usually not important.
                    
                    "fix_kappa"   :    0, # Estimate the ratio of transitions to transversions
                    "kappa"       :    2, # Initial value for transition/transversion ratio estimate

                    "cleandata"   :    0, # No automatic removal of gaps
                    
                    "method"      :    0, # Optimisation method for branch lengths. Never need to change.
                    "Small_Diff"  :    0.5e-8 # Threshold for comparing floating point values. Can lower if numerical instabilty occurs.
                    
                    
    }

    baseml_dict = { "verbose"     :    1, # show detailed output
                    "model"       :    modstr, # Specify the model of DNA substitution. See PAML manual for the format.
                    "nhomo"       :    0, # The model of DNA substitution does not change throughout the tree.
                    "alpha"       :    0.1, # Gives a diffuse distribution of rates among sites allowing considerable variation.
    
    }

    codeml_dict = { "model"       :    2, # This tells PAML to fit different nonsynonymous/synonymous substitution rates to different branches of the tree, according to labels ($1, $2 etc.) inserted into the tree file. Set to 0 for equal rates across the tree.
                    "seqtype"     :    1, # Specifies that the alignment is a codon sequence. Seqtype=2 is for amino acids ("aaml").
                    "CodonFreq"   :    2, # Codon frequency parameters in the codon model are calculated from averages per position - good balance of accuracy & not too many parameters
                    "icode"       :    0 if args.geneticCode == 'standard' else 1, # Integer for the genetic code used to link base and codon evolution. Other codes are available - see PAML manual.
                    "NSsites"     :    0, # This is used to search for specific sites under selection - we are comparing genomic rates so no need for this
                    "aaDist"      :    0, # Assume that all amino acid substitutions happen at equal rates. Fine for simple analyses.

                    "fix_omega"   :    0, # Omega is the parameter representing the ratio of nonsynonymous to synonymous substitution rate - that is, the strength of selection. fix_omega = 0 tells PAML to estimate this value.
                    "omega"       :    0.4, # Initial estimate of omega.


                    "alpha"       :    0, # Gamma distribution of rates among sites is not available for codeml, so set this to 0
                 
    }
    
    if args.program == "baseml":
        ctl_dict.update(baseml_dict)
    elif args.program == "codeml":
        ctl_dict.update(codeml_dict)
    
    with open(ctl_outname, 'w') as f:
    
        for k in ctl_dict.keys():

            f.writelines('\t%s\t=\t%s' % (k, str(ctl_dict[k]) + '\n'))
    
    print(args.program + " control file written to " + ctl_outname)

    if not args.doNotRun and os.path.exists(ctl_outname):
        completed = subprocess.run([args.program, ctl_outname])
        sys.exit(completed.returncode)
    

### ---------------------------END OF FILE--------------------###