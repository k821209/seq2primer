from __future__ import print_function
import sys
sys.path.append('/ref/analysis/pipelines/')
import primer3
import kang

file_seq = sys.argv[1] # fasta file 
Outfile  = file_seq+'.primers.txt'


dicFa = kang.Fasta2dic(file_seq)

def get_opt_cloningprimer_pair(seq):
    oligo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=20, dv_conc=1.5, dntp_conc=0.8, dna_conc=50, max_nn_length=60)
    # Change if you have personal condition for PCR 
    ## mv_conc   : The millimolar (mM) concentration of monovalent salt cations (usually KCl) in the PCR.
    ## dv_conc   : The millimolar concentration of divalent salt cations (usually MgCl^(2+)) in the PCR.
    ## dntp_conc : The millimolar concentration of the sum of all deoxyribonucleotide triphosphates.
    ## dna_conc  : A value to use as nanomolar (nM) concentration of each annealing oligo over the course the PCR. 
    ## max_nn_length : longest seq length for primer Tm analysis
    primerF_list = []
    primerR_list = []
    for i in range(16,22):
        primerF = seq[0:i]
        primerR = kang.rev_comp(seq[-i-3:-3]) # stop codon excluded
        if 56 <= oligo_calc.calcTm(primerF) <= 63:
            primerF_list.append(primerF)
        else:
            pass
            #print (primerF, oligo_calc.calcTm(primerF) )
        if 56 <= oligo_calc.calcTm(primerR) <= 63:
            primerR_list.append(primerR)

    primerF_list_Tm = [oligo_calc.calcTm(x) for x in primerF_list]
    primerR_list_Tm = [oligo_calc.calcTm(x) for x in primerR_list]


    comp_tm = []
    for i,tmF in enumerate(primerF_list_Tm):
        for j,tmR in enumerate(primerR_list_Tm):
            if abs(tmF-tmR) <= 3:  
                comp_tm.append([i,j,tmF,tmR])
    comp_tm.sort(key=lambda x : x[2]+x[3])

    try:
        ixF,ixR = comp_tm[0][0],comp_tm[0][1]
        return primerF_list[ixF],primerF_list_Tm[ixF],primerR_list[ixR],primerR_list_Tm[ixR]
    except IndexError:
        return None,None,None,None

with open(Outfile,'w') as f:
    for key in dicFa:
        print(key,'\t'.join(map(str,get_opt_cloningprimer_pair(dicFa[key]))),sep='\t',file=f)
