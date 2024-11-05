#Tutorial: https://www.biotite-python.org/latest/examples/gallery/sequence/homology/avidin_alignment.html#sphx-glr-examples-gallery-sequence-homology-avidin-alignment-py
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics
import biotite.sequence.io.fasta as fasta
import sys
import matplotlib.pyplot as plt

##Comparing Papain and RIM13 cystiene peptidase from Cryptococcus neoformans

#Setting paths to the two fasta files that will be compared
fasta_1 = sys.argv[1]
fasta_2 = sys.argv[2]

#Parsing throught the two fasta files
rim13_fasta = fasta.FastaFile.read(fasta_1)
papain_fasta= fasta.FastaFile.read(fasta_2)

#Creating 'ProteinSequence' classes from the sequences in the two fasta files
for name, sequence in rim13_fasta.items():
    rim13_seq = seq.ProteinSequence(sequence)

for name, sequence in papain_fasta.items():
    papain_seq = seq.ProteinSequence(sequence)

#Initiate matrix for the alignment of the two protein sequences
matrix = align.SubstitutionMatrix.std_protein_matrix()

#Peform sequence alignmnet with a affine gap penalty of -10,-1
alignments = align.align_optimal(rim13_seq,papain_seq, matrix, gap_penalty= (-10,-1), terminal_penalty=False)

#plot the sequence alignment
fig = plt.figure(figsize=(10.0, 4.5))
ax = plt.subplot(111)
graphics.plot_alignment_similarity_based(axes=ax, 
                                         alignment=alignments[0],
                                         matrix=matrix, 
                                         labels=["Rim13", "Papain"], 
                                         show_numbers=True, 
                                         show_line_position=True)
fig.tight_layout
plt.show()