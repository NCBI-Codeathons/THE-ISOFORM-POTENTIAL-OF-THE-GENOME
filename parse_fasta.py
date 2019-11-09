from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


# All ORF start with Methionine
record_dict = SeqIO.index("files/lncipedia_4_0_hg19_orf_max.fasta", "fasta")
test_max =[]
records_all = {}
record_ids_max = {}
for record_key in record_dict.keys():
    length_orf =  len(record_dict[record_key].seq)
    test_max.append(length_orf)
    record_id = record_dict[record_key].id.split('_')[0]
    try:
        records_all[record_id].append(length_orf)
    except KeyError:
        records_all[record_id] = [length_orf]

for k,v in records_all.items():
    record_ids_max[k] = max(v)

plt.hist(np.log2(test_max), 25, normed=0, facecolor='#C6D7E9',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.xlim([2,14])
plt.axvline(np.log2(np.median(test_max)), color='k', linestyle='dashed', linewidth=2)
#plt.text(np.log2(np.median(test)+1), 50000, 'median=24 AA', fontsize=8, fontweight='bold')
plt.suptitle('Range of Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('log2(Peptide length)',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/range_of_universe_peptide_lengths.png')
plt.clf()

plt.hist(np.log2(record_ids_max.values()), 25, normed=0, facecolor='#C6D7E9',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.xlim([2,14])
plt.axvline(np.log2(np.median(record_ids_max.values())), color='k', linestyle='dashed', linewidth=2)
#plt.text(np.log2(np.median(record_ids.values())+1), 5700, 'median=60 AA', fontsize=8, fontweight='bold')
plt.suptitle('Longest Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('log2(Peptide length)',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/longest_of_universe_peptide_lengths.png')
plt.clf()



# All ORF start with Methionine < 25kDa
record_dict = SeqIO.index("files/lncipedia_4_0_hg19_orf.fasta", "fasta")
test =[]
records = {}
record_ids = {}
for record_key in record_dict.keys():
    length_orf =  len(record_dict[record_key].seq)
    test.append(length_orf)
    record_id = record_dict[record_key].id.split('_')[0]
    try:
        records[record_id].append(length_orf)
    except KeyError:
        records[record_id] = [length_orf]

for k,v in records.items():
    record_ids[k] = max(v)

plt.hist(test, 50, normed=0, facecolor='#E0A8A6',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.xlim([0,250])
plt.axvline(np.median(test), color='k', linestyle='dashed', linewidth=2)
#plt.text(25, 120000, 'median=24 AA', fontsize=8, fontweight='bold')
plt.suptitle('Range of Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('Peptide length',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/range_of_peptide_lengths.png')
plt.clf()

plt.hist(record_ids.values(), 45, normed=0, facecolor='#E0A8A6',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.xlim([0,250])
plt.ylim([0,8000])
plt.axvline(np.median(record_ids.values()), color='k', linestyle='dashed', linewidth=2)
#plt.text(60, 1600, 'median=60 AA', fontsize=8, fontweight='bold')
plt.suptitle('Longest Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('Peptide length',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/longest_of_peptide_lengths.png')
plt.clf()


# All Protein Coding genes
record_dict = SeqIO.index("files/gencode.v19.annotation_protein_coding_orf.fasta", "fasta")
test_max =[]
records_all = {}
record_ids_max = {}
for record_key in record_dict.keys():
    length_orf =  len(record_dict[record_key].seq)
    test_max.append(length_orf)
    record_id = record_dict[record_key].id.split('_')[0]
    try:
        records_all[record_id].append(length_orf)
    except KeyError:
        records_all[record_id] = [length_orf]

for k,v in records_all.items():
    record_ids_max[k] = max(v)

plt.hist(test_max, 50, normed=0, facecolor='#DC8E60',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.ylim([0,350000])
plt.xlim([0,250])
plt.axvline(np.median(test_max), color='k', linestyle='dashed', linewidth=2)
#plt.text(26, 120000, 'median=25 AA', fontsize=8, fontweight='bold')
plt.suptitle('Range of Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('Peptide length',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/range_of_peptide_lengths_pcc.png')
plt.clf()

plt.hist(record_ids_max.values(), 45, normed=0, facecolor='#DC8E60',linewidth=3,edgecolor='black')
plt.tick_params(axis='y', direction='out')
plt.tick_params(axis='x', direction='out')
plt.tick_params(top='off', right='off')
plt.xlim([0,250])
plt.ylim([0,7000])
plt.axvline(np.median(record_ids_max.values()), color='k', linestyle='dashed', linewidth=2)
#plt.text(60, 1600, 'median=96 AA', fontsize=8, fontweight='bold')
plt.suptitle('Longest Peptide lengths ORF', fontsize=8, fontweight='bold')
plt.xlabel('Peptide length',fontsize=8, fontweight='bold')
plt.ylabel('Frequency',fontsize=8, fontweight='bold')
plt.savefig('figures/longest_of_peptide_lengths_pcc.png')
plt.clf()
