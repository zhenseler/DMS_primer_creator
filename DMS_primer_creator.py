"""This script creates primers for either PALS mutagenesis (Kitzman 2015) or One-Pot 
	Saturation mutagenesis (Wrenbeck 2016) from an input sequence. Can be used to either
	make programmed codon mutations, or 'NNN' degenerate primers."""


from sys import argv

try:
	x = argv[1]
except IndexError:
	print " \n Function: creates primers for One-pot Saturation Mutagenesis (Wrenbeck 2016)"+\
		'\n Usage: \n python DMS_primer creator.py argv1 - argv7 \n \n' +\
		' argv1 = input fasta filepath \n' +\
		' argv2 = mode (options: programmed or degenerate) \n' +\
		' argv3 = codon table filepath (codons must be in column one, no header) \n' +\
		' argv4 = number of first nucleotide to mutagenize \n' +\
		' argv5 = number of last nucleotide to mutagenize \n' +\
		' argv6 = sequence for 5\'-overhang \n' +\
        ' argv7 = sequence for 3\'-overhang of the primer (do not reverse complement) \n' +\
		' argv8 = output filepath, ending with .txt \n'
	quit()


input_fasta_file = argv[1]
mode = str(argv[2])
codon_table = argv[3]
mutagenesis_area_first_nt = argv[4]
mutagenesis_area_last_nt = argv[5]
fiveprime_overhang = argv[6]
threeprime_overhang = argv[7]
deletion_fiveprime_overhang = argv[8]
deletion_threeprime_overhang = argv[9]
output_folder = argv[10]

deletion_primer = ''
sequence = ''
tm_list = []
GC_list = []
primer_flank = 9
min_tm = 0
extra_nt = 0
tm = 0

for l in open(input_fasta_file,'U'):
	if not l.startswith('>'):
		sequence += l.strip('\n')

if (int(mutagenesis_area_last_nt) - int(mutagenesis_area_first_nt) + 1) % 3 != 0:
	raise ValueError('\n'+'\n'+'\t'+'Area of mutagenesis must be evenly divisible by 3.'+\
		'Check first/last nt'+'\n')
	
if int(mutagenesis_area_first_nt) < 30:
	raise ValueError('\n'+'\n'+'\t'+'Needs to be at least 30 nt flanking the 5-prime ' +\
		'end of mutagenesis region'+'\n')
if len(sequence) - int(mutagenesis_area_last_nt) < 30:
	raise ValueError('\n'+'\n'+'\t'+'Needs to be at least 30 nt flanking the 3-prime ' +\
	 	'end of mutagenesis region'+'\n')


out = open(output_folder,"w")
out.close()
out = open(output_folder,"a")


# The following is for programmed mutations, including codon deletions

if mode == 'programmed':

	out.write('\t'.join(('Primer','Length','Codon Mutated','Tm','% GC','\n')))

	codon_list = []
	for l in open(codon_table,'U'):
		codon_list.append(l.split()[0])
		
	for i, nt in enumerate(sequence):
		if i+1 >= int(mutagenesis_area_first_nt) and i+1 <= int(mutagenesis_area_last_nt):
			current_nt = i-int(mutagenesis_area_first_nt)+2
			if current_nt % 3 == 1:
				for opt_codon in codon_list:
					mismatch = 0
					codon = sequence[i:i+3]
					for nuc in [0,1,2]:
						if opt_codon[nuc] != codon[nuc]:
							mismatch += 1
				
			# From codon, starts extending primer on both ends until tm > 78oC

					while tm < 78:
						primer_flank += 1
						primer = sequence[i-primer_flank:i] + opt_codon +\
						 sequence[i+3:i+primer_flank+3]
						GC = primer.count('G')+primer.count('C')
						GC_percent = 100*(float(GC)/len(primer))
						tm = 81.5 + 0.41 * GC_percent - (675./len(primer))\
						 - 100*(float(mismatch)/len(primer))


					#If no GC clamp, looks 3' up to 5 nt for one, and extends primer 

					if primer.endswith('A') or primer.endswith('T'):
						if ('C' in sequence[i+3+primer_flank:i+3+primer_flank+5]) or\
						 ('G' in sequence[i+3+primer_flank:i+3+primer_flank+5]):
							for extra_nt in [0,1,2,3,4,5]:
								if (sequence[i+3+primer_flank+extra_nt] == 'C') or\
								 (sequence[i+3+primer_flank+extra_nt] == 'G'):
									primer = sequence[i-primer_flank:i] + opt_codon +\
									 sequence[i+3:i+3+primer_flank+extra_nt+1]
									tm = 81.5 + 0.41 * GC_percent - (675./len(primer))\
									- 100*(float(mismatch)/len(primer))
									GC = primer.count('G')+primer.count('C')
									GC_percent = 100*(float(GC)/len(primer))
									break

					output = '\t'.join((primer, str(len(primer)), codon, \
					str(round(tm,2)), str(round(GC_percent,2)), '\n')) 
					out.write(output)

					tm_list = []
					GC_list = []
					tm = 0
					primer_flank = 9
					extra_nt = 0

				# Following is to make the codon deletion primer
				while tm < 78:
					primer_flank += 1
					primer = sequence[i-primer_flank:i] + sequence[i+3:i+primer_flank+3]
					GC = primer.count('G')+primer.count('C')
					GC_percent = 100*(float(GC)/len(primer))
					tm = 81.5 + 0.41 * GC_percent - (675./len(primer)) -\
					 100*(float(mismatch)/len(primer))
	
				if primer.endswith('A') or primer.endswith('T'):
						if ('C' in sequence[i+3+primer_flank:i+3+primer_flank+5]) or\
						 ('G' in sequence[i+3+primer_flank:i+3+primer_flank+5]):
							for extra_nt in [0,1,2,3,4,5]:
								if (sequence[i+3+primer_flank+extra_nt] == 'C') or\
								 (sequence[i+3+primer_flank+extra_nt] == 'G'):
									primer = sequence[i-primer_flank:i] + opt_codon +\
									 sequence[i+3:i+3+primer_flank+extra_nt+1]
									tm = 81.5 + 0.41 * GC_percent - (675./len(primer))\
									- 100*(float(mismatch)/len(primer))
									GC = primer.count('G')+primer.count('C')
									GC_percent = 100*(float(GC)/len(primer))
									break

				output = '\t'.join((primer, str(len(primer)), codon, \
				str(round(tm,2)), str(round(GC_percent,2)), '\n'))
				out.write(output)
				
				
# The following is for programmed mutations using the most optimal codons (user provided)

if mode == 'degenerate':

	out.write('\t'.join(('Primer','Length','Codon Mutated','Avg Tm','Min Tm','Max Tm',\
	'Avg % GC','Min % GC','Max % GC','\n')))

	# Starting from mutated codon, extends primer outwards until minimum Tm is 76oC

	for i, nt in enumerate(sequence):
		if i+1 >= int(mutagenesis_area_first_nt) and i+1 <= int(mutagenesis_area_last_nt):
			current_nt = i-int(mutagenesis_area_first_nt)+2
			if current_nt % 3 == 1:
				while min_tm < 76:
					primer_flank += 1
					primer = sequence[i-primer_flank:i] + 'NNN' +\
					 sequence[i+3:i+primer_flank+3]
					min_GC = primer.count('G')+primer.count('C')
					min_GC_percent = 100*(float(min_GC)/len(primer))
					min_tm = 81.5 + 0.41 * min_GC_percent - \
					(675./len(primer)) - 100*(3./len(primer))
				primer_before_GC = primer
			

				# If no GC clamp, looks for a G or C5 nt down 3' end
				# and extends primer to there
			
				if primer.endswith('A') or primer.endswith('T'):
					if ('C' in sequence[i+3+primer_flank:i+3+primer_flank+5]) or\
					 ('G' in sequence[i+3+primer_flank:i+3+primer_flank+5]):
						for extra_nt in [0,1,2,3,4,5]:
							if (sequence[i+3+primer_flank+extra_nt] == 'C') or\
							 (sequence[i+3+primer_flank+extra_nt] == 'G'):
								primer = sequence[i-primer_flank:i] + 'NNN' +\
								 sequence[i+3:i+3+primer_flank+extra_nt+1]
								break

				max_GC = primer.count('G')+primer.count('C')+3
				max_GC_percent = 100*(float(max_GC)/len(primer))

				GC = primer.count('G')+primer.count('C')

				for GC_count_in_codon in [0,1,2,3]:
					GC_percent = 100*(float(GC + GC_count_in_codon)/len(primer))
					GC_list.append(GC_percent)
					for MM_count_in_codon in [0,1,2,3]:
						tm = 81.5 + 0.41 * min_GC_percent - (675./len(primer)) -\
						 100*(float(MM_count_in_codon)/len(primer))
						tm_list.append(tm)
				average_tm = sum(tm_list)/len(tm_list)
				average_GC = sum(GC_list)/len(GC_list)
			
				max_tm = 81.5 + 0.41 * max_GC_percent - (675./len(primer))
	
				codon = sequence[i:i+3]
				primer = str(fiveprime_overhang) + primer + str(threeprime_overhang)
			
				output = '\t'.join((primer, str(len(primer)), codon, \
				str(round(average_tm,2)), str(round(min_tm,2)), str(round(max_tm,2)), \
				str(round(average_GC,2)), str(round(min_GC_percent,2)), \
				str(round(max_GC_percent,2)), '\n')) 
				out.write('\n')
				out.write(output)
				deletion_primer = primer.replace('N','')
				deletion_primer = deletion_primer.replace(fiveprime_overhang,deletion_fiveprime_overhang)
				deletion_primer = deletion_primer.replace(threeprime_overhang,deletion_threeprime_overhang)
                out.write(deletion_primer)

                tm_list = []
                GC_list = []
                min_tm = 0
                primer_flank = 9
                extra_nt = 0
                deletion_primer = ''

out.close()