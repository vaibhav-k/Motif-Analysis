import pandas as pd
import os
import csv
import re

def analyse_variable_length_sequences(sequences, biggest):
	analysis = []
	total = len(sequences)
	for i in range(biggest-1):
		agct = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
		total_counted = 0
		for j in range(total):
			if len(sequences[j]) <= i:
				continue
			letter = sequences[j][i]
			agct[letter] += 1
			total_counted += 1
		for key, val in agct.items():
			agct[key] *= 100
			agct[key] /= total_counted
		analysis.append(agct)
	return analysis

if __name__ == "__main__":
	motif1 = "AGCT"
	motif2 = "AAAG"
	path = "Disease - All.txt"
	all_sequences = []
	gene_count = {}
	with open(path, 'r') as file:
		for line in file:
			count = []
			for i in range(31):
				count.append(0)
			gene_id = line.split('\t')[0]
			gene_sequence = line.split('\t')[1]
			sequences = re.findall(r"(?<=AGCT).*?(?=AAAG)", gene_sequence)
			for seq in sequences:
				if len(seq) > -1 and len(seq) < 31:
					count[len(seq)] += 1
				all_sequences.append(seq)
			gene_count[gene_id] = count
	frequency = []
	count = []
	for dist in range(1,31):
		spacer_count = 0
		spacer_frequency = []
		for seq in all_sequences:
			if len(seq) == dist:
				spacer_frequency.append(seq)
				spacer_count += 1
		frequency.append(spacer_frequency)
		count.append(spacer_count)
	overall_analysis = analyse_variable_length_sequences(all_sequences, 31)
	with open("AGCT-AAAG "+path.split('.txt')[0]+", flank.csv", 'w') as file:
		pos=1
		file.write("Position,A,G,C,T\n")
		for item in overall_analysis:
			file.write("%d," % (pos))
			pos += 1
			for key, val in item.items():
				file.write("%s: %f, " % (key, val))
			file.write('\n')
	'''
	print(gene_count)
	with open("AAAG-AGCT gene spacer count.csv", 'w') as file:
		for key, val in gene_count.items():
			file.write(key+",")
			for v in val:
				file.write(str(v) + ',')
			file.write('\n')
	#print(frequency)
	df = pd.DataFrame({
		"Spacer Frequency": count
		})
	df.to_csv("AGCT-AAAG spacer frequency Disease - All.csv", index=False)
	'''