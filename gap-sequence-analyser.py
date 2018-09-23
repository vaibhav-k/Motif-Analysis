import pandas as pd
import os
import csv

def find_motifs(motif, genome):
	indices = []
	index = -1
	while True:
		index = genome.find(motif, index + 1)
		if index == -1:
			break
		indices.append(index)
	return indices

def find_distances_between_motifs(positions, motif_length):
	diff = []
	for pos in range(1,len(positions)):
		diff.append(positions[pos] - positions[pos-1] - motif_length)
	return diff

def find_start_end_positions(positions, diff, distance):
	start = []
	end = []
	for i in range(len(diff)):
		if(diff[i] == distance):
			start.append(positions[i])
			end.append(positions[i+1])
	return start, end

def find_seq_between_motifs(genome, motif, start, end):
	sequences = []
	for num in range(len(start)):
		seq = ""
		for i in range(start[num], end[num]):
			seq += genome[i]
			seq = seq.split(motif)[0]
		sequences.append(seq)
	return sequences

def filter_sequence_by_threshold(sequences, at, gc):
	filtered_sequences_by_threshold = []
	for seq in sequences:
		agct = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
		for j in range(len(seq)):
			letter = seq[j]
			agct[letter] += 1
		for key, val in agct.items():
			agct[key] *= 100
			agct[key] /= len(seq)
		at_content = agct['A'] + agct['T']
		gc_content = agct['G'] + agct['C']
		if at_content >= at and gc_content >= gc:
			filtered_sequences_by_threshold.append(seq)
	return filtered_sequences_by_threshold

def analyse_sequences(sequences, distance):
	analysis = []
	for i in range(distance):
		agct = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
		for j in range(len(sequences)):
			letter = sequences[j][i]
			agct[letter] += 1
		analysis.append(agct)
	return analysis

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

def count_genes_with_motif(path, motif, distance):
	count = 0
	spacer_sequences = []
	gene_ids = []
	gene_sequences = []
	with open(path,'r') as data:
		for line in data:
			gene_id = line.split('\t')[0]
			gene_sequence = line.split('\t')[1]
			positions = find_motifs(motif, gene_sequence)
			diff = find_distances_between_motifs(positions, len(motif))
			start, end = find_start_end_positions(positions, diff, distance)
			sequences = find_seq_between_motifs(gene_sequence, motif, start, end)
			if len(sequences):
				gene_ids.append(gene_id)
				gene_sequences.append(gene_sequence)
				for sequence in sequences:
					spacer_sequences.append(sequence)
				count += 1
	#filtered_sequences_by_threshold = filter_sequence_by_threshold(spacer_sequences, 40, 25)
	return count, spacer_sequences, gene_ids, gene_sequences

def find_spacers_in_genes(path, motif, distance):
	counter = {}
	with open(path, 'r') as  file:
		for line in file:
			gene_id = line.split("\t")[0]
			gene_sequence = line.split("\t")[1]
			positions = find_motifs(motif, gene_sequence)
			diff = find_distances_between_motifs(positions, len(motif))
			sequence_count = []
			for dist in range(distance+1):
				sequence_count.append(0)
				start, end = find_start_end_positions(positions, diff, dist)
				if len(start):
					sequence_count[dist] = len(start)
			counter[gene_id] = sequence_count
	with open("test_counter.csv", "w") as outfile:
   		writer = csv.writer(outfile)
   		writer.writerow(counter.keys())
   		writer.writerows(zip(*counter.values()))
	print ('saving is complete')
	return counter

if __name__ == "__main__":
	motif = "AAAG"
	path = "Disease - All.txt"
	gap_length = []
	gene_count = []
	spacer_count = []
	all_spacer_sequences = []
	filtered_spacer_sequences = []
	smallest = 1
	biggest = 31
	for distance in range(smallest, biggest):
		if distance == 17:
			continue
		count, spacer_sequences, gene_ids, gene_sequences = count_genes_with_motif(path, motif, distance)
		gap_length.append(distance)
		gene_count.append(count)
		spacer_count.append(len(spacer_sequences))
		for seq in spacer_sequences:
			all_spacer_sequences.append(seq)
		analysis = analyse_sequences(spacer_sequences, distance)
		newpath = path.split('.txt')[0]
		if not os.path.exists(newpath):
			os.makedirs(newpath)
		''''
		print("Writing gene IDs and there sequences for distance = %d" % distance)
		with open(path.split('.txt')[0]+"/%s,%d-genes.txt" % (motif, distance), 'w') as file:
			for i in range(len(gene_ids)):
				file.write(gene_ids[i]+' \n'+gene_sequences[i]+'\n')
		print("Writing spacer sequences for distance = %d" % distance)
		with open(path.split('.txt')[0]+"/%s,%d-spacer_sequences.txt" % (motif, distance), 'w') as file:
			total_sequences = len(spacer_sequences)
			for sequence in spacer_sequences:
				file.write(sequence)
				file.write('\n')
		print("Writing spacer sequence analysis for distance = %d" % distance)
		with open(path.split('.txt')[0]+"/%s,%d-spacer_base_analyses.csv" % (motif, distance), 'w') as file:
			pos=1
			file.write("Position,A,G,C,T\n")
			for item in analysis:
				file.write("%d," % (pos))
				pos += 1
				for key, val in item.items():
					if total_sequences > 0:
						file.write("%s: %f, " % (key, 100*val/total_sequences))
				file.write('\n')
	print("Writing spacer frequencies")
	df = pd.DataFrame({
		"Gap length": gap_length,
		"Gene count": gene_count,
		"Spacer sequence count": spacer_count
		})
	df.to_csv(path.split('.txt')[0]+".csv", index=False)
	print("Writing flanking analysis")
	overall_analysis = analyse_variable_length_sequences(all_spacer_sequences, biggest)
	with open(path.split('.txt')[0]+"/%s-flank.csv" % (motif), 'w') as file:
		pos=1
		file.write("Position,A,G,C,T\n")
		for item in overall_analysis:
			file.write("%d," % (pos))
			pos += 1
			for key, val in item.items():
				file.write("%s: %f, " % (key, val))
			file.write('\n')
	'''
	counter = find_spacers_in_genes(path, motif, 30)
	print(counter)