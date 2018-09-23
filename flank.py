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

def analyse_sequences(sequences, total):
	analysis = []
	agct = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
	for j in range(len(sequences)):
		letter = sequences[j][i]
		agct[letter] += 1
	analysis.append(agct)
	for item in analysis:
		for key, val in item.items():
			item[key] *= 100
			item[key] /= total
	return analysis

if __name__ == "__main__":
	motif = "AAAG"
	motif_length = len(motif)
	path = "Disease - All.txt"
	distance = 30
	flank_left = []
	flank_right = []
	total_left = 0
	total_right = 0
	with open(path, 'r') as file:
		for line in file:
			gene_id = line.split('\t')[0]
			gene_sequence = line.split('\t')[1]
			positions = find_motifs(motif, gene_sequence)
			diff = find_distances_between_motifs(positions, len(motif))
			for dist in range(2,distance+1):
				start, end = find_start_end_positions(positions, diff, dist)
				for i in range(len(start)):
					if gene_sequence[start[i] - 1] in {'A','G','C','T'}:
						total_left += 1
						flank_left.append(gene_sequence[start[i] - 1])
					if gene_sequence[end[i] + motif_length] in {'A','G','C','T'}:
						total_right += 1
						flank_right.append(gene_sequence[end[i] + motif_length])
	print("Left flank")
	print(analyse_sequences(flank_left, total_left))
	print("Right flank")
	print(analyse_sequences(flank_right, total_right))