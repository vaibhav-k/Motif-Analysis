from Bio import SeqIO
import re
import pandas as pd
import os

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

def find_start_end_positions(positions, diff, dist):
	start = []
	end = []
	for i in range(len(diff)):
		if(diff[i] == dist):
			start.append(positions[i])
			end.append(positions[i+1])
	return start, end

def find_frequency_of_motifs_in_genome(paths, motif):
	frequency = 0
	for path in paths:
		with open(path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				genome = record.seq
				positions = find_motifs(motif, genome)
				frequency += len(positions)
	return frequency

def find_spacer_frequency_of_motifs_in_genome(paths, motif, distance):
	spacer_frequency = 0
	motif_length = len(motif)
	for path in paths:
		with open(path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				genome = record.seq
				positions = find_motifs(motif, genome)
				diff = find_distances_between_motifs(positions, motif_length)
				start, end = find_start_end_positions(positions, diff, distance)
				spacer_frequency += len(start)
	return spacer_frequency

def find_total_and_spacer_frequency(paths, motif, distance):
	frequency = 0
	spacer_frequency = []
	for i in range(distance+1):
		spacer_frequency.append(0)
	motif_length = len(motif)
	for path in paths:
		with open(path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				genome = record.seq
				positions = find_motifs(motif, genome)
				frequency += len(positions)
				diff = find_distances_between_motifs(positions, motif_length)
				for i in range(distance+1):
					print(i)
					start, end = find_start_end_positions(positions, diff, i)
					spacer_frequency[i] += len(start)
	return frequency, spacer_frequency

if __name__ == "__main__":
	paths = []
	for file in os.listdir("hv-genome/"):
		paths.append("hv-genome/" + file)

	motif = "AAAG"
	frequency, spacer_frequency = find_total_and_spacer_frequency(paths, motif, 30)
	print("The motif %s occurs %d times in the genome." % (motif, frequency))
	df = pd.DataFrame({
		"Spacer frequency": spacer_frequency,
		})
	print(df)
	df.to_csv("AAAG spacer frequency in HV genome.csv")
