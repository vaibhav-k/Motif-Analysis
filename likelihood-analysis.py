from Bio import SeqIO
import re
import pandas as pd
import os

def change_to_upper(file):
	with open(file, "r") as inputFile:
		content = inputFile.read()
	with open("PromoterSequences.txt", "w") as outputFile:
		outputFile.write(content.upper())

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

def load_promoter_sequences(path):
	promoters = {}
	with open(path, "r") as ins:
		for line in ins:
			proid = line.split('\t')[0]
			proid = proid.split('.')[0]
			seq = line.split('\t')[1]
			promoters[proid] = seq
	return promoters

def extract_promoter_sequences(path, promoters):
	os.chdir(path)
	files = os.listdir(".")
	all_data = []
	for x in files:
		data = {}
		with open(x, "r") as ids:
			for line in ids:
				if line.split('\n')[0] in promoters:
					proid = line.split('\n')[0]
					data[proid] = promoters[proid].split('\n')[0]
		all_data.append(data)
	os.chdir("../environmental-conditions-sequences")
	for i in range(len(files)):
		with open(files[i], "w") as file:
			data = all_data[i]
			for key, val in data.items():
				file.write("%s\t%s\n" % (key, val))

def find_frequency_of_motifs_genes(paths, motif):
	frequency = []
	names = []
	for path in paths:
		count = 0
		names.append(path.split("/")[1].split(".txt")[0])
		with open(path, "r") as handle:
			for line in handle:
				positions = find_motifs(motif, line)
				count += len(positions)
			frequency.append(count)
	df = pd.DataFrame({
		"Condition": names,
		"AAAG frequency": frequency
		})
	df.to_csv("AAAG in conditions.csv",  index = False)
	return frequency

def find_spacer_frequency_of_motifs_genes(paths, motif, distance):
	spacer_frequency = []
	motif_length = len(motif)
	for path in paths:
		count = 0
		with open(path, "r") as handle:
			for line in handle:
				positions = find_motifs(motif, line)
				diff = find_distances_between_motifs(positions, motif_length)
				start, end = find_start_end_positions(positions, diff, distance)
				count += len(start)
			spacer_frequency.append(count)
	return spacer_frequency

def find_frequency_of_motifs(paths, motif):
	frequency = 0
	for path in paths:
		with open(path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				genome = record.seq
				positions = find_motifs(motif, genome)
				frequency += len(positions)
	return frequency

def find_spacer_frequency_of_motifs(paths, motif, distance):
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

if __name__ == "__main__":
	#paths = ["at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa", "at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa", "at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa", "at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa", "at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa"]
	motif = "AAAG"

	#frequency = find_frequency_of_motifs(paths, motif)
	found_frequency_AAAG = 910702
	''''
	spacer_frequency = []
	for distance in range(0,31):
		occurences = find_spacer_frequency_of_motifs(paths, motif, distance)
		spacer_frequency.append(occurences)
		print("%d\t%d" % (distance, occurences))
	df = pd.DataFrame({
		"Spacer frequency": spacer_frequency,
		})
	df.to_csv("AAAG spacer frequency in at-genome.csv")
	'''
	#change_to_upper("Promoter Seq of A. thaliana genes.txt")
	'''
	promoters_path = "PromoterSequences.txt"
	promoters = load_promoter_sequences(promoters_path)
	newpath = "environmental-conditions-sequences"
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	extract_promoter_sequences("environmental-conditions/", promoters)
	newpaths = []
	names = []
	for file in os.listdir("environmental-conditions-sequences-without-id/"):
		newpaths.append("environmental-conditions-sequences-without-id/" + file)
		names.append(file.split(".txt")[0])
	#find_frequency_of_motifs_genes(newpaths, motif)
	df = pd.DataFrame({
		"Conditions": names
		})
	spacer_frequency = []
	for distance in range(0,31):
		occurences = find_spacer_frequency_of_motifs_genes(newpaths, motif, distance)
		spacer_frequency.append(occurences)
		df["%d" % distance] = occurences
	df.to_csv("AAAG spacer frequency conditions.csv", index = False)
	'''
