from Bio import SeqIO
import os
import pandas as pd

def read_motif():
	motif = input("Enter a motif: ").upper()
	return motif

def find_motifs(motif, genome):
	indices = []
	index = -1
	while True:
		index = genome.find(motif, index + 1)
		if index == -1:
			break
		indices.append(index)
	return indices

def take_distance():
	distance = 0
	while True:
		distance = input("What is the length of sequence between motifs? ")

		if distance.isnumeric():
			distance = int(distance)
			break

		print("Please enter a positive integer.")
	return distance

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

def find_seq_between_motifs(genome, motif, start, end):
	sequences = []
	for num in range(len(start)):
		seq = ""
		for i in range(start[num], end[num]):
			seq += genome[i]
			seq = seq.split(motif)[0]
		sequences.append(seq)
	return sequences

def analyse_sequences(sequences, distance):
	analysis = []
	for i in range(distance):
		agct = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
		for j in range(len(sequences)):
			letter = sequences[j][i]
			agct[letter] += 1
		analysis.append(agct)
	return analysis

def save_to_txt(sequences, motif, distance, proid):
	newpath = motif+", "+str(distance)
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	print("There is/are %d instances of motifs having a gap of %d in promoter: %s." % (len(sequences), distance, proid))
	analysis = analyse_sequences(sequences, distance)
	with open("%s, %d/id=%s.txt" % (motif, distance, proid), 'w') as file:
		for i in range(len(sequences)):
			file.write(proid + "\t" + sequences[i])
			file.write("\n")
		pos=1
		for item in analysis:
			file.write("At position %d, " % pos)
			pos += 1
			for key, val in item.items():
				file.write("%s: %d, " % (key, val))
			file.write('\n')

def save_combined(all_sequences, motif, distance, promoters):
	newpath = motif+", "+str(distance)
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	df = pd.DataFrame({
		"Promoter ID": promoters,
		"Gap sequences": all_sequences
		})
	df.index += 1
	df.to_csv("%s, %d/%s, %d-sequences.csv" % (motif, distance, motif, distance))
	with open("%s, %d/%s, %d-frequency.csv" % (motif, distance, motif, distance), 'w') as file:
		#file.write(str(len(all_sequences))+"\n")
		analysis = analyse_sequences(all_sequences, distance)
		file.write("Position,A,G,C,T\n")
		pos=1
		'''
		for item in all_sequences:
			file.write("%s\n" % item)
		'''
		for item in analysis:
			''''
			file.write("At position %d, " % pos)
			'''
			file.write("%d," % (pos))
			pos += 1
			for key, val in item.items():
				file.write("%d," % (val))
			file.write('\n')

def find_frequency(paths):
	motif = "AAAG"
	motif_length = len(motif)
	for path in paths:
		with open(path, "r") as handle:
			all_sequences = []
			freq = []
			gap = []
			for dist in range(0,31):
				freq.append(0)
				gap.append(dist)
			for record in SeqIO.parse(handle, "fasta"):
				promoter = record.seq
				proid = record.id

				positions = find_motifs(motif, promoter)
				diff = find_distances_between_motifs(positions, motif_length)

				for dist in range(0,31):
					start, end = find_start_end_positions(positions, diff, dist)
					sequences = find_seq_between_motifs(promoter, motif, start, end)
					if len(sequences):
						freq[dist] += 1
			data = pd.DataFrame({
				"Heat Down regulated": freq
				})
			data.to_csv("heatdownreg.csv")

def analyse_promoters(paths):
	#motif = read_motif()
	motif = "AAAG"
	motif_length = len(motif)
	distance = take_distance()
	for path in paths:
		with open(path, "r") as handle:
			all_sequences = []
			promoters = []
			for record in SeqIO.parse(handle, "fasta"):
				promoter = record.seq
				proid = record.id

				positions = find_motifs(motif, promoter)
				diff = find_distances_between_motifs(positions, motif_length)
				start, end = find_start_end_positions(positions, diff, distance)
				sequences = find_seq_between_motifs(promoter, motif, start, end)

				if len(sequences):
					#save_to_txt(sequences, motif, distance, proid)
					for sequence in sequences:
						all_sequences.append(sequence)
						promoters.append(proid)

			save_combined(all_sequences, motif, distance, promoters)

if __name__ == "__main__":
	paths = ["gene-expression-data/HEAT DOWNREG PROM.txt"]
	analyse_promoters(paths)
	#find_frequency(paths)