from Bio import SeqIO
import re

def take_input():
	motif = input("Enter a motif: ").upper()

	while True:
		verbose = input("Do you want to see the positions of the motif? (y/n) ").upper()

		if verbose == "Y" or verbose == "N":
			break

		print("Please enter y for yes and n for no.")

	return motif, verbose

def find_indices_of_motifs(positions):
	print(positions)

def find_motifs(motif, genome):
	"""
	Returns list of indices where motif begins in genome

	>>> find_substring('me', "The cat says meow, meow")
	[13, 19]
	"""
	indices = []
	index = -1  # Begin at -1 so index + 1 is 0
	while True:
		# Find next index of motif, by starting search from index + 1
		index = genome.find(motif, index + 1)
		if index == -1:
			break  # All occurrences have been found
		indices.append(index)
	return indices

def ask_to_continue():
	cont = ""
	while True:
		cont = input("Do you want to find sequences between motifs? (y/n) ").upper()

		if cont == "Y" or cont == "N":
			break

		print("Please enter y for yes and n for no.")
	return cont

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

def ask_to_save():
	save = ""
	while True:
		save = input("Do you want to save the sequences and the frequencies of bases to a text file? (y/n) ").upper()

		if save == "Y" or save == "N":
			break

		print("Please enter y for yes and n for no.")

	return save

def save_to_txt(sequences, motif, distance, analysis):
	thefile = open('motif=%s, distance=%d, number=%d.txt' % (motif, distance, len(sequences)), 'w')

	thefile.write("The %d sequences are:\n\n" % len(sequences))

	for item in sequences:
		thefile.write("%s\n" % item)

	pos = 1
	thefile.write("\nThe frequency of bases at each position are:\n\n")

	for item in analysis:
		thefile.write('\n')
		thefile.write("At position %d, " % pos)
		pos += 1
		for key, val in item.items():
			thefile.write("%s: %d, " % (key, val))

def run_motif_finder(paths):
	for path in paths:
		with open(path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				genome = record.seq

				motif, verbose = take_input()

				motif_length = len(motif)

				positions = find_motifs(motif, genome)

				if verbose == "Y":
					find_indices_of_motifs(positions)

				print("\n----------------------------------------------------------")
				print("The motif %s occurs %d times in the genome supplied." % (motif, len(positions)))
				print("----------------------------------------------------------\n")

				cont = ask_to_continue()

				if cont == "N":
					break

				distance = take_distance()

				diff = find_distances_between_motifs(positions, motif_length)

				start, end = find_start_end_positions(positions, diff, distance)
				sequences = find_seq_between_motifs(genome, motif, start, end)

				print("\n----------------------------------------------------------")
				print("There are %d sequences that are %dbp long between motifs." % (len(sequences), distance))
				print("----------------------------------------------------------\n")

				analysis = analyse_sequences(sequences, distance)

				save = ask_to_save()

				if save == "Y":
					save_to_txt(sequences, motif, distance, analysis)

if __name__ == "__main__":
	paths = ["at-genome/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa"]
	run_motif_finder(paths)
