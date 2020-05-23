import sys
from Bio import SeqIO
import statistics as stat

def vista_statistics():
	pos_len = []
	for seq_record in SeqIO.parse("VistaDataset/human_positive.fasta", "fasta"):
		pos_len.append(len(seq_record.seq))

	neg_len = []
	for seq_record in SeqIO.parse("VistaDataset/human_negative.fasta", "fasta"):
		neg_len.append(len(seq_record.seq))

	print("Positive Samples: ", len(pos_len))
	print("Average Length: ", stat.mean(pos_len))
	print("Standard Deviation: ", stat.stdev(pos_len))
	
	print("Negative Samples: ", len(neg_len))
	print("Average Length: ", stat.mean(neg_len))
	print("Standard Deviation: ", stat.stdev(neg_len))

vista_statistics()

# def read_vista():
# 	seq_0 = ""
# 	chromosome_id = 0
# 	start_idx = 0
# 	end_idx = 0
# 	for seq_record in SeqIO.parse("VistaDataset/human_positive.fasta", "fasta"):
# 	    id = seq_record.id
# 	    chromosome_id = id.split(":")[0]
# 	    chromosome_id = chromosome_id.split("|")[1]

# 	    range_idx = id.split(":")[1]
# 	    range_idx = range_idx.split("-")
# 	    start_idx = range_idx[0]
# 	    end_idx = range_idx[1]

# 	    seq_0 = seq_record.seq
# 	    seq_0 = seq_0.lower()

# 	    print(chromosome_id)
# 	    print(seq_record.id)
# 	    print(start_idx)
# 	    print(end_idx)
# 	    print(seq_0)
# 	    break

# 	return seq_0


# def read_genome(seq):
# 	for seq_record in SeqIO.parse("HumanGenome/genomic_data.fna", "fasta"):
# 		description = seq_record.description
# 		sequence = seq_record.seq.lower()
# 		# description = description.split(",")[0].split()
# 		# chromosome_id = description[4]
		
# 		if " 16 " in description:
# 			print(seq_record.description)
# 			print(len(seq_record.seq))
# 			if seq in sequence:
# 				print("----------------------- FOUND -----------------------")
# 			# break

# read_genome(seq)