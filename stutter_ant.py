#!/usr/bin/env python3

'''SCRIPT DESCRIPTION
'''

# Parameters
import argparse
parser = argparse.ArgumentParser()
requireparser = parser.add_argument_group('required arguments')
require.parser.add_argument('-gb', '--genbank', help='Input genbank file', required=True)
require.parser.add_argument('-n', '--gene_name', help='Gene name in genbank file to search for', required=True)
parser.add_argument('-f', '--force', help='Force overwrite output files', default = 'True')
args = parser.parse_args()

# Import libraries
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
import os
import gzip
import pandas as pd

gb_file = args.genbank
gene = args.gene_name

##### Force overwrite output files
if args.force == "True":
	# remove files to overwrite data before starting loop
	outputfiles = ["single_copy_rpoc.fa", \
		"short_rpoc.fa", \
		"long_rpoc.fa", \
		"multi_copy_rpoc.fa", \
		"multi_long_rpoc.fa", \
		"multi_short_rpoc.fa"]

	for f in outputfiles:
		if os.path.isfile(f):
			os.remove(f)

######### SINGLE COPY PROCESSING
# first pull all single copy genes into a single fasta file
for gb_record in SeqIO.parse(gzip.open(gb_file, "rt"), "genbank"):
	for feature in gb_record.features:
		if feature.type == "CDS" and 'gene' in feature.qualifiers:
			if gene in feature.qualifiers['gene'][0] and "_" not in feature.qualifiers['gene'][0]:
				start = feature.location.start
				end = feature.location.end

				# is the single copy gene within expected size range?
				if 4118 + 902 >= len(gb_record.seq[start:end]) >= 4118 - 902:
					with open("single_copy_rpoc.fa", "a") as f:
						print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "single_copy", \
							gb_record.name, \
							feature.qualifiers['product'][0], \
							len(gb_record.seq[start:end]), \
							start, end, \
							gb_record.seq[start:end]), file = f)
				elif 4118 + 902 <= len(gb_record.seq[start:end]):
					with open("long_rpoc.fa", "a") as f:
						print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "single_copy", \
							gb_record.name, \
							feature.qualifiers['product'][0], \
							len(gb_record.seq[start:end]), \
							start, end, \
							gb_record.seq[start:end]), file = f)			
				elif 4118 - 902 >= len(gb_record.seq[start:end]):
					with open("short_rpoc.fa", "a") as f:
						print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "single_copy", \
							gb_record.name, \
							feature.qualifiers['product'][0], \
							len(gb_record.seq[start:end]), \
							start, end, \
							gb_record.seq[start:end]), file = f)	

######### MULTI COPY PROCESSING
multi = []

for gb_record in SeqIO.parse(gzip.open(gb_file, "rt"), "genbank"):
	for feature in gb_record.features:
		if feature.type == "CDS" and 'gene' in feature.qualifiers:
			if gene in feature.qualifiers['gene'][0] and "_" in feature.qualifiers['gene'][0]:
				multi += [[gb_record.name, \
				feature.qualifiers['gene'][0], \
				feature.location.start, \
				feature.location.end, \
				len(gb_record.seq[feature.location.start:feature.location.end])]]

	# convert to pandas dataframe
	df = pd.DataFrame(multi, columns = ['ContigID', 'annotated_as', 'start', 'end', 'length'])
	# group into length and coordinates for each contig
	df1 = df.groupby("ContigID")["length"].sum().reset_index().merge(df.groupby("ContigID")["start"].min().reset_index())
	df = df1.merge(df.groupby("ContigID")["end"].max().reset_index())

	# get actual length from coordinates
	df.insert(4, "coordinate_len", df['end'] - df['start'])

# Use length profiles to pull genes by coordinate values
for i in range(len(df)):
	for gb_record in SeqIO.parse(gzip.open(gb_file, "rt"), "genbank"):
		if gb_record.id == df.loc[i].at["ContigID"]:
			seq = gb_record.seq 
			start = df.loc[i].at["start"]
			end = df.loc[i].at["end"]

			# is the multi copy gene within expected size range?
			if 4118 + 902 >= len(gb_record.seq[start:end]) >= 4118 - 902:
				with open("multi_copy_rpoc.fa", "a") as f:
					print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "multi_copy", \
						gb_record.name, \
						feature.qualifiers['product'][0], \
						len(gb_record.seq[start:end]), \
						start, end, \
						gb_record.seq[start:end]), file = f)
			# long multi copy reads
			elif 4118 + 902 <= len(gb_record.seq[start:end]):
				with open("multi_long_rpoc.fa", "a") as f:
					print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "multi_copy", \
						gb_record.name, \
						feature.qualifiers['product'][0], \
						len(gb_record.seq[start:end]), \
						start, end, \
						gb_record.seq[start:end]), file = f)			
			elif 4118 - 902 >= len(gb_record.seq[start:end]):
				with open("multi_short_rpoc.fa", "a") as f:
					print(">%s|%s|%s|%s|%i|%s:%s\n%s" % (feature.qualifiers['gene'][0], "multi_copy", \
						gb_record.name, \
						feature.qualifiers['product'][0], \
						len(gb_record.seq[start:end]), \
						start, end, \
						gb_record.seq[start:end]), file = f)		










