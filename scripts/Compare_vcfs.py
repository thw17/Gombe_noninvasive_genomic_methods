# 2020
# Tim Webster
# University of Utah

from __future__ import print_function
import argparse
import collections
import gzip
import pandas as pd


def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--vcf1", required=True,
		help="")

	parser.add_argument(
		"--vcf2", required=True,
		help="")

	parser.add_argument(
		"--indel_window", type=int, default=5,
		help="")

	parser.add_argument(
		"--output_file", required=True,
		help="Path to and name of output file")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	removed_list = []
	dict1 = collections.OrderedDict()
	# gt, dp, mq, gq
	dict1["header"] = [
		"chrom", "position", "ref", "alt", "alt2",
		"gt1", "dp1", "mq1", "gq1",
		"gt2", "dp2", "mq2", "gq2"]

	print("Processing vcf1")
	with gzip.open(args.vcf1, 'rt', encoding='utf-8') as f:
		for line in f:
			# Ignore header
			if line[0] == "#":
				continue
			record = line.split()
			chrom = record[0]
			# Looking for biallelic snps
			if len(record[3]) != 1:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + len(record[3]) + 1 + args.indel_window):
					removed_list.append((chrom, i))
				continue

			alts = record[4].split(",")
			alts = [x for x in alts if x != "<NON_REF>"]
			if len(alts) not in [1, 2]:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + 1 + args.indel_window):
					removed_list.append((chrom, i))
				continue
			lengths = [len(x) for x in alts]
			if all(x == 1 for x in lengths) is False:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + 1 + args.indel_window):
					removed_list.append((chrom, i))

				continue

			location = (record[0], record[1])
			qual = record[5]
			ref = record[3]
			if len(alts) == 1:
				alt = alts[0]
				alt2 = "."
			else:
				alt = alts[0]
				alt2 = alts[1]

			info = record[7].split(";")
			if info == ["."]:
				mq = "."
				dp = "."
			else:
				dp = None
				mq = None
				raw_mq = None
				for entry in info:
					if entry[0:3] == "DP=":
						dp = entry.split("=")[1]
					if entry[0:3] == "MQ=":
						mq = entry.split("=")[1]
					elif entry[0:7] == "RAW_MQ=":
						raw_mq = entry.split("=")[1]

			if raw_mq is not None:
				mq = (float(raw_mq) / float(dp)) ** (1.0 / 2.0)

			if mq is None:
				print("No MQ. Record:\n{}".format(record))
				continue

			format = record[9].split(":")
			if len(format) < 4:
				print("No FORMAT. Record:\n{}".format(record))
				continue
			gt = format[0]
			gq = format[3]

			idx1 = int(gt[0])
			idx2 = int(gt[2])

			allele_list = [ref] + alts
			gt = "{}/{}".format(allele_list[idx1], allele_list[idx2])

			dict1[location] = [
				location[0], int(location[1]), ref, alt, alt2, gt, dp, mq, gq]

	print("Processing vcf2")
	with gzip.open(args.vcf2, 'rt', encoding='utf-8') as f:
		for line in f:
			# Ignore header
			if line[0] == "#":
				continue
			record = line.split()
			chrom = record[0]
			# Looking for biallelic snps
			if len(record[3]) != 1:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + len(record[3]) + 1 + args.indel_window):
					removed_list.append((chrom, i))
				continue

			alts = record[4].split(",")
			alts = [x for x in alts if x != "<NON_REF>"]
			if len(alts) not in [1, 2]:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + 1 + args.indel_window):
					removed_list.append((chrom, i))
				continue
			lengths = [len(x) for x in alts]
			if all(x == 1 for x in lengths) is False:
				for i in range(int(record[1]) - args.indel_window, int(record[1]) + 1 + args.indel_window):
					removed_list.append((chrom, i))
				continue

			location = (record[0], record[1])
			qual = record[5]
			ref = record[3]
			if len(alts) == 1:
				alt = alts[0]
				alt2 = "."
			else:
				alt = alts[0]
				alt2 = alts[1]

			info = record[7].split(";")
			if info == ["."]:
				mq = "."
				dp = "."
			else:
				dp = None
				mq = None
				raw_mq = None
				for entry in info:
					if entry[0:3] == "DP=":
						dp = entry.split("=")[1]
					if entry[0:3] == "MQ=":
						mq = entry.split("=")[1]
					elif entry[0:7] == "RAW_MQ=":
						raw_mq = entry.split("=")[1]

			if raw_mq is not None:
				mq = (float(raw_mq) / float(dp)) ** (1.0 / 2.0)

			if mq is None:
				print("No MQ. Record:\n{}".format(record))
				continue

			format = record[9].split(":")
			if len(format) < 4:
				print("No FORMAT. Record:\n{}".format(record))
				continue
			gt = format[0]
			gq = format[3]

			idx1 = int(gt[0])
			idx2 = int(gt[2])

			allele_list = [ref] + alts
			gt = "{}/{}".format(allele_list[idx1], allele_list[idx2])

			# location is present in vcf1 and vcf2
			if location in dict1:
				dict1[location] += [gt, dp, mq, gq]
				if dict1[location][4] == ".":
					if dict1[location][3] == alt:
						dict1[location][4] = alt2
					else:
						dict1[location][4] = alt
			# location is present in vcf2, but not 1
			else:
				dict1[location] = [
					location[0], int(location[1]), ref, alt, alt2,
					".", ".", ".", ".", gt, dp, mq, gq]

	# take care of positions present in vcf1, but not 2
	for key in dict1:
		if dict1[key][0] == "chrom":
			continue
		if len(dict1[key]) != 13:
			dict1[key] = dict1[key] + [".", ".", ".", "."]

	# Remove filtered sites (indel window)
	removed_list = set(removed_list)
	for key in dict1:
		if key in removed_list:
			del dict1[key]

	df = pd.DataFrame.from_dict(dict1, orient='index')
	df.columns = df.iloc[0]
	new_columns = [
		"chrom", "position", "ref", "alt", "alt2",
		"gt1", "gt2", "dp1", "dp2", "mq1", "mq2", "gq1", "gq2"]
	df = df[new_columns]
	df = df.drop(["header"])
	print(df)
	df = df.sort_values(["chrom", "position"])
	df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
	main()
