# 2020
# Tim Webster
# University of Utah

from __future__ import print_function
import argparse


def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input", required=True,
		help="")

	parser.add_argument(
		"--dp", required=True, type=float,
		help="")

	parser.add_argument(
		"--mq", required=True, type=float,
		help="")

	parser.add_argument(
		"--gq", required=True, type=float,
		help="")

	parser.add_argument(
		"--sample1_name", required=True,
		help="")

	parser.add_argument(
		"--sample2_name", required=True,
		help="")

	parser.add_argument(
		"--output", required=True,
		help="Path to and name of output file")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()
	print(type(args.dp))
	passing1 = 0
	passing2 = 0
	passingboth = 0
	shared = 0
	dropout1 = 0
	dropout2 = 0
	ambiguous = 0
	with open(args.input, "r") as f:
		for line in f:
			one_pass = False
			two_pass = False
			line1 = line.strip().split(",")
			if "position" in line1:
				continue
			one = [line1[x] for x in [5, 7, 9, 11]]
			two = [line1[x] for x in [6, 8, 10, 12]]

			if "." not in one:
				one[1] = float(one[1])
				one[2] = float(one[2])
				one[3] = float(one[3])
				if one[1] >= args.dp:
					if one[2] >= args.mq:
						if one[3] >= args.gq:
							one_pass = True
							passing1 += 1
			if "." not in two:
				two[1] = float(two[1])
				two[2] = float(two[2])
				two[3] = float(two[3])
				if two[1] >= args.dp:
					if two[2] >= args.mq:
						if two[3] >= args.gq:
							two_pass = True
							passing2 += 1

			if one_pass is True and two_pass is True:
				passingboth += 1
				if one[0] == two[0]:
					shared += 1
				else:
					one_alleles = [one[0][0], one[0][2]]
					two_alleles = [two[0][0], two[0][2]]

					# 1 is heterozygous
					if one_alleles[0] != one_alleles[1]:
						# 2 is heterozygous
						# ambigous because the two samples contain different alleles
						if two_alleles[0] != two_alleles[1]:
							ambiguous += 1
							print("Ambiguous case 1: {}".format(line))
							print(one_alleles)
							print(two_alleles)
						# 2 is homozygous
						else:
							# 2's allele is present in 1
							if two_alleles[0] in one_alleles:
								dropout2 += 1
								print("Dropout 2: {}".format(line))
							# 2 has a different allele than the two present in 1
							else:
								ambiguous += 1
								print("Ambiguous case 2: {}".format(line))
								print(one_alleles)
								print(two_alleles)
					# 1 is homozygous
					else:
						# 2 is heterozygous
						if two[0][0] != two[0][2]:
							if one[0][0] in two_alleles:
								dropout1 += 1
								print("Dropout 1: {}".format(line))
							else:
								ambiguous += 1
								print("Ambiguous case 3: {}".format(line))
								print(one_alleles)
								print(two_alleles)
						# 2 is homozygous, must be different allele
						else:
							ambiguous += 1
							print("Ambiguous case 4: {}".format(line))
							print(one_alleles)
							print(two_alleles)

	with open(args.output, "w") as o:
		o.write("Passing {}: {}\n".format(args.sample1_name, passing1))
		o.write("Passing {}: {}\n".format(args.sample2_name, passing2))
		o.write("Passing both: {}\n".format(passingboth))
		o.write("Shared: {}\n".format(shared))
		o.write("Dropout {}: {}\n".format(args.sample1_name, dropout1))
		o.write("Dropout {}: {}\n".format(args.sample2_name, dropout2))
		o.write("Ambiguous: {}\n".format(ambiguous))


if __name__ == "__main__":
	main()
