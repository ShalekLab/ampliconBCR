import os
import pandas as pd
import json
import argparse

def main(args):
	
	align_heavy = args.heavy.split(',')
	print(align_heavy)
	align_light = args.light.split(',')
	print(align_light)
	align_combined = list()

	for heavy in align_heavy:
		sample = os.path.basename(heavy).replace('_corr_align.txt', '')
		for light in align_light:
			if sample == os.path.basename(light).replace('_corr_align.txt', ''):
				align_combined.append({
					"sample": sample,
					"heavy": heavy,
					"light": light
				})
				break
		else:
			raise Exception("Couldn't find corresponding light results for heavy object")
	print("aligncombined", align_combined)
	for results in align_combined:
		print("result:", results)
		light = pd.read_csv(results["light"], sep='\t', header=0)
		heavy = pd.read_csv(results["heavy"], sep='\t', header=0)
		combined = pd.concat([light, heavy])
		output_filename = os.path.join(args.output_path, results["sample"] + "_combined_corr_align.txt")
		combined.to_csv(output_filename, sep='\t', header=True, index=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-o', "--output-path", 
		type=str,
		help="The path to which the files from this program will output."
	)
	parser.add_argument(
		"-H", "--heavy", 
		type=str,
		help="Comma-delimited list of heavy files"
	)
	parser.add_argument(
		"-L", "--light", 
		type=str,
		help="Comma-delimited list of light files"
	)
	args = parser.parse_args()
	main(args)