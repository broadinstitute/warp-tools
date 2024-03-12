import glob
import argparse
import scipy
from scipy.io import mmwrite

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process input files and generate a matrix.")
    parser.add_argument("input_files", type=str, help="Path to the files containing input files.")
    parser.add_argument("output_file", type=str, help="Name of the output file for the final matrix.")
    return parser.parse_args()

def process_files(input_files):
    file_list = input_files
    file1 = file_list[0]

    matrix_dict = {}
    for file in file_list:
        with open(file, 'r') as f:
            for _ in range(3):
                next(f)
            for line in f:
                values = line.strip().split()[:3]
                key = "-".join(values[0:2])
                if key in matrix_dict:
                    matrix_dict[key].append(values[2])
                else:
                    matrix_dict[key] = [values[2]]

    elements = len(matrix_dict.keys())
    print("Printing size of elements: ", elements)

    with open(file1, 'r') as f, open("new_mtx.mtx", 'w') as outfile:
        for line in f:
            if line.startswith('%'):
                outfile.write(line)
            else:
                value1, value2 = line.strip().split()[0:2]
                outfile.write(" ".join([str(value1), str(value2), str(elements)]))
                break

    with open("new_mtx.mtx", 'a') as outfile:
        with open("new_mtx.mtx", 'r') as infile:
            for _ in range(3):
                next(infile)

        outfile.write('\n')

        for key, value in matrix_dict.items():
            split_key = key.split('-')
            line = ' '.join(split_key + value) + '\n'
            outfile.write(line)

    coo_mat = scipy.io.mmread('new_mtx.mtx')
    mmwrite(args.output_file, coo_mat)

if __name__ == "__main__":
    args = parse_arguments()
    process_files(args.input_directory)
