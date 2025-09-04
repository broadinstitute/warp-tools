import argparse
import scipy
from scipy.io import mmwrite

def process_files(input_files,output_file):
    # Read the first file of the input files to get the header later
    file1 = input_files[0]
    print(file1)   
    output_name=output_file

    matrix_dict = {}
    for file in input_files:
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

    print("Writing matrix market header, number of rows, columns and elements")
    with open(file1, 'r') as f, open("new_mtx.mtx", 'w') as outfile:
        for line in f:
            if line.startswith('%'):
                outfile.write(line)
            else:
                value1, value2 = line.strip().split()[0:2]
                outfile.write(" ".join([str(value1), str(value2), str(elements)]))
                break

    print("Writing the keys and values of dictionary to outfile")
    with open("new_mtx.mtx", 'a') as outfile:
        with open("new_mtx.mtx", 'r') as infile:
            for _ in range(3):
                next(infile)

        outfile.write('\n')

        for key, value in matrix_dict.items():
            split_key = key.split('-')
            line = ' '.join(split_key + value) + '\n'
            outfile.write(line)

    print("Read output mtx file into matrix with scipy to correct formatting")
    coo_mat = scipy.io.mmread('new_mtx.mtx')
    print("Writing output mtx file in correct format")
    mmwrite(output_name, coo_mat)

def main():
    parser = argparse.ArgumentParser(description='Process matrix files')
    # Adding parser argument using nargs to allow for array
    parser.add_argument('input_files', nargs='+', help='List of input matrix files')
    parser.add_argument('output_file', help='Output file name')

    args = parser.parse_args()
    process_files(args.input_files, args.output_file)

if __name__ == "__main__":
    main()