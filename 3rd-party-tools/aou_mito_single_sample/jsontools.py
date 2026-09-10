import json
import argparse
import os

def parse_var(s):
    """
    Parse a key, value pair, separated by '='
    That's the reverse of ShellArgs.

    On the command line (argparse) a declaration will typically look like:
        foo=hello
    or
        foo="hello world"
    """
    items = s.split('=')
    key = items[0].strip() # we remove blanks around keys, as is logical
    if len(items) > 1:
        # rejoin the rest:
        value = '='.join(items[1:])
    return (key, value)


def parse_vars(items, type):
    """
    Parse a series of key-value pairs and return a dictionary
    """
    d = {}

    if items:
        for item in items:
            key, value = parse_var(item)
            d[key] = type(value)
    return d

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, required=True)
parser.add_argument("--set",
                        metavar="KEY=VALUE",
                        nargs='+',
                        help="Set a number of key-value pairs "
                             "(do not put spaces before or after the = sign). "
                             "If a value contains spaces, you should define "
                             "it with double quotes: "
                             'foo="this is a sentence". Note that '
                             "values are always treated as strings.")
parser.add_argument("--set-int",
                        metavar="KEY=VALUE",
                        nargs='+',
                        help="Set a number of key-value pairs "
                             "(do not put spaces before or after the = sign). "
                             "If a value contains spaces, you should define "
                             "it with double quotes: "
                             'foo="this is a sentence". Note that '
                             "values are always treated as ints.")
parser.add_argument("--set-float",
                        metavar="KEY=VALUE",
                        nargs='+',
                        help="Set a number of key-value pairs "
                             "(do not put spaces before or after the = sign). "
                             "If a value contains spaces, you should define "
                             "it with double quotes: "
                             'foo="this is a sentence". Note that '
                             "values are always treated as ints.")

if __name__ == "__main__":
    args = parser.parse_args()
    data_in = parse_vars(args.set, type=str)
    data_in_int = parse_vars(args.set_int, type=int)
    data_in_float = parse_vars(args.set_float, type=float)
    data_in.update(data_in_int)
    data_in.update(data_in_float)
    keys_in = list(data_in.keys())
    
    if os.path.exists(args.path):
        with open(args.path, 'r') as json_file:    
            file_of_interest = json.load(json_file)
        k_in_file = list(file_of_interest.keys())
        new_in_old = all([x in k_in_file for x in keys_in])
        old_in_new = all([x in keys_in for x in k_in_file])
        if not new_in_old or not old_in_new:
            raise ValueError('ERROR: json file contains different keys from input.')
        
        for k in keys_in:
            file_of_interest[k].append(data_in[k])
    else:
        file_of_interest = {k: [data_in[k]] for k in keys_in}
    
    with open(args.path, 'w') as new_file:
        json.dump(file_of_interest, new_file)