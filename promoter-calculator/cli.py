"""
Command line interface for the promoter-calculator.
"""

import argparse
import os
import csv
from copy import copy
from .wrapper import promoter_calculator
from dataclasses import asdict

def print_promo_calculator(output):
    print("#------------------------------------------------------#")
    print("Highest expressing promotors: \n")
    if not output:
        print("No promotors found.")
        return
    if len(output) > 25:
        num_results_to_print = 25
    else:
        num_results_to_print = len(output)
    for i in range(0, num_results_to_print):
        print(f"Highest Tx rate {i+1}: {output[i]['Tx_rate']} | Length: {output[i]['length']} | Position: {output[i]['UP_position']}{output[i]['strand']} | Fr {output[i]['fraction']}")


def main():
    """ CLI entry point """
    parser = argparse.ArgumentParser(description=f'Salis Lab Promoter Calculator')

    parser.add_argument(
        '-i', '--input',
        action='store',
        metavar='str/filepath',
        dest='i',
        required=False,
        type=str,
        help="Input DNA/RNA sequence.",
    )

    parser.add_argument(
        '-v', '--verbosity',
        action='store',
        metavar='int',
        dest='v',
        required=False,
        type=int,
        help="Verbosity level. 0 = silent, 1 = progress bar, 2 = detailed progress.",
        default=1
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        metavar='filepath',
        dest='o',
        required=False,
        type=str,
        help="Output file path.",
    )

    parser.add_argument(
        '-j', '--threads',
        action='store',
        metavar='int',
        dest='j',
        required=False,
        type=int,
        help="Number of threads to use. (Min 1)",
        default=1
    )

    parser.add_argument('-c',
        action='store_true',
        dest='c',
        required=False,
        help='Indecate the input is a circular sequence.',
        default=False)
    

    options = parser.parse_args()

    # Show help if no input is provided
    if not options.i:
        parser.print_help()
        return

    # Show help if threads is invalid
    if options.j < 1:
        parser.print_help()
        return

    if os.path.isfile(options.i):
        print("Input file is a file. Reading sequence from file.")
        with open(options.i, 'r') as f:
            sequence_lines = f.readlines()
            sequence = []
            for line in sequence_lines:
                line = line.strip()
                if ">" in line:
                    continue
                sequence.append(line)
            sequence = "".join(sequence)
            print(f"Sequence length: {len(sequence)}")
    else:
        sequence = options.i

    output = promoter_calculator(sequence, verbosity = options.v, threads = options.j, circular=options.c)
    if len(output) == 0:
        print("No promotors found.")
        return

    if not options.o:
        print_promo_calculator(output)
    else:
        # Make a CSV with all the hits
        with open(options.o, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(asdict(output[0]).keys())
            for row in output:
                writer.writerow(asdict(row).values())

if __name__ == "__main__":
    main()
