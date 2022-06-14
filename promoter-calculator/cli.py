"""
Command line interface for the promoter-calculator.
"""

import argparse
import os
from .wrapper import promoter_calculator

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

    options = parser.parse_args()
    if options.i is None:
        print("No input file provided. Exiting.")
        return
    if os.path.isfile(options.i):
        print("Input file is a file. Reading sequence from file.")
        with open(options.i, 'r') as f:
            sequence = f.read()
            sequence = f.read()
    else:
        sequence = options.i
    print("Input sequence:", sequence)

    output = promoter_calculator(sequence)
    print_promo_calculator(output)

if __name__ == "__main__":
    main()
