"""
 Wrapper for Promoter_Calculator.py
"""

import time
import math
from copy import copy
from .promoter_calculator import Promoter_Calculator, PromoCalcResults
from dataclasses import dataclass

#@dataclass(frozen=True)
class PromoCalcResult(PromoCalcResults):
    pass


def promoter_calculator(sequence, quiet=False):
    begin = time.time()
    #sequence = "".join([random.choice(['A','G','C','T']) for x in range(100000)])
    sequence = sequence.upper()
    output = []

    # Helper function to get raw output
    def _run_calculator(target, start_pos, end_pos, fraction):
        calc = Promoter_Calculator()
        calc.run(target, TSS_range = [0, len(target)])
        rev_results = calc.output()['Reverse_Predictions_per_TSS']
        fwd_results = calc.output()['Forward_Predictions_per_TSS']
        calculator_result = []
        for i in fwd_results.keys():
            fwd_results[i]['TSS_name'] = f"Fwd{fwd_results[i]['TSS']}"
            fwd_results[i]['TSS'] = start_pos + fwd_results[i]['TSS']
            fwd_results[i]['UP_position'][0] = start_pos + fwd_results[i]['UP_position'][0]
            fwd_results[i]['UP_position'][1] = start_pos + fwd_results[i]['UP_position'][1]
            fwd_results[i]['hex35_position'][0] = start_pos + fwd_results[i]['hex35_position'][0]
            fwd_results[i]['hex35_position'][1] = start_pos + fwd_results[i]['hex35_position'][1]
            fwd_results[i]['spacer_position'][0] = start_pos + fwd_results[i]['spacer_position'][0]
            fwd_results[i]['spacer_position'][1] = start_pos + fwd_results[i]['spacer_position'][1]
            fwd_results[i]['hex10_position'][0] = start_pos + fwd_results[i]['hex10_position'][0]
            fwd_results[i]['hex10_position'][1] = start_pos + fwd_results[i]['hex10_position'][1]
            fwd_results[i]['disc_position'][0] = start_pos + fwd_results[i]['disc_position'][0]
            fwd_results[i]['disc_position'][1] = start_pos + fwd_results[i]['disc_position'][1]
            fwd_results[i]['drop'] = False
            if type(fraction) != str:
                fwd_results[i]['fraction'] = fraction + 1
            fwd_results[i]['strand'] = '+'
            calculator_result.append(copy(fwd_results[i]))
        for i in rev_results.keys():
            rev_results[i]['TSS_name'] = f"Rev{rev_results[i]['TSS']}"
            rev_results[i]['TSS'] = end_pos - rev_results[i]['TSS']
            rev_results[i]['UP_position'][0] = end_pos - rev_results[i]['UP_position'][0]
            rev_results[i]['UP_position'][1] = end_pos - rev_results[i]['UP_position'][1]
            rev_results[i]['hex35_position'][0] = end_pos - rev_results[i]['hex35_position'][0]
            rev_results[i]['hex35_position'][1] = end_pos - rev_results[i]['hex35_position'][1]
            rev_results[i]['spacer_position'][0] = end_pos - rev_results[i]['spacer_position'][0]
            rev_results[i]['spacer_position'][1] = end_pos - rev_results[i]['spacer_position'][1]
            rev_results[i]['hex10_position'][0] = end_pos - rev_results[i]['hex10_position'][0]
            rev_results[i]['hex10_position'][1] = end_pos - rev_results[i]['hex10_position'][1]
            rev_results[i]['disc_position'][0] = end_pos - rev_results[i]['disc_position'][0]
            rev_results[i]['disc_position'][1] = end_pos - rev_results[i]['disc_position'][1]
            rev_results[i]['strand'] = '-'
            rev_results[i]['drop'] = False
            if type(fraction) != str:
                rev_results[i]['fraction'] = fraction + 1
            calculator_result.append(copy(rev_results[i]))
        for result in calculator_result:
            result['length'] = len(result['promoter_sequence'])

        to_keep = []
        length = len(calculator_result)
        for i in range(length):
            if not quiet:
                percent_complete = i/length*100
                if percent_complete % 10 == 0:
                    print(f"{percent_complete:.2f}%")
            promotor = calculator_result.pop(0)
            if 'drop' not in promotor.keys():
                promotor['drop'] = False
            for i, other_promotor in enumerate(calculator_result):
                if promotor['fraction'] != other_promotor['fraction'] and promotor['fraction'] != 'junction' and other_promotor['fraction'] != 'junction':
                    continue
                if promotor['strand'] != other_promotor['strand']:
                    continue
                for position_type in ['UP_position', 'hex35_position', 'spacer_position', 'hex10_position', 'disc_position']:
                    if promotor[position_type] == other_promotor[position_type]:
                        if promotor['dG_total'] < other_promotor['dG_total']:
                            promotor['drop'] = True
                            break
                        if promotor['dG_total'] > other_promotor['dG_total']:
                            other_promotor['drop'] = True
                            break
                        else:
                            promotor['drop'] = True
                            break
            if promotor['drop'] == False:
                to_keep.append(promotor)
        return to_keep

    # Run raw output of calculator, chunk if necessary
    fraction_length = 2500
    if len(sequence) >= 5*fraction_length:
        if not quiet:
            print('Large input detected. Chunking sequence into smaller pieces.')
        fractions = math.floor(len(sequence)/fraction_length)
        remainder = len(sequence) % fraction_length
        if remainder > 0:
            fractions += 1
        else:
            # remainder onto the last run
            pass
        for i in range(fractions):
            if i == fractions-1:
                end_pos = len(sequence)
            else:
                end_pos = (i+1)*fraction_length+100
            if i == 0:
                start_pos = 0
            else:
                start_pos = i*fraction_length-100
            if not quiet:
                print(f"Fraction {i+1} from {start_pos} to {end_pos}")
            target_sequence = sequence[start_pos:end_pos]
            if not quiet:
                print("Fraction:", len(target_sequence), "Sequence: ", len(sequence))

            result = _run_calculator(target_sequence, start_pos, end_pos, i)
            output.extend(result)

    else:
        result = _run_calculator(sequence, 0, len(sequence), 0)
        output.extend(result)

    # If the DNA is circular, examine the junction
    ciruclar = False
    if ciruclar:
        if not quiet:
            print("Circular DNA detected. Examining junction.")
        # Get the junction
        junction = sequence[-100:] + sequence[:100]
        result = _run_calculator(junction, len(sequence)-100, len(sequence)+100, "junction")
        # Filter out promotors that don't cross the junction
        for result in result:
            if result['strand'] == '+':
                bounds = [result['UP_position'][0], result['disc_position'][1]]
            if result['strand'] == '-':
                bounds = [result['disc_position'][0], result['UP_position'][1]]
            if bounds[0] <= 100 and bounds[1] >= 101:
                output.extend(result)

    if not quiet:
        print("Removing duplicates.")
    output = [i for n, i in enumerate(output) if i not in output[n + 1:]] # remove duplicates
    # Find the best results
    output.sort(key=lambda x: x['Tx_rate'])
    output = output[::-1]
    return output
