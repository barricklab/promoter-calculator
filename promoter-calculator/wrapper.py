"""
 Wrapper for Promoter_Calculator.py
"""

import time
import math
import importlib
from copy import copy
from .promoter_calculator import Promoter_Calculator, PromoCalcResults
from dataclasses import dataclass
from typing import Callable
from .util import timer


if importlib.util.find_spec("progress") is not None:
    from progress.bar import Bar
else:
    Bar = None

def _progress_callback(tasks: int, description:str) -> Callable:
    # Update the progress bar
    if Bar is not None:
        bar = Bar(f'{description}: ', max=tasks)
    else:
        return None
    max = tasks
    itterations = 0

    def _callback():
        nonlocal itterations
        itterations += 1
        nonlocal max
        if itterations == max:
            bar.next()
            bar.finish()
            return None

        bar.next()
        return None
    return _callback

#@dataclass(frozen=True)
class PromoCalcResult(PromoCalcResults):
    pass


def promoter_calculator(sequence, threads=1, verbosity=1, callback=None, circular=False):
    begin = time.time()
    #sequence = "".join([random.choice(['A','G','C','T']) for x in range(100000)])
    sequence = sequence.upper()
    output = []


    # Helper function to get raw output
    def _run_calculator(target, start_pos, end_pos, fraction, callback=None):
        calc = Promoter_Calculator(threads=threads, verbosity=verbosity)
        calc.run(target, TSS_range = [0, len(target)], callback=callback)
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

        # Filter out overlapping promoters

        calculator_result.sort(key=lambda x: x['TSS'])
        SEARCH_RANGE = 10  

        to_keep = []
        length = len(calculator_result)
        if length == 0:
            return [], None, None
        average = sum([x['dG_total'] for x in calculator_result]) / length

        # Legacy overlap filter
        for i in range(length):
            promotor = calculator_result.pop(0)
            if 'drop' not in promotor.keys():
                promotor['drop'] = False
            for _, other_promoter in enumerate(calculator_result):
                if other_promoter['TSS'] > promotor['TSS'] + SEARCH_RANGE:
                    continue
                elif promotor['fraction'] != other_promoter['fraction'] and promotor['fraction'] != 'junction' and other_promoter['fraction'] != 'junction':
                    continue
                elif promotor['strand'] != other_promoter['strand']:
                    continue
                for position_type in ['UP_position', 'hex35_position', 'spacer_position', 'hex10_position', 'disc_position']:
                    if promotor[position_type] == other_promoter[position_type]: # Keep the more negative dG total, drop the current promoter for ties
                        if promotor['dG_total'] > other_promoter['dG_total']: 
                            promotor['drop'] = True
                            break
                        if promotor['dG_total'] < other_promoter['dG_total']:
                            other_promoter['drop'] = True
                            break
                        else:
                            promotor['drop'] = True
                            break
            if promotor['drop'] == False:
                to_keep.append(promotor)

        return to_keep, length, average
        


    # Run raw output of calculator, chunk if necessary
    fraction_length = 2500


    if len(sequence) >= 5*fraction_length and False:  # Disabled code
        if verbosity >= 2:
            print('Large input detected. Chunking sequence into smaller pieces.')
        fractions = math.floor(len(sequence)/fraction_length)
        remainder = len(sequence) % fraction_length
        if remainder > 0:
            fractions += 1
        else:
            # remainder onto the last run
            pass
        
        
        if not callback:
            if verbosity >= 1:
                callback = _progress_callback(description = "Running Promoter Calculator", tasks = (len(sequence) + 2*100*(fractions-1))*2)
            else:
                callback = lambda: None
        

        run_number = 0
        
        sample_hits = []
        sample_averages = []

        for i in range(fractions):
            if i == fractions-1:
                end_pos = len(sequence)
            else:
                end_pos = (i+1)*fraction_length+100
            if i == 0:
                start_pos = 0
            else:
                start_pos = i*fraction_length-100
            
            target_sequence = sequence[start_pos:end_pos]

            result, num_hits, average = _run_calculator(target_sequence, start_pos, end_pos, i, callback=callback)
            sample_hits.append(num_hits)
            sample_averages.append(average)
            output.extend(result)
            run_number += len(target_sequence)

        # Calculate weighted average
        num_hits = sum(sample_hits)
        average = sum([x*y for x,y in zip(sample_hits, sample_averages)]) / num_hits

    else:
        if not callback:
            if verbosity >= 1:
                callback = _progress_callback(description = "Running Promoter Calculator", tasks = len(sequence)*2)
            else:
                callback = lambda: None

        result, num_hits, average = _run_calculator(sequence, 0, len(sequence), 0, callback=callback)
        output.extend(result)

    # If the DNA is circular, examine the junction
    junction_results = {}
    if circular:
        if verbosity >= 2:
            print("Circular DNA detected. Examining junction.")
        # Get the junction
        junction = sequence[-100:] + sequence[:100]
        result, num_hits, average = _run_calculator(junction, len(sequence)-100, len(sequence)+100, "junction")
        # Filter out promotors that don't cross the junction
        for result in result:
            if result['strand'] == '+':
                bounds = [result['UP_position'][0], result['disc_position'][1]]
            if result['strand'] == '-':
                bounds = [result['disc_position'][0], result['UP_position'][1]]
            if bounds[0] <= 100 and bounds[1] >= 101:
                # Subtract 100 from the location to make the position appropriate relative to the rest of the sequence
                result['UP_position'][0] -= 100
                result['UP_position'][1] -= 100
                result['hex35_position'][0] -= 100
                result['hex35_position'][1] -= 100
                result['spacer_position'][0] -= 100
                result['spacer_position'][1] -= 100
                result['hex10_position'][0] -= 100
                result['hex10_position'][1] -= 100
                result['disc_position'][0] -= 100
                result['disc_position'][1] -= 100
                result['TSS'] -= 100
                if result['TSS'] < 0:
                    result['TSS'] += len(sequence)
                junction_results[f"{result['TSS']}{result['strand']}"] = result
    
        # Replace TSS with the actual position in the sequence
        for i in range(len(output)):
            position = output[i]['TSS']
            strand = output[i]['strand']
            if f"{position}{strand}" not in junction_results.keys():
                continue  # Ignore positions that are not hit by circularization
            # Take the more negative dG_total
            if output[i]['dG_total'] > junction_results[f"{position}{strand}"]['dG_total']:
                output[i] = junction_results[f"{position}{strand}"]

    # Find the best results
    output.sort(key=lambda x: x['Tx_rate'])
    output = output[::-1]
    return output
