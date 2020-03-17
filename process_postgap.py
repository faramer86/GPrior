#!/usr/bin/env python3.6

from modules import *
import GPriorPostgapImportProcessingFunctions as import_fun
import MLAutomationFunctions as proc_fun
import AdditionalFeatures as proc_fts
warnings.filterwarnings("ignore", category=FutureWarning)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='./gprior.py',
                            usage='%(prog)s [options]',
                            description="""
                            -//-
                            """)
    parser.add_argument('-i', '--input', required=True,
                        help='Path to file/folder with postgap output')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to output')
    parser.add_argument('--merge', action='store_true', required=False,
                        help='-//-')
    parser.add_argument('--aggregate', action='store_true', required=False,
                        help='-//-')
    
    args = parser.parse_args()

    if args.merge:
        try:
            input_file = import_fun.import_postgap_file(args.input, COL_DTYPES)
        except ZeroDivisionError:
            print('\nThere are no .tsv files in repository!')
            sys.exit()
        input_file.to_csv(Path(args.output), sep='\t')
        sys.exit("Done!")

    if args.aggregate:
        try:
            input_file = import_fun.import_merged_files(args.input)
            input_file.to_csv(f'{args.output}', sep='\t')
        except ZeroDivisionError:
            print('\nThere are no .tsv files in repository!')
            sys.exit()
        sys.exit("Done!")
