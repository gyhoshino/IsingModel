from fn_ising import unpickle_to_csv
from sys import argv, exit
from os import path

if len(argv) < 2:
    exit('fatal: requires an input directory containing the pickle/ directory.')

if not path.isdir(argv[1]):
    exit(f'fatal: input "{argv[1]}" is not a directory.')

if not path.isdir(f'{argv[1]}/pickles'):
    exit(f'fatal: directory "pickles" is not present in {argv[1]}/')

unpickle_to_csv(argv[1])
