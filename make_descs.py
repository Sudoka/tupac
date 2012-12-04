#!/usr/bin/python -tt
import re
from sys import argv as args
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import cPickle as pickle
import time

def main():
  if '-h' in args or '--h' in args or '--help' in args:
    print 'Usage is as follows: ./makefile -i infile -o outfile -s stepsize -t runtime_in_seconds'
    print 'Good default values are -s 500 and -t 120'
    sys.exit()

  vals = [args[args.index(flag) + 1] for flag in ['-i', '-o', '-s', '-t']]
  infile = vals[0]
  outfile = vals[1]
  step = int(vals[2])
  runtime = int(vals[3]) 
 
  seq_dict = pickle.load(open(infile, 'r'))

  #offset for chunks of sequences to blast
  o = 0
  start = time.time()
  seq_dict_keys = seq_dict.keys(); n = len(seq_dict_keys)
  descs = []

  #logic that queries as many queries possible in time given, and if we run out of queries before allotted time, end querying
  while time.time()-start < runtime and o+step < n:
    sample_seqs = [(k, seq_dict[k]) for k in seq_dict_keys[o:o+step]]
    descs += get_descriptions(sample_seqs)
    o += step
    print 'Time Elapsed:', int(time.time()-start), 'seconds'
  pickle.dump(descs, open(outfile, 'w'))


# given list of tuples of form [('gi|xxx|xxxx..', 'AA_SEQUENCE'), ... ]
# and performs queries based on this.
def get_descriptions(sample_seqs):
  query_string = ''
  for seq in sample_seqs:
    query_string += seq[0]+'\n'+seq[1]+'\n'

  blast_handle = NCBIWWW.qblast('tblastn', 'nr', query_string, entrez_query='scenedesmus dimorphus')
  blast_handle.seek(0)
  records = NCBIXML.parse(blast_handle)
  descs = []
  i = 0
  for record in records:
    if len(record.alignments) > 0:
      for align in record.alignments:
        desc = [sample_seqs[i][0], align.hit_id]

        frames = [hsp.frame[1] for hsp in align.hsps]
        if valid_align(frames):
          desc.append(plus_or_minus(frames[0]))
        else:
          desc.append('/')

        query_coverage = float(sum([len(hsp.sbjct) for hsp in align.hsps])) / len(sample_seqs[i][1])
        desc.append(query_coverage)

        #list of tuples of form ( (query_start, query_end), (sbjct_start, sbjct_end), (query, match, sbjct, frame) ) sorted by query_start
        hsp_info = sorted([((hsp.query_start, hsp.query_end), (hsp.sbjct_start, hsp.sbjct_end), (hsp.query, hsp.match, hsp.sbjct, hsp.frame[1])) for hsp in align.hsps], key= lambda t: t[0][0])
        desc.append(hsp_info)

    else:
      desc = [sample_seqs[i][0], ' ', ' ', 0.0, [], []]

    descs.append(desc)
    i += 1
  return descs

def plus_or_minus(n):
  if n < 0:
    return '-'
  return '+'

#given list of frame values, returns true if all are < 0 or all are > 0
#this ensures the alignment is either sense or antisense derived
def valid_align(x):
  return all(item > 0 for item in x) or all(item < 0 for item in x)

# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()

