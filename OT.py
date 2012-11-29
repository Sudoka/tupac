#!/usr/bin/python -tt
import re
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def main():
  seq = '''MAPVADPKSKRYDRQIRIWGTHGQQRLESCCICLLNCGPTGSETLKNLVLGGIASFTIVDGGKVEARDLGNNFLVSASNLGEPRAKVVTELLQELNESVS
GSYVEEVPEVIIADNPAFFNGFDLVIATQLREQDAVVLDGICRASGRARLLLVRSYGLVGYLRASLPEHRIVESKPDSQLDDLRLNAPWPELSSFAASFA
LEQLDEVAHAHVPYVVLLLQAAARWRAGHGGGLPATSADKAAFKAAVGGMRRTADGVPLPSENFDEALKAAFHVWTPYAIPSEVRALLADDAASLSGAGG
GAGGAGGGGGGGGGLGPGSDDFWVLVAALRAFVEGEGGGCLPLEGSIPDMHATTDDAALYVLLRAADRFFAQTGRYPGDVTPGDPSDDIPLLRQAAQQAG
RGVCVVAVCVRACAGMCVHVLIDTGVVPGGSGGSASGVGFNPRKSPDSSSGSGAAAAAAVSEDLIAEFCRFGAAELHVVAAFMGGVAAQEAIKLVTRQFV
PLAGSLVYNAMSATTTVLEL'''.replace('\n', '')

  blast_handle = NCBIWWW.qblast('tblastn', 'nr', seq, entrez_query='scenedesmus dimorphus')
  blast_handle.seek(0)
  records = NCBIXML.parse(blast_handle)
  for record in records:
    for alignment in record.alignments:
      print len(alignment.hsps)
      for hsp in alignment.hsps:
        print hsp.query
        print hsp.match
        print hsp.sbjct

# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()

