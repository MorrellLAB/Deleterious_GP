#!/usr/bin/env python
"""Associate the PolyPhen2 predictions with the output from SNP_Effect_Predictor
and print to stdout. Relies on the transcript ID, aa1, aa2, and original
position data to be consistent."""

import sys

effect_table = sys.argv[1]
pph_table = sys.argv[2]

#   Read the PPH2 table, and store the information in a dictionary.
pph_data = {}
with open(pph_table, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            #   The PPH2 fields are strangely formatted - padded with
            #   whitespace and separated by tabs. We strip the newline, then
            #   split on tabs, then strip the whitespace off each field
            tmp = [s.strip() for s in line.strip().split('\t')]
            #   Then, we build a key out of the transcript ID, position, and
            #   the amino acid states. These are fields 1, 2, 3, and 4,
            #   respectively. The prediction is field 15. Not a very clean key,
            #   but it's the only way to uniquely key a SNP with its prediction.
            k = (tmp[0], tmp[1], tmp[2], tmp[3])
            v = tmp[14]
            pph_data[k] = v


#   Then iterate through the effects table and print out the PPH prediction at
#   the end of the table.
with open(effect_table, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            #   Print the first line with a header for PPH2 prediction
            print line.strip() + '\tPPH2'
        else:
            tmp = line.strip().split()
            #   Rebuild the key for looking up PPH2 predictions. Be aware that
            #   the transcript name has to hvae full stops replaced with
            #   underscores.
            # tid = tmp[4].replace('.', '_')
            #   We don't need to do the replacement anymore!
            tid = tmp[4]
            k = (tid, tmp[11], tmp[8], tmp[9])
            if k in pph_data:
                towrite = line.strip() + '\t' + pph_data[k]
            else:
                towrite = line.strip() + '\t-'
            print towrite
