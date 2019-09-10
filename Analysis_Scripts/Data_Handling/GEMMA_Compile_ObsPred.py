#!/usr/bin/env python
"""Simple script to join the observed and predicted phenotypic values from the
GEMMA cross-validation. Takes two arguments:
    1) Observed phenotypes .FAM
    2) Directory with predicted phenotypes
"""

import sys
import gzip
import os


def get_preds(pd):
    """Crawl the predictions directory and reutrn a path to each one."""
    fp = os.path.abspath(os.path.expanduser(pd))
    p = [
        os.path.join(fp, x)
        for x
        in os.listdir(fp)
        if x.endswith('prdt.txt.gz')
        ]
    return p



def main(obs, pred):
    """Main function."""
    predictions = get_preds(pred)
    # Read the observed phenotypes into the 
    phen = []
    with open(obs, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            phen.append(tmp[-1])
    print('Observed Predicted')
    for p in predictions:
        prd = []
        with gzip.open(p, 'rt') as f:
            for line in f:
                prd.append(line.strip())
        for ob, pr in zip(phen, prd):
            print(ob + ' ' + pr)
    return


main(sys.argv[1], sys.argv[2])
