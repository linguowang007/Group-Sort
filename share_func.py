import itertools
import os
import gzip
import numpy as np
from collections import Counter
from functools import lru_cache


def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    in_num1, in_num2 = itertools.tee(iterable)
    next(in_num2, None)
    return zip(in_num1, in_num2)


def de_duplicates(input_list, key=None):
    # de_duplicates('1 2 4 2 5 4 6'.split()) --> 1 2 4 5 6
    seen = set()
    dedup_list = []
    for item in input_list:
        key_value = item if key is None else key(item)
        if key_value not in seen:
            seen.add(key_value)
            dedup_list.append(item)
    return dedup_list


def get_chr_order(*, sam_file, chunk_pos_list):  # chunk-file, start-pos
    os.system(f'samtools view -H {sam_file} > {sam_file}-for-sort.sam')
    with open(f'{sam_file}-sort.sh', 'wt') as h1:
        for chunk, pos1, in chunk_pos_list:
            print(f'bgzip -b {pos1} {chunk} | head -1 >> {sam_file}-for-sort.sam', file=h1)
    os.system(rf'''sh {sam_file}-sort.sh && rm {sam_file}-sort.sh
            samtools sort -O SAM {sam_file}-for-sort.sam | \
            samtools view | \
            perl -nale 'print $F[2]' > {sam_file}-for-sort.sam.chr''')
    with open(f'{sam_file}-for-sort.sam.chr') as h1:
        result = {i.strip(): num for num, i in enumerate(h1)}
    os.system(f'rm {sam_file}-for-sort.sam {sam_file}-for-sort.sam.chr')
    
    return result


def get_uniq_pos_acc_len(lines):
    for pos, pos_group in itertools.groupby(lines, lambda x: x[1]):
        total_length = sum(i[-1] for i in pos_group)
        yield pos, total_length


def get_keep_bin_unmap(filename):
    total = Counter()
    acc = 0
    keep_bin = []
    with gzip.open(filename, 'rt') as h1:
        m1 = ((i.split('\t', maxsplit=4), len(i)) for i in h1)
        m1 = ((i[2], int(i[3]), b) for i, b, in m1)
        for chrom, chrom_group in itertools.groupby(m1, lambda x: x[0]):
            d1 = np.fromiter(itertools.chain.from_iterable(get_uniq_pos_acc_len(chrom_group)), dtype=np.uint32)
            len1 = d1[1::2]
            total[chrom] = int(len1.sum())
            if chrom != '*':
                d1.tofile(f'{filename}-{chrom}-sort-acc.bin')
                keep_bin.append([filename, chrom, acc, f'{filename}-{chrom}-sort-acc.bin'])
                acc += int(len1.sum())
    unmapped = None
    if total['*']:
        off_base = total['*']
        tmp1 = sum(val for key, val in total.items() if key != '*')
        unmapped = f'{filename}-unmapped.gz'
        os.system(f'bgzip -b {tmp1} -s {off_base} {filename} | bgzip > {unmapped}')
    return keep_bin, unmapped


@lru_cache(maxsize=1)
def get_pos_len(bin_file):
    d1 = np.fromfile(bin_file, dtype=np.uint32)
    return d1
