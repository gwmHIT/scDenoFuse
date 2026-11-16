from __future__ import print_function
from __future__ import division
import sys
from bisect import bisect_right
from heapq import heappop, heappush
from pathlib import Path
import pysam
sys.path.insert(0, str(Path(__file__).resolve().parent))
from scfusion.orientation import GeneStrands


# ***** readme *****
# This code wants to verify that fake fusion is more likely to be homo


class GeneSpanIndex:
    """Disjoint lookup intervals retaining the original first-gene preference."""

    def __init__(self, genes):
        intervals = sorted((bounds[0][0], order, bounds[0][1], gene)
                           for order, (gene, bounds) in enumerate(genes.items()))
        boundaries = sorted({coordinate for start, _, end, _ in intervals
                             for coordinate in (start, end + 1)})
        active = []
        cursor = 0
        segments = []
        for start, following in zip(boundaries, boundaries[1:]):
            while cursor < len(intervals) and intervals[cursor][0] <= start:
                _, order, end, gene = intervals[cursor]
                heappush(active, (order, end, gene))
                cursor += 1
            while active and active[0][1] < start:
                heappop(active)
            if not active:
                continue
            gene = active[0][2]
            end = following - 1
            if segments and segments[-1][2] == gene and segments[-1][1] + 1 == start:
                segments[-1] = (segments[-1][0], end, gene)
            else:
                segments.append((start, end, gene))
        self.starts = [segment[0] for segment in segments]
        self.segments = segments

    def first_gene(self, position):
        index = bisect_right(self.starts, position) - 1
        if index >= 0 and position <= self.segments[index][1]:
            return self.segments[index][2]
        return None


def GCcontent(str1, str2):
    l1 = len(str1)
    l2 = len(str2)
    count = 0
    for i in range(l1):
        if str1[i].upper() == 'G' or str1[i].upper() == 'C':
            count += 1
    res1 = count / l1
    count = 0
    for i in range(l2):
        if str2[i].upper() == 'G' or str2[i].upper() == 'C':
            count += 1
    res2 = count / l2
    return [res1, res2]


def getTCGAstring(chromo, start, end):
    """Original one-based inclusive sequence lookup, without whole-genome RAM."""
    contig = chromo if chromo in reference.references else chromo.removeprefix('chr')
    return reference.fetch(contig, start - 1, end)


def onethread(halfwindowwidth, fusionfile):
    geneset = {}
    strands = GeneStrands() if len(sys.argv) > 4 else None
    for line in refannotfile:
        if line.startswith('#'):
            continue
        info = line.split('\t')
        if strands is not None:
            strands.observe(info)
        if info[2] == 'exon':
            chromo = info[0]
            start = int(info[3])
            end = int(info[4])
            geneinfo = info[8].split('; ')
            genename = ''
            genetype = ''
            genetypeFound = False
            if len(geneinfo) > 4:
                for item in geneinfo:
                    if item.find('gene_name') > -1:
                        genename = item[11:-1]
                    if item.find('gene_type') > -1:
                        genetype = item[11:-1]
                        genetypeFound = True
                    if item.find('gene_biotype') > -1:
                        genetype = item[14:-1]
                        genetypeFound = True
            else:
                continue
            if genename == '':
                continue
            if genetype != 'protein_coding' and genetype != 'processed_transcript' and genetypeFound:
                continue
            if chromo not in geneset:
                geneset[chromo] = {}
            if genename not in geneset[chromo]:
                geneset[chromo][genename] = [[start, end], [[start, end]]]
            else:
                found = False
                for item in geneset[chromo][genename][1]:
                    if start <= item[1] and end >= item[0]:
                        found = True
                        break
                if found:
                    continue
                geneset[chromo][genename][1].append([start, end])
                geneset[chromo][genename][0][0] = min(geneset[chromo][genename][0][0], start)
                geneset[chromo][genename][0][1] = max(geneset[chromo][genename][0][1], end)
    refannotfile.close()
    if strands is not None:
        strands.write(sys.argv[4])
    for i in geneset:
        for j in geneset[i]:
            geneset[i][j][1].sort()

    gene_indexes = {chrom: GeneSpanIndex(genes) for chrom, genes in geneset.items()}
    for line in fusionfile:
        if line[0] == '#':
            continue
        try:
            info = line.split('\t')
            gene1 = info[0]
            gene2 = info[1]
            chr1 = info[3].replace(' ', '')
            chr2 = info[4].replace(' ', '')
            if not chr1.startswith('chr'):
                chr1 = 'chr' + chr1
            if not chr2.startswith('chr'):
                chr2 = 'chr' + chr2
            pos1 = int(info[5].replace(' ', ''))
            pos2 = int(info[6].replace(' ', ''))
            exoninterval1 = -1
            exoninterval2 = -1
            gene1 = gene_indexes[chr1].first_gene(pos1) or gene1
            gene2 = gene_indexes[chr2].first_gene(pos2) or gene2
            if gene1 not in geneset[chr1] or gene2 not in geneset[chr2]:
                pass
            else:
                for i in range(len(geneset[chr1][gene1][1])):
                    if geneset[chr1][gene1][1][i][0] <= pos1 <= geneset[chr1][gene1][1][i][1]:
                        exoninterval1 = i
                        break
                for i in range(len(geneset[chr2][gene2][1])):
                    if geneset[chr2][gene2][1][i][0] <= pos2 <= geneset[chr2][gene2][1][i][1]:
                        exoninterval2 = i
                        break
            # 获取外显子的序列
            if exoninterval2 == -1 or exoninterval1 == -1:
                res1d = getTCGAstring(chr1, pos1, pos1 + halfwindowslength)
                res1u = getTCGAstring(chr1, pos1 - halfwindowslength, pos1 - 1)
                res2d = getTCGAstring(chr2, pos2, pos2 + halfwindowslength)
                res2u = getTCGAstring(chr2, pos2 - halfwindowslength, pos2 - 1)
            else:
                res1d = ''
                res2d = ''
                res1u = ''
                res2u = ''
                last1 = 0
                last2 = 0
                if geneset[chr1][gene1][1][exoninterval1][1] - pos1 + 1 >= halfwindowslength:
                    res1d = getTCGAstring(chr1, pos1, pos1 + halfwindowslength - 1)
                else:
                    res1d = getTCGAstring(chr1, pos1, geneset[chr1][gene1][1][exoninterval1][1])
                    last1 = geneset[chr1][gene1][1][exoninterval1][1]
                if geneset[chr2][gene2][1][exoninterval2][1] - pos2 + 1 >= halfwindowslength:
                    res2d = getTCGAstring(chr2, pos2, pos2 + halfwindowslength - 1)
                else:
                    res2d = getTCGAstring(chr2, pos2, geneset[chr2][gene2][1][exoninterval2][1])
                    last2 = geneset[chr2][gene2][1][exoninterval2][1]
                fi1 = exoninterval1
                fi2 = exoninterval2
                while len(res1d) < halfwindowslength:
                    fi1 += 1
                    if fi1 >= len(geneset[chr1][gene1][1]):
                        break
                    if last1 >= geneset[chr1][gene1][1][fi1][0]:
                        break
                    if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(
                            res1d) >= halfwindowslength:
                        res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0],
                                               geneset[chr1][gene1][1][fi1][0] + halfwindowslength - len(res1d) - 1)
                    else:
                        res1d += getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0], geneset[chr1][gene1][1][fi1][1])
                        last1 = geneset[chr1][gene1][1][fi1][1]
                while len(res2d) < halfwindowslength:
                    fi2 += 1
                    if fi2 >= len(geneset[chr2][gene2][1]):
                        break
                    if last2 >= geneset[chr2][gene2][1][fi2][0]:
                        break
                    if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(
                            res2d) >= halfwindowslength:
                        res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0],
                                               geneset[chr2][gene2][1][fi2][0] + halfwindowslength - len(res2d) - 1)
                    else:
                        res2d += getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0], geneset[chr2][gene2][1][fi2][1])
                        last2 = geneset[chr2][gene2][1][fi2][1]
                if - geneset[chr1][gene1][1][exoninterval1][0] + pos1 + 1 >= halfwindowslength:
                    res1u = getTCGAstring(chr1, pos1 - halfwindowslength + 1, pos1)
                else:
                    res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][exoninterval1][0], pos1)
                    last1 = geneset[chr1][gene1][1][exoninterval1][0]
                if - geneset[chr2][gene2][1][exoninterval2][0] + pos2 + 1 >= halfwindowslength:
                    res2u = getTCGAstring(chr2, pos2 - halfwindowslength + 1, pos2)
                else:
                    res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][exoninterval2][0], pos2)
                    last2 = geneset[chr2][gene2][1][exoninterval2][0]
                fi1 = exoninterval1
                fi2 = exoninterval2
                while len(res1u) < halfwindowslength:
                    fi1 -= 1
                    if fi1 < 0:
                        break
                    if last1 <= geneset[chr1][gene1][1][fi1][1]:
                        break
                    if geneset[chr1][gene1][1][fi1][1] - geneset[chr1][gene1][1][fi1][0] + 1 + len(
                            res1u) >= halfwindowslength:
                        res1u = getTCGAstring(chr1,
                                              geneset[chr1][gene1][1][fi1][1] - halfwindowslength + len(res1u) + 1,
                                              geneset[chr1][gene1][1][fi1][1]) + res1u
                    else:
                        res1u = getTCGAstring(chr1, geneset[chr1][gene1][1][fi1][0],
                                              geneset[chr1][gene1][1][fi1][1]) + res1u
                        last1 = geneset[chr1][gene1][1][fi1][0]
                while len(res2u) < halfwindowslength:
                    fi2 -= 1
                    if fi2 < 0:
                        break
                    if last2 <= geneset[chr2][gene2][1][fi2][1]:
                        break
                    if geneset[chr2][gene2][1][fi2][1] - geneset[chr2][gene2][1][fi2][0] + 1 + len(
                            res2u) >= halfwindowslength:
                        res2u = getTCGAstring(chr2,
                                              geneset[chr2][gene2][1][fi2][1] - halfwindowslength + len(res2u) + 1,
                                              geneset[chr2][gene2][1][fi2][1]) + res2u
                    else:
                        res2u = getTCGAstring(chr2, geneset[chr2][gene2][1][fi2][0],
                                              geneset[chr2][gene2][1][fi2][1]) + res2u
                        last2 = geneset[chr2][gene2][1][fi2][0]
            string1 = res1u.lower() + res1d.upper()
            string2 = res2u.lower() + res2d.upper()
            gccontent = GCcontent(string1, string2)
            print(line.rstrip() + '\t-\t-\t0\t' + str(gccontent[0]) + '\t' + str(gccontent[1]) + '\t0\t0\t' + string1 + '\t' + string2)
        except:
            sys.stderr.write('Bad Line:' + line)


def main():
    onethread(halfwindowwidth, fusionfile)


if __name__ == "__main__":
    halfwindowwidth = 10
    halfwindowslength = 100
    fusionfile = open(sys.argv[1])
    reffilepath = sys.argv[2]
    refannotfile = open(sys.argv[3])
    reference = pysam.FastaFile(reffilepath)
    main()
    reference.close()
