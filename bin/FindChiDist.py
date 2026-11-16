"""Original breakpoint clustering and calibrated expression; no per-read index.

The greedy local/global merge and representative selection intentionally retain
the reference implementation's order and thresholds. Retained-side directions
are counted inline without SAM rescans or per-read tracking. This is an
internal subprocess of scFusion.py.
"""
import csv
import hashlib
import os
import sys
from functools import lru_cache
import numpy

DIRECTION_PAIRS = (('+', '+'), ('+', '-'), ('-', '+'), ('-', '-'))


def ParsePosition(key):
    positions = key.split(',')
    return int(positions[0]), int(positions[1])


ParsePosition = lru_cache(maxsize=8192)(ParsePosition)


def LoadDirections(path):
    """Geometry keys only; absent legacy sidecars yield unknown directions."""
    if not os.path.isfile(path):
        return {}
    directions = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (row['Gene1'], row['Gene2'], row['EvidenceKey'])
            direction = (row['direction1'], row['direction2'])
            if key in directions and directions[key] != direction:
                directions[key] = ('?', '?')
            else:
                directions[key] = direction
    return directions


def DirectionVector(direction):
    counts = [0, 0, 0, 0, 0]
    index = DIRECTION_PAIRS.index(direction) if direction in DIRECTION_PAIRS else 4
    counts[index] = 1
    return counts


def MergeDirections(target, source):
    for index in range(5):
        target[index] += source[index]


def SwapDirections(counts):
    return [counts[0], counts[2], counts[1], counts[3], counts[4]]


def NewCluster(positions, average, cell, directions):
    # Keep the legacy indices: member coordinates=0, average=1, cells=2,
    # representative coordinates=3. Five tiny direction counters are index 4.
    return [positions, average, [cell], None, list(directions)]


def RepresentativeDirection(counts):
    maximum = max(counts[:4])
    if maximum == 0:
        return ('?', '?', 'UNKNOWN')
    winners = [index for index, count in enumerate(counts[:4]) if count == maximum]
    if len(winners) != 1:
        return ('?', '?', 'TIE')
    return DIRECTION_PAIRS[winners[0]] + ('RESOLVED',)


def AddCount(genename, q=1.0):
    genename = genename.split(';')
    for gene in genename:
        gene = gene.split('.')[0]
        if gene in ExprDic:
            if ExprDic[gene][1] == i:
                ExprDic[gene][0][-1] += q
            else:
                for k1 in range(ExprDic[gene][1]+1, i):
                    ExprDic[gene][0].append(0)
                ExprDic[gene][0].append(q)
                ExprDic[gene][1] = i
        else:
            ExprDic[gene] = [[], i]
            for k1 in range(start, i):
                ExprDic[gene][0].append(0)
            ExprDic[gene][0].append(q)


def LoadChimericExpressionCount(cellindex):
    """Load mappability-passed, pre-hotspot chimeric expression once."""
    path = os.path.join(ChimericExpressionDir, str(cellindex) + '.gene_counts.tsv')
    if not os.path.isfile(path):
        raise IOError('Missing chimeric expression count file: ' + path)
    with open(path) as countfile:
        for line in countfile:
            if line.startswith('Gene\t'):
                continue
            info = line.rstrip().split('\t')
            if len(info) < 2 or info[0] == '':
                continue
            AddCount(info[0], float(info[1]))



filedir = sys.argv[1]
start, last = int(sys.argv[2]), int(sys.argv[3])
ExprDir = sys.argv[4]
prefix = sys.argv[6] if len(sys.argv) > 6 and sys.argv[6] != '.' else ''
ChimericExpressionDir = os.environ['SCFUSION_CHIMERIC_EXPR_DIR']
FusionMatrix, EncompassingMatrix, FusionPos, ExprDic = {}, {}, {}, {}
PairSpanning = {}
HomologyDic = {}
good_dist = 1
with open(sys.argv[5]) as homology:
    for line in homology:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\r\n').split('\t')
        HomologyDic['\t'.join(fields[index] for index in (3, 4, 5, 6))] = fields[12:14]

for i in range(start, last + 1):
    infile = open(os.path.join(filedir, str(i) + '_FusionSupport.txt'))
    directions = LoadDirections(os.path.join(filedir, str(i) + '_FusionSupport.directions.tsv'))
    with open(os.path.join(filedir, str(i) + '_FusionSupport.read_length.txt')) as metadata:
        readlength = int(metadata.read().strip())
    with open(os.path.join(ExprDir, str(i) + '.rpkm.txt')) as exprfile:
        lastgene, totalrc = '', 0
        for line in exprfile:
            fields = line.split('\t')
            gene, count = fields[0], int(fields[1])
            if gene == lastgene:
                totalrc += count
            else:
                if lastgene:
                    AddCount(lastgene, totalrc)
                lastgene, totalrc = gene, count
        if lastgene:
            AddCount(lastgene, totalrc)
    LoadChimericExpressionCount(i)
    for line in infile:
        info = line.split('\t')
        gene1 = info[0]
        gene2 = info[1]
        chromo1 = info[4]
        chromo2 = info[5]
        encompass = int(info[2])
        pair = tuple(sorted((gene1, gene2)))
        PairSpanning.setdefault(pair, {})[i] = max(PairSpanning.get(pair, {}).get(i, 0), encompass)
        if gene1 + '\t' + gene2 in EncompassingMatrix:
            if encompass > 0:
                EncompassingMatrix[gene1 + '\t' + gene2].append(encompass)
        elif gene2 + '\t' + gene1 in EncompassingMatrix:
            if encompass > 0:
                EncompassingMatrix[gene2 + '\t' + gene1].append(encompass)
        else:
            if encompass > 0:
                EncompassingMatrix[gene1 + '\t' + gene2] = [encompass]
            else:
                EncompassingMatrix[gene1 + '\t' + gene2] = []
        splitcount = int(info[3])
        if splitcount == 0:
            continue
        splitreadinfo = info[6].split(';')
        # Deduplicate deterministically: the following greedy clustering must
        # not depend on Python's process-randomized set iteration order.
        splitreadinfo = sorted(set(splitreadinfo))
        subscore = []
        subgrp = []
        subdirections = []
        for item in splitreadinfo: #扫描每个split-read的情况
            if len(item) > 3:
                iteminfo = item.split('+')
                if not 6 <= int(iteminfo[1]) <= readlength - 6:
                    continue
                pos = iteminfo[0].split(',')
                position_key = pos[0] + ',' + pos[1]
                pos0, pos1 = ParsePosition(position_key)
                ssscore = 1
                direction_counts = DirectionVector(directions.get((gene1, gene2, item), ('?', '?')))
                # subgrp是一个LIST，存单个融合每个聚类每个断点的count，配合subscore使用
                if len(subgrp) == 0:
                    subgrp = [{position_key: 1}]
                    subscore.append(ssscore)
                    subdirections.append(direction_counts)
                else:
                    fin = False
                    # 两层循环是因为同时比聚类和断点两层
                    for j in range(len(subgrp)): 
                        for key in subgrp[j]:
                            grpos0, grpos1 = ParsePosition(key)
                            # good_dist = 1
                            if abs(pos0 - grpos0) <= 3 * good_dist and abs(pos1 - grpos1) <= 3 * good_dist:
                                fin = True
                                if position_key in subgrp[j]:
                                    subgrp[j][position_key] += 1
                                else: # 这个是进入到这个聚类了，但是以前没有这个断点，加一个断点
                                    subgrp[j][position_key] = 1
                                subscore[j] += ssscore
                                MergeDirections(subdirections[j], direction_counts)
                                break
                        if fin:
                            break
                    if not fin: # 不再任何聚类的任何断点，则加一个聚类
                        subgrp.append({position_key: 1})
                        subscore.append(ssscore)
                        subdirections.append(direction_counts)
        if len(subgrp) == 0:
            continue
        rawscore = []
        for j in subscore:
            if j > 0.5:
                rawscore.append(j)
            else:
                rawscore.append(0) # 这不会被执行到
        AvePos = [] # 存储该融合每个聚类的平均左右坐标和总split-read数
        for mdict in subgrp:  # calculate average pos for each cluster
            left = 0
            right = 0
            cc = 0
            for key in mdict:
                pos0, pos1 = ParsePosition(key)
                left += pos0 * mdict[key]
                right += pos1 * mdict[key]
                cc += mdict[key]
            left = int(left / cc)
            right = int(right / cc)
            AvePos.append([left, right, cc])
        if gene1 + '\t' + gene2 in FusionMatrix:
            added = list(numpy.zeros(len(FusionMatrix[gene1 + '\t' + gene2])))
            for j in range(len(subgrp)): # 将该细胞的每个聚类合并到总全局聚类中
                # 记录这个细胞这个融合的这个聚类，距离已有的聚类（之前细胞检测出的）
                # 的距离列表
                dist = []
                for k in FusionPos[gene1 + '\t' + gene2][1]:
                    d1 = abs(k[1][0] - AvePos[j][0])
                    d2 = abs(k[1][1] - AvePos[j][1])
                    if d1 <= good_dist and d2 <= good_dist:
                        dist.append(d1 + d2)
                    else:
                        dist.append(10 * good_dist)
                if min(dist) >= 3 * good_dist:
                    # 如果距离最近的已有聚类都很远，则创建一个新聚类，这里的评分仅仅是split的支持数量
                    FusionPos[gene1 + '\t' + gene2][1].append(NewCluster(subgrp[j], AvePos[j], i, subdirections[j]))
                    FusionMatrix[gene1 + '\t' + gene2].append([rawscore[j]])
                    # This newly created cluster already contains a
                    # contribution from the current cell.  A later local
                    # subgroup from the same cell that maps to this cluster
                    # must be summed into the last entry, not appended as a
                    # fictitious second supporting cell.
                    added.append(1)
                else: # 和现有聚类合并
                    smallindex = dist.index(min(dist)) # 
                    if added[smallindex] == 0:
                        FusionMatrix[gene1 + '\t' + gene2][smallindex].append(rawscore[j])
                        added[smallindex] = 1
                    else:
                        FusionMatrix[gene1 + '\t' + gene2][smallindex][-1] += (rawscore[j])
                    old = FusionPos[gene1 + '\t' + gene2][1][smallindex][1]
                    newcc = AvePos[j][2] + old[2]
                    newleft = int((old[2] * old[0] + AvePos[j][2] * AvePos[j][0]) / newcc)
                    newright = int((old[2] * old[1] + AvePos[j][2] * AvePos[j][1]) / newcc)
                    FusionPos[gene1 + '\t' + gene2][1][smallindex][1] = [newleft, newright, newcc]
                    FusionPos[gene1 + '\t' + gene2][1][smallindex][2].append(i)
                    MergeDirections(FusionPos[gene1 + '\t' + gene2][1][smallindex][4], subdirections[j])
            gene = gene1 + '\t' + gene2
        elif gene2 + '\t' + gene1 in FusionMatrix:
            added = list(numpy.zeros(len(FusionMatrix[gene2 + '\t' + gene1])))
            for j in range(len(subgrp)):
                dist = []
                for k in FusionPos[gene2 + '\t' + gene1][1]:
                    d1 = abs(k[1][0] - AvePos[j][1])
                    d2 = abs(k[1][1] - AvePos[j][0])
                    if d1 <= good_dist and d2 <= good_dist:
                        dist.append(d1 + d2)
                    else:
                        dist.append(10 * good_dist)
                if min(dist) >= 3 * good_dist:
                    FusionPos[gene2 + '\t' + gene1][1].append(
                        NewCluster(subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], i,
                                   SwapDirections(subdirections[j])))
                    FusionMatrix[gene2 + '\t' + gene1].append([rawscore[j]])
                    # Same current-cell accounting rule for reversed gene
                    # order (see the forward-order branch above).
                    added.append(1)
                else:
                    smallindex = dist.index(min(dist))
                    if added[smallindex] == 0:
                        FusionMatrix[gene2 + '\t' + gene1][smallindex].append(rawscore[j])
                        added[smallindex] = 1
                    else:
                        FusionMatrix[gene2 + '\t' + gene1][smallindex][-1] += (rawscore[j])
                    old = FusionPos[gene2 + '\t' + gene1][1][smallindex][1]
                    newcc = AvePos[j][2] + old[2]
                    newleft = int((old[2] * old[0] + AvePos[j][2] * AvePos[j][1]) / newcc)
                    newright = int((old[2] * old[1] + AvePos[j][2] * AvePos[j][0]) / newcc)
                    FusionPos[gene2 + '\t' + gene1][1][smallindex][1] = [newleft, newright, newcc]
                    FusionPos[gene2 + '\t' + gene1][1][smallindex][2].append(i) # 这是关键
                    MergeDirections(FusionPos[gene2 + '\t' + gene1][1][smallindex][4],
                                    SwapDirections(subdirections[j]))
            gene = gene2 + '\t' + gene1
        else:
            FusionMatrix[gene1 + '\t' + gene2] = []
            for item in rawscore:
                FusionMatrix[gene1 + '\t' + gene2].append([item])
            FusionPos[gene1 + '\t' + gene2] = [[chromo1, chromo2], []]
            for k in range(len(subgrp)):
                FusionPos[gene1 + '\t' + gene2][1].append(NewCluster(subgrp[k], AvePos[k], i, subdirections[k]))
            gene = gene1 + '\t' + gene2


    infile.close()

for genename in ExprDic: # 添加那些末尾不表达该基因的细胞表达量
    if ExprDic[genename][1] < last:
        for k in range(ExprDic[genename][1]+1, last+1):
            ExprDic[genename][0].append(0)
for gene in FusionMatrix:
    temp = FusionPos[gene][1]
    for i in range(len(temp)):
        keylist = sorted(temp[i][0], key=temp[i][0].__getitem__, reverse=True)
        possum = 0
        countlist = []
        for key in keylist:
            possum += temp[i][0][key]
            countlist.append(int(temp[i][0][key]))
        if possum > 0 and max(countlist) / possum >= 0.4:
            for key in temp[i][0]:
                if key == keylist[0]:
                    pos0, pos1 = ParsePosition(keylist[0])
                    if abs(pos0 - temp[i][1][0]) + abs(pos1 - temp[i][1][1]) < abs(pos1 - temp[i][1][0]) + abs(
                            pos0 - temp[i][1][1]):
                        FusionPos[gene][1][i][3] = [pos0, pos1]
                    else:
                        FusionPos[gene][1][i][3] = [pos1, pos0]
                    break
        else:
            FusionPos[gene][1][i][3] = [-1, -1]

templines = []
for gene in FusionMatrix:
    for i in range(len(FusionPos[gene][1])):
        if FusionPos[gene][1][i][3][0] > -0.5:
            genepart = gene.split('\t')
            usecell = FusionPos[gene][1][i][2]
            # Cells are processed in ascending order and FusionMatrix keeps
            # counts in that order. Keep the cell list deterministic and
            # positionally aligned with SplitReadCountsByCell.
            usecell = sorted(set(usecell))
            split_counts_by_cell = FusionMatrix[gene][i]
            if len(usecell) != len(split_counts_by_cell):
                raise RuntimeError(
                    'Fusion cluster has inconsistent supporting cells/counts: ' +
                    gene.replace('\t', '--') + ' cluster=' + str(i) +
                    ' cells=' + str(sorted(usecell)) +
                    ' counts=' + str(split_counts_by_cell))
            if int(sum(split_counts_by_cell)) != int(FusionPos[gene][1][i][1][2]):
                raise RuntimeError(
                    'Fusion cluster has inconsistent total split count: ' +
                    gene.replace('\t', '--') + ' cluster=' + str(i) +
                    ' per_cell=' + str(split_counts_by_cell) +
                    ' cluster_total=' + str(FusionPos[gene][1][i][1][2]))
            try:
                aaaaa = []
                thisgenename = genepart[0].split('.')[0]
                for cell in usecell:
                    aaaaa.append(ExprDic[thisgenename][0][cell-start])
            except:
                aaaaa = list(numpy.zeros(last-start+1, int))
                continue
            try:
                bbbbb = []
                thisgenename = genepart[1].split('.')[0]
                for cell in usecell:
                    bbbbb.append(ExprDic[thisgenename][0][cell-start])
            except:
                bbbbb = list(numpy.zeros(last-start+1, int))
                continue
            if FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1] + '\t' + str(FusionPos[gene][1][i][3][0]) + '\t' \
                    + str(FusionPos[gene][1][i][3][1]) in HomologyDic:
                gccontent = HomologyDic[FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1]
                                                                + '\t' + str(FusionPos[gene][1][i][3][0]) + '\t' + str(
                    FusionPos[gene][1][i][3][1])]
            elif FusionPos[gene][0][1] + '\t' + FusionPos[gene][0][0] + '\t' + str(FusionPos[gene][1][i][3][1]) + '\t' \
                    + str(FusionPos[gene][1][i][3][0]) in HomologyDic:
                gccontent = HomologyDic[FusionPos[gene][0][1] + '\t' + FusionPos[gene][0][0]
                                                                + '\t' + str(FusionPos[gene][1][i][3][1]) + '\t' + str(
                    FusionPos[gene][1][i][3][0])]
            else:
                continue
            engene = gene
            if gene not in EncompassingMatrix:
                engene = genepart[1] + '\t' + genepart[0]
            pair_counts = PairSpanning[tuple(sorted(genepart))]
            pair_cells = sorted(cell for cell, count in pair_counts.items() if count > 0)
            fields = [
                genepart[0], genepart[1], len(split_counts_by_cell), str(split_counts_by_cell),
                len(EncompassingMatrix[engene]), str(EncompassingMatrix[engene]),
                FusionPos[gene][0][0], FusionPos[gene][0][1],
                FusionPos[gene][1][i][3][0], FusionPos[gene][1][i][3][1], i,
                gccontent[0], gccontent[1], str(aaaaa), str(bbbbb),
            ]
            identity = tuple(str(fields[index]) for index in (0, 1, 6, 7, 8, 9, 10))
            digest = hashlib.sha1('\t'.join(identity).encode('utf-8')).hexdigest()[:12]
            cid = '{}--{}|{}:{}|{}:{}|cluster{}|{}'.format(
                fields[0], fields[1], fields[6], fields[8], fields[7], fields[9], fields[10], digest)
            fields.extend([cid, str(usecell), str(pair_cells),
                           str([pair_counts[cell] for cell in pair_cells])])
            direction_counts = FusionPos[gene][1][i][4]
            direction1, direction2, direction_status = RepresentativeDirection(direction_counts)
            fields.extend([direction1, direction2, *direction_counts, direction_status])
            templines.append(fields)

HEADERS = (
    'Gene1', 'Gene2', 'CellSupport', 'SplitReadCountsByCell',
    'SpanningSupportingCellCount', 'SpanningReadCountsByCell',
    'Chrom1', 'Chrom2', 'Breakpoint1', 'Breakpoint2', 'ClusterIndex',
    'GC1', 'GC2', 'Gene1ExpressionBySplitCells', 'Gene2ExpressionBySplitCells',
    'CandidateID', 'SupportingCells', 'PairSpanningCells', 'PairSpanningCounts',
    'direction1', 'direction2', 'DirectionCountPP', 'DirectionCountPM',
    'DirectionCountMP', 'DirectionCountMM', 'DirectionAmbiguousCount', 'DirectionStatus')
print('\t'.join(HEADERS))
for fields in sorted(templines, key=lambda row: len(row[3]), reverse=True):
    print('\t'.join(map(str, fields)))
