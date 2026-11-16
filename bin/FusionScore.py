from __future__ import print_function
from __future__ import division
import sys
import math
import os
from functools import lru_cache


# Original broad clustering retained only to reproduce the GC calibration lookup.
# This code calculate fusion score for each cell and
# add rawscore for every cell


def ParsePosition(key):
    positions = key.split(',')
    return int(positions[0]), int(positions[1])


ParsePosition = lru_cache(maxsize=8192)(ParsePosition)


def AddCount(genenum, q=1.0):
    genenum = genenum.split('.')[0]
    if genenum in ExprDic:
        ExprDic[genenum] += q
    else:
        ExprDic[genenum] = q


def LoadChimericExpressionCount(cellindex):
    """Load the pre-hotspot chimeric contribution exactly once."""
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
start = sys.argv[2]
last = sys.argv[3]
ExprDir = sys.argv[4]
FusionMatrix = {}
FusionPos = {}
ExprDic = {}
good_dist = 20
ChimericExpressionDir = os.environ['SCFUSION_CHIMERIC_EXPR_DIR']
for i in range(int(start), int(last) + 1):
    try:
        infile = open(filedir + '/' + str(i) + '_FusionSupport.txt')
        exprfile = open(ExprDir + '/' + str(i) + '.rpkm.txt')
    except:
        continue
    lastgene = ''
    totalrc = 0
    ExprDic = {}
    for line in exprfile:
        info = line.split('\t')
        gene = info[0]
        rc = int(info[1])
        if gene == lastgene:
            totalrc += rc
            continue
        if lastgene != '':
            ExprDic[lastgene] = totalrc
        lastgene = gene
        totalrc = rc
    exprfile.close()
    ExprDic[lastgene] = totalrc  # 各个基因的表达计数
    LoadChimericExpressionCount(i)
    for line in infile:
        info = line.split('\t')
        gene1 = info[0]
        gene2 = info[1]
        chromo1 = info[4]
        chromo2 = info[5]
        encompass = int(info[2])
        splitcount = int(info[3])
        if splitcount + encompass == 0:
            continue
        splitreadinfo = info[6].split(';')
        splitreadinfo = list(set(splitreadinfo))        # remove duplication
        encmpscore = encompass
        subscore = []  # 和下面的subgrp配合起来使用的
        try:
            scale1 = ExprDic[gene1]
            scale2 = ExprDic[gene2]
        except:
            sys.stderr.write('!!!' + gene1 + '\t' + gene2)
            continue
        scale = math.sqrt(math.log(1+scale1, 2) * math.log(1+scale2, 2))
        subgrp = []
        for item in splitreadinfo:  # splitreadinfo是多个断点信息的集合，遍历单个read
            if len(item) > 3:  # item格式：132814563,62592244+18+18,83
                iteminfo = item.split('+')  # iteminfo[1]是断点
                pos = iteminfo[0].split(',')
                position_key = pos[0] + ',' + pos[1]
                pos0, pos1 = ParsePosition(position_key)
                badmap = iteminfo[2].split(',')
                editscore = max(50 - abs(50 - int(iteminfo[1])) - int(badmap[0]) - int(badmap[1]), 0) # 编辑距离？
                if editscore > 12:
                    ssscore = 12 + (editscore - 12) / 10  # 分数过大会衰减
                else:
                    ssscore = editscore
                if len(subgrp) == 0:
                    subgrp = [{position_key: 1}]  #对一个融合的多个read，逐个遍历read并分成多个聚类，评分和记录支持数。
                    subscore.append(ssscore)
                else:  # 对一个融合的各个支持read的断点进行聚类
                    fin = False
                    for j in range(len(subgrp)):# 下面是双重循环，每个聚类的每个断点都要检查
                        for key in subgrp[j]:
                            grpos0, grpos1 = ParsePosition(key)
                            if abs(pos0 - grpos0) <= 20 and abs(pos1 - grpos1) <= 20:
                                # 这行判断断点两侧在参考基因组上坐标是否和聚类内断点很近,<20bp
                                fin = True
                                # subgrp[j][pos[0] + ',' + pos[1]]会记录每个聚类的每个断点
                                if position_key in subgrp[j]:
                                    subgrp[j][position_key] += 1
                                else: 
                                    subgrp[j][position_key] = 1
                                subscore[j] += ssscore  # 匹配上聚类内任何断点，都给聚类加上分数
                                break
                        if fin:
                            break
                    if not fin: # 如果不在任何聚类范围内，则新建一个聚类。
                        subgrp.append({position_key: 1})
                        subscore.append(ssscore)
        rawscore = []
        for j in subscore: # 为每个聚类评分，叫做rawscore
            if j > 0.5 or encmpscore > 0: # 评分条件很宽松，只要不是超多错配以及没有一条encompass
                rawscore.append(pow(1 + (j + encmpscore) / (7 + scale), 1.3))
            else:
                rawscore.append(0)
        AvePos = []
        for mdict in subgrp:  # calculate average pos for each cluster
            left = 0
            right = 0
            cc = 0
            for key in mdict:
                pos0, pos1 = ParsePosition(key)
                left += pos0 * mdict[key]
                right += pos1 * mdict[key]
                cc += mdict[key]
            left = left / cc
            right = right / cc
            AvePos.append([left, right, cc])
        if gene1 + '\t' + gene2 in FusionMatrix:
            # gene1 + '\t' + gene2和gene2 + '\t' + gene1都出现是为了合并二者为一个条目
            for j in range(len(subgrp)):  # 这是BUG，多了一层循环
                if len(FusionPos[gene1 + '\t' + gene2][1]) == 0 and len(subgrp) > 0:
                    # 如果该基因对还没有记录
                    FusionMatrix[gene1 + '\t' + gene2] = rawscore
                    for k in range(len(subgrp)):
                        # FusionPos[gene1 + '\t' + gene2][1]是该基因对下的所有断点信息，而FusionPos[gene1 + '\t' + gene2][0]是两个染色体信息
                        FusionPos[gene1 + '\t' + gene2][1].append([subgrp[k], AvePos[k], [i], 1])
                        #这里的[i]是细胞编号，重要
                else:
                    # FusionPos是全局的记录，所以包含所有细胞，而上面的循环是单个细胞单个融合的计算，如果之前的细胞也检测到该融合，则走这一步
                    # FusionMatrix只记录该融合的分数
                    for subpos in AvePos:
                        dist = []  # 求每个集群和矩阵记录的所有细胞的该融合集群的距离
                        for k in FusionPos[gene1 + '\t' + gene2][1]:
                            d1 = abs(k[1][0] - subpos[0])
                            d2 = abs(k[1][1] - subpos[1])
                            if d1 <= good_dist and d2 <= good_dist:
                                dist.append(d1 + d2) # 符合20bp阈值，记录距离
                            else:
                                dist.append(10 * good_dist) # 太远了，标记一个很大的惩罚值
                        if len(dist) > 0:
                            if min(dist) >= 3 * good_dist:
                            # 如果发现当前断点离所有已有的点都很远，说明这可能是同一个基因对的另一种融合形式
                            # （比如不同的剪接位点）。代码会创建一个新的记录。
                                FusionPos[gene1 + '\t' + gene2][1].append([subgrp[j], AvePos[j], [i], 1])
                                FusionMatrix[gene1 + '\t' + gene2].append(rawscore[j])
                            else:
                                smallindex = dist.index(min(dist))
                                FusionMatrix[gene1 + '\t' + gene2][smallindex] += rawscore[j]
                                old = FusionPos[gene1 + '\t' + gene2][1][smallindex][1]
                                newcc = subpos[2] + old[2]
                                newleft = (old[2] * old[0] + subpos[2] * subpos[0]) / newcc
                                newright = (old[2] * old[1] + subpos[2] * subpos[1]) / newcc
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][1] = [newleft, newright, newcc]
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][2].append(i)
                                FusionPos[gene1 + '\t' + gene2][1][smallindex][3] += 1
            gene = gene1 + '\t' + gene2
        elif gene2 + '\t' + gene1 in FusionMatrix:
            for j in range(len(subgrp)):
                if len(FusionPos[gene2 + '\t' + gene1][1]) == 0 and len(subgrp) > 0:
                    FusionMatrix[gene2 + '\t' + gene1] = rawscore
                    for k in range(len(subgrp)):
                        FusionPos[gene2 + '\t' + gene1][1].append(
                            [subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], [i], 1])
                else:
                    for subpos in AvePos:
                        dist = []
                        for k in FusionPos[gene2 + '\t' + gene1][1]:
                            d1 = abs(k[1][0] - subpos[1])
                            d2 = abs(k[1][1] - subpos[0])
                            if d1 <= good_dist and d2 <= good_dist:
                                dist.append(d1 + d2)
                            else:
                                dist.append(10 * good_dist)
                        if min(dist) >= 3 * good_dist:
                            FusionPos[gene2 + '\t' + gene1][1].append(
                                [subgrp[j], [AvePos[j][1], AvePos[j][0], AvePos[j][2]], [i], 1])
                            FusionMatrix[gene2 + '\t' + gene1].append(rawscore[j])
                        else: # 合并进入这个聚类
                            smallindex = dist.index(min(dist))
                            FusionMatrix[gene2 + '\t' + gene1][smallindex] += rawscore[j]
                            old = FusionPos[gene2 + '\t' + gene1][1][smallindex][1]
                            newcc = subpos[2] + old[2]
                            newleft = (old[2] * old[0] + subpos[2] * subpos[1]) / newcc
                            newright = (old[2] * old[1] + subpos[2] * subpos[0]) / newcc
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][1] = [newleft, newright, newcc]
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][2].append(i)
                            FusionPos[gene2 + '\t' + gene1][1][smallindex][3] += 1
            gene = gene2 + '\t' + gene1
        else:
            FusionMatrix[gene1 + '\t' + gene2] = rawscore
            FusionPos[gene1 + '\t' + gene2] = [[chromo1, chromo2], []]
            for k in range(len(subgrp)):
                FusionPos[gene1 + '\t' + gene2][1].append([subgrp[k], AvePos[k], [i], 1])
            gene = gene1 + '\t' + gene2
    infile.close()
    exprfile.close()
for gene in FusionMatrix:  
    # 这里是为FusionPos[gene][1]的LIST添加上第五个元素，表明最终决定的聚类的断点两侧坐标
    # 逻辑是必须要有一个断点在集群内占40%以上的支持，否则标记FusionPos[gene][1][i].append([-1, -1])
    temp = FusionPos[gene][1]
    for i in range(len(temp)):# 对每个聚类
        keylist = sorted(temp[i][0], key=temp[i][0].__getitem__, reverse=True)
        possum = 0
        countlist = []
        for key in keylist:
            possum += temp[i][0][key]
            countlist.append(int(temp[i][0][key])) # 每个聚类的每个位置有多少read支持都记录了
        if possum > 0 and max(countlist) / possum >= 0.4:
            for key in temp[i][0]:
                if key == keylist[0]:
                    pos0, pos1 = ParsePosition(keylist[0])
                    if abs(pos0 - temp[i][1][0]) + abs(pos1 - temp[i][1][1]) < abs(pos1 - temp[i][1][0]) + abs(
                            pos0 - temp[i][1][1]):
                    # 自动纠正 pos0 和 pos1 的左右顺序
                        FusionPos[gene][1][i].append([pos0, pos1])
                    else:
                        FusionPos[gene][1][i].append([pos1, pos0])
                    break
        else:
            FusionPos[gene][1][i].append([-1, -1])
templines = []
for gene in FusionMatrix:
    for i in range(len(FusionPos[gene][1])):
        if FusionPos[gene][1][i][4][0] > -0.5:
            genepart = gene.split('\t')
            templine = genepart[0] + '\t' + genepart[1] + '\t' + str(
                FusionMatrix[gene][i]) + '\t' + FusionPos[gene][0][0] + '\t' + FusionPos[gene][0][1] + '\t' \
                       + str(FusionPos[gene][1][i][4][0]) + '\t' + str(FusionPos[gene][1][i][4][1]) + '\t' + str(i) \
                       + '\t' + str(FusionPos[gene][1][i][3])
            templine += '\n'
            templines.append(templine)
linedic = {}
for line in templines:
    info = line.split('\t')
    score = float(info[2])
    linedic[line] = score
for key in sorted(linedic, key=linedic.__getitem__, reverse=True):
    print(key, end='')
