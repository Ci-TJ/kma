# intron ops: tools for intron operations
# Copyright (C) 2015 Harold Pimentel
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import os
import sys

import pysam

from gtf_parser import Transcript

def print_trans_names(trans):
    """
    trans:
    ENSG00000223972.5
        ENST00000456328.2
        ENST00000450305.2
    """
    if len(trans) == 0:
        return

    print trans[0].gene_id
    for t in trans:
        print '\t', t.transcript_id

def reduce_to_gene(trans_list, process_by_gene = None):
    """ 
    Performs 'process_by_gene' on a list of transcripts which all have the
    same gene names. returns the gene -> transcript mapping and returns the gene
    -> result mapping in a tuple. *IMPORTANT* Assumes that transcripts are
    sorted by position
    """
    ##“把 transcript 列表按 gene_id 分组，并对每个基因调用 process_by_gene()。返回 gene→transcript 和 gene→结果 的两个字典。注意：transcript 列表必须事先按位置排序。”

    cur_gene = None

    ret_dict = {} ## 存储 gene → process_by_gene(transcripts) 的结果
    gene_to_trans = {} ## 存储 gene → transcript 列表 的映射

    # makes no assumption about overlaps
    for trans in trans_list:
        cur_gene = trans.gene_id
        gene_list = gene_to_trans.get(cur_gene, []) ## 获取该基因已有的 transcript 列表，如果没有则返回空列表
        # XXX: need to decide what best design is here
        # if len(gene_list) > 0:
            # if gene_list[0].refname != trans.refname:
            #     print >> sys.stderr, 'ERROR: Gene showed up several times, discarding mismatch'
            #     continue
        gene_list.append(trans)
        gene_to_trans[cur_gene] = gene_list

    if process_by_gene is not None:
        for gene in gene_to_trans:
            ret_dict[gene] = process_by_gene( gene_to_trans[gene] )
        return (gene_to_trans, ret_dict)

    return gene_to_trans

class Intron:
    def __init__(self, ref, start, stop, extend_start = 0, extend_stop = 0):
        self.refname = ref
        self.coords = (start, stop)
        self.extension = (extend_start, extend_stop)

    def __getitem__(self, i):
        if i == 0:
            return self.coords[i] - self.extension[i]
        elif i == 1:
            return self.coords[i] + self.extension[i]

        # should raise an out of bounds exception
        return self.coords[i]

    @property
    def start(self):
        return self.coords[0]

    @property
    def stop(self):
        return self.coords[1]

    def __eq__(self, other):
        return self.refname == other.refname and self.coords == other.coords

    def __repr__(self):
        return '{0}:{1}-{2}'.format(self.refname, self[0],
                                    self[1])

    def to_string_noext(self):
        return '{0}:{1}-{2}'.format(self.refname, self.coords[0],
                                    self.coords[1])

    @staticmethod
    def from_string(intron_string, extend_start = 0, extend_stop = 0):
        ref, rest = intron_string.split(":")
        start, stop = rest.split("-")

        return Intron(ref, int(start) + extend_start, int(stop) - extend_stop,
                      extend_start, extend_stop)

def get_introns(trans, extend = 0):
    """ Given a transcript, return a set of introns with (start, stop)
    locations. If extend > 0, then extends the intron on the left and right
    side.  """
    ## 给定一个 transcript，返回它的 intron 列表（每个 intron 是一个区间） # extend > 0 时，会在 intron 左右两侧各延伸 extend 个碱基；文章中默认延伸25nt！！！
    if extend < 0:
        raise Exception("Non-sensical value for extend (must be >= 0)")

    introns = []
    for i in xrange(len(trans.exons) - 1):
        intron = Intron(trans.refname, trans.exons[i][1],
                        trans.exons[i + 1][0], extend, extend)
        introns.append( intron )

    return introns

def intron_all_trans(trans_list):
    """ Given a list of transcripts, return a sorted list of introns that every
    transcript shares. """
    ## 给定一个 transcript 列表，返回所有 transcript 共同拥有的 intron（交集） # 也就是说：找出所有转录本共享的 intron 区间
    if len(trans_list) == 0:
        return []
    all_introns = map(get_introns, trans_list)
    print all_introns
    # TODO: fixme ## 下面这段代码其实是错误的（KMA 作者自己也标注了 fixme） # 因为 introns.coords 并不是标准属性，但我们保持原样注释
    all_introns = [introns.coords for introns in all_introns]
    all_introns = [set(introns) for introns in all_introns]
    all_introns = list(reduce(set.intersection, all_introns))
    all_introns.sort()

    return all_introns

def intron_intersection(i1, i2):
    """ Given two introns, finds their maximal overlap """
    left = max(i1[0], i2[0])
    right = min(i1[1], i2[1])
    if right <= left:
        return None
    return (left, right)

# @profile
def transcript_union(trans_list):
    """ Given a list of transcripts, return a transcript that is the 'union' of
    them. That is, take the union overlap region of every exon."""
    # 给定多个 transcript，返回一个新的 transcript， # 其 exon 是所有 transcript 的 exon 的“并集”（union）。 # 注意：这里的 union 是“区间合并”，不是数学意义上的集合并集。

    ## 取出所有 transcript 的所有 exon，形成一个平铺列表
    all_exons = [exon for trans in trans_list for exon in trans.exons]
    all_exons = sorted(list(set(all_exons)))

    exon_union = []
    candidate = all_exons[0]
    ## 从第二个 exon 开始遍历
    for it in xrange(1, len(all_exons)):
        cur_exon = all_exons[it]
        if candidate[1] < cur_exon[1] and cur_exon[0] < candidate[1]:
            candidate = (candidate[0], cur_exon[1])
        elif candidate[1] <= cur_exon[0]:
            if candidate[1] == cur_exon[0]:
                candidate = (candidate[0], cur_exon[1])
            else:
                exon_union.append( candidate )
                candidate = cur_exon
    exon_union.append(candidate)

    ## 构造一个新的 Transcript 对象，作为 union transcript
    t0 = trans_list[0]
    strand = '-' if t0.is_reverse else '+'
    t = Transcript(t0.gene_id,
                   t0.refname,
                   strand,
                   t0.frame,
                   t0.gene_id_attributes,
                   t0.gene_id,
                   t0.score,
                   t0.source)
    t.exons = exon_union

    return t

def intron_trans_compat(intron_list, trans_list):
    """ Under the assumption that the intron_list is derived from the
    union gene from the trans_list, returns a dictionary where intron =>
    list of transcript ids that it is compatible with """
    
    """
    intron_list
    [
      Intron(refname="chr1", start=12227, end=12612),
      Intron(refname="chr1", start=12721, end=13220),
      ...
    ]

    trans_list:
    [
      Transcript(
          transcript_id="ENST00000456328.2",
          refname="chr1",
          front_coordinate=11868,
          end_coordinate=14409,
          exons=[(11868,12227),(12612,12721),(13220,14409)]
      ),
      Transcript(
          transcript_id="ENST00000450305.2",
          refname="chr1",
          front_coordinate=12010,
          end_coordinate=13670,
          exons=[(12010,12057),(12179,12227),(12612,12721)]
      ),
      ...
    ]
    """
    
    intron_to_trans = {} ## 判断 transcript 是否覆盖 intron # 条件 1：trans 起点 < intron 起点 # 条件 2：trans 终点 > intron 终点 # 条件 3：同一条染色体
    for intron in intron_list:
        key = str(intron)
        for trans in trans_list:
            if trans.front_coordinate < intron[0] and \
                    intron[1] < trans.end_coordinate and \
                    intron.refname == trans.refname:
                matches = intron_to_trans.get(key, [])
                matches.append(trans.transcript_id)
                intron_to_trans[key] = matches

    return intron_to_trans

def get_overlapping_transcripts(transcripts):
    """ Iterate through a list of transcripts. Every yield returns a list of
    transcripts that overlap in genomic coordinates. """

    transcripts = sorted(transcripts, key = lambda x: (x.refname,
                                                       x.front_coordinate,
                                                       x.end_coordinate))
    cur_end = None
    cur_ref = None
    trans_list = []
    for t in transcripts:
        if cur_end is not None and t.front_coordinate < cur_end and \
                t.refname == cur_ref:
            cur_end = max(cur_end, t.end_coordinate)
        else:
            if cur_end is not None:
                yield trans_list
            cur_end = t.end_coordinate
            cur_ref = t.refname
            trans_list[:] = []
        trans_list.append(t)

    yield trans_list

def discard_overlapping_introns(transcripts, extend = 0):
    """ Given a dictionary which maps from gene_to_union (reduce_to_gene(trans,
    transcript_union)), and gene_to_intron, remove introns that overlap with
    other genes. Side effects: gene_to_introns has an updated set of introns """
    ## 这个函数的目的： # 对每个“重叠的 transcript block”，计算它们的 union transcript， # 再生成 introns，然后删除那些跨越不同基因之间重叠区域的 introns。 
    # # 换句话说： # 如果两个基因在基因组上重叠，那么它们的 intron 可能跨基因边界， # 这些 intron 是不可信的，需要丢弃。 # # 返回：gene_id → 过滤后的 intron 列表
    gene_to_introns = {}

    for overlap_trans in get_overlapping_transcripts( transcripts ):
        # get the gene union then the introns from this union
        # ------------------------------------------------------------ 
        # 第一步：对这个簇按 gene 分组，并计算每个 gene 的 union transcript 
        # ------------------------------------------------------------ 
        # reduce_to_gene(overlap_trans, transcript_union) # 返回： 
        # g2t: gene → transcript 列表 # g2u: gene → union transcript（由 transcript_union() 生成）
        g2t, g2u = reduce_to_gene( overlap_trans, transcript_union )
        ## 第二步：对每个 gene 的 union transcript 生成 intron 列表
        g2i = { gene: get_introns(g_union, extend) for (gene, g_union) in \
               g2u.iteritems() }

        ## 第三步：把所有 gene 的 union transcript 按坐标排序 # 用于后面检测基因之间的重叠区间
        t_unions = sorted(g2u.values(), key = lambda x: (x.front_coordinate,
                              x.end_coordinate))
        intersection = []
        ## 第四步：找出所有 union transcript 之间的重叠区间
        for j in xrange(len(t_unions)):
            for k in xrange(j + 1, len(t_unions)):
                if t_unions[k].front_coordinate < t_unions[j].end_coordinate:
                    l_isect = t_unions[k].front_coordinate
                    r_isect = min(t_unions[j].end_coordinate,
                                  t_unions[k].end_coordinate)
                    intersection.append( (l_isect, r_isect) )
                else:
                    break

        # given a set of candidate intersections, again, find the union of
        # intersection
        # print 'the intersection', intersection
        intersection = unionize_regions(intersection)
        # print 'the intersection', intersection

        # we have a list of intersection regions, now need to find introns
        # that cross intersection
        ## 第五步：把所有重叠区间合并（unionize） # 例如 [(100,150),(120,180)] → [(100,180)]
        for cur_gene in t_unions:
            valid_introns = []
            for intron in g2i[cur_gene.gene_id]:
                valid_intron = True
                for invalid_region in intersection:
                    # print 'int: ', intron, ' ir ', invalid_region
                    if (intron[1] <= invalid_region[1] and \
                            invalid_region[0] < intron[1]) or \
                            (invalid_region[1] <= intron[1] and \
                             intron[0] < invalid_region[1]):
                        valid_intron = False
                        break
                if valid_intron:
                    valid_introns.append( intron )

            # print valid_introns
            introns = gene_to_introns.get(cur_gene.gene_id, [])
            introns.extend( valid_introns )
            gene_to_introns[cur_gene.gene_id] = introns

    return gene_to_introns

def unionize_regions(regions):
    """ Given a set of (sorted) regions, take the union of the overlaps """

    cand = None
    unionized = []
    for reg in regions:
        if cand is None:
            cand = reg
        elif reg[0] <= cand[1]:
            cand = (cand[0], max(reg[1], cand[1]))
        else:
            unionized.append(cand)
            cand = reg
    if cand is not None:
        unionized.append(cand)

    return unionized

def intron_all_junction_left(trans_list):
    """ Given a list of transcripts, return a sorted list of intronic regions
    that every transcript shares. """
    ## 输入：多个 transcript # 输出：所有 transcript 共同拥有的 intron 左侧交集区域 # # 注意：这是“左侧 junction 兼容 intron”，不是完整 intron 交集 # 逻辑非常 ad-hoc，作者自己也承认不准确。

    if len(trans_list) == 0:
        return []

    trans_list = sorted(trans_list, key = lambda x: x.front_coordinate)


    all_introns = [get_introns(trans) for trans in trans_list]

    if len(all_introns) == 1:
        return all_introns

    # TODO: get the union of all introns and see if it still works -- should be ### 作者承认：应该先做 intron union 再求交集，但他没做！！！
    # more accurate

    candidates_l = all_introns[0]
    mark_for_removal = []

    # iterate each transcript
    for i in xrange(1, len(all_introns)):
        # XXX: might need to move this into the next loop
        c_start = 0
        for intron in all_introns[i]:
            for j in xrange(c_start, len(candidates_l)):
                cand = candidates_l[j]
                intersection = intron_intersection(intron, cand)
                ## 情况 1：两个 intron 左端点相同
                if intron[0] == cand[0]:
                    candidates_l[j] = intersection
                    c_start += 1
                ## 情况 2：有重叠，但左端点不同
                elif intersection is not None:
                    ## 如果当前 transcript 的起点在 cand 内部 # 缩短 cand 的右端，使其左侧对齐
                    if trans_list[i].front_coordinate > cand[0]:
                        candidates_l[j] = (cand[0],
                                           trans_list[i].front_coordinate)
                    else:
                        mark_for_removal.append(j)
                ## 情况 3：intron 覆盖 cand 的左端点（左侧 intron 更大）# 这种情况不应该从左侧包含 → 删除 cand
                elif intron[0] < cand[0] and intron[1] > cand[0]:
                    # there is an overlap from an intron further on the left
                    # side. shouldn't include it from the left side but from the
                    # right
                    mark_for_removal.append(j)
                ## 情况 4：cand 完全在 intron 右侧 → 后面不会再有重叠
                elif cand[0] > intron[1]:
                    break

        candidates_l = [ c for idx, c in enumerate(candidates_l)
                        if idx not in mark_for_removal]
        mark_for_removal[:] = []

    # cleanup the introns that overlap an exonic region
    for trans in trans_list:
        for j, cand in enumerate(candidates_l):
            if trans.compatible(cand[0]) is not None or \
                trans.compatible(cand[1] - 1) is not None:
                mark_for_removal.append(j)

    candidates_l = [ c for idx, c in enumerate(candidates_l)
                    if idx not in mark_for_removal]

    return candidates_l

# mapping: dictionary where keys are gene names and values is a list of
# transcripts intersection: dictionary where keys are gene names and values are
# intronic regions that are common amonst all isoforms
# def intron_retained_transcripts(mapping, intersection):


class IntronCoverage:
    def __init__(self, ref = None, coords = None, coverage = 0,
            support = (0, 0)):
        self.ref = ref
        self.coords = coords
        self.coverage = coverage
        self.support = support

def junction_support( ps_handle, ref, intron, read_len ):
    # TODO: test me
    ## 功能：计算 intron 左右两侧 splice junction 的支持 read 数 # # ps_handle: pysam.AlignmentFile 对象（BAM 文件句柄） 
    # ref: 染色体名称 # intron: (start, end) 坐标 # read_len: read 长度（用于定义窗口大小）

    left_start = (intron[0] - 1) - read_len + 1
    left_end = (intron[0] - 1) + read_len - 2

    try:
        ## 遍历 BAM 中落在左侧窗口的 reads # read.overlap(a,b) 返回 read 在区间 [a,b] 的重叠长度 # read.rlen 是 read 的总长度 
        # # 条件：read.rlen == read.overlap(...) # 意味着 read 完整覆盖窗口（即 read 完全跨越 splice junction）
        left_count = sum(read.rlen == read.overlap(left_start, left_end)
                for read in ps_handle.fetch(ref, left_start, left_end))
    except ValueError:
        ## 如果窗口越界（负坐标），直接返回 0
        return (0, 0)

    right_start = (intron[1] - 1) - read_len + 2
    right_end = (intron[1] - 1) + read_len - 1

    right_count = sum(read.rlen == read.overlap(right_start, right_end)
            for read in ps_handle.fetch(ref, right_start, right_end))

    return (left_count, right_count)



def compute_coverage( ps_handle, ref, intron, read_len ):
    # TODO: test me
    ## 功能：计算 intron 区域的覆盖度（coverage） 
    # # 覆盖度定义为： # 覆盖 intron 扩展区间的“完整覆盖 read 数” / 区间长度 
    # # ps_handle: pysam.AlignmentFile（BAM 文件句柄） # ref: 染色体名称 # intron: (start, end) 坐标 # read_len: read 长度，用于定义扩展窗口
    left_start = intron[0] - read_len + 1
    right_end = (intron[1] - 1) + read_len - 1

    try:
        ## 遍历 BAM 中落在 [left_start, right_end] 的 reads # # 条件 1：read.rlen == read.overlap(left_start, right_end) # → read 必须完整覆盖整个窗口 
        # # 条件 2：(read.aend - read.pos) == read.rlen # → read 没有 soft-clip / indel，是真正的全长匹配 # # 满足两个条件的 read 才计数
        count = sum(read.rlen == read.overlap(left_start, right_end) and
                    (read.aend - read.pos) == read.rlen for read in
                    ps_handle.fetch(ref, left_start, right_end))
        cov = float(count) / (right_end - left_start + 1)
    except ValueError:
        return 0.0

    return cov

# bam_fname: a BAM file name
#
def bam_to_measurable(bam_fname, gene_to_trans, gene_to_introns):
    # TODO: test me
    ## 功能： # 给定 BAM 文件、gene→transcript 映射、gene→introns 映射， 
    # 计算每个 intron 的 coverage 和 junction support， 
    # 并返回： # 1. 每个基因覆盖度最高的 intron # 2. 每个基因所有 intron 的测量结果（IntronCoverage 对象）

    bam_handle = pysam.Samfile(bam_fname, 'rb')
    tmp = bam_handle.next()
    read_len = tmp.rlen

    gene_to_max_introns = {}
    gene_to_measurable_introns = {}

    for gene_name, all_introns in gene_to_introns.iteritems():
        if len(all_introns) == 0:
            continue
        max_intron = None
        max_cov = 0.0

        measurable_introns = []

        ref = gene_to_trans[gene_name][0].refname
        for intron in all_introns:
            cur_cov = None
            junc_supp = junction_support(bam_handle, ref, intron, read_len) ## 计算左右剪接点支持度（left_count, right_count）
            cov = compute_coverage(bam_handle, ref, intron, read_len)
            cur_cov = IntronCoverage(ref, intron, cov, junc_supp)

            ## 更新最大覆盖度 intron
            if cov > max_cov:
                max_intron = cur_cov
                max_cov = cov
            measurable_introns.append( cur_cov )

            # For now compute for every intron
            # if sum(junc_supp) > 0:
            #     cov = compute_coverage(ref, intron)
            #     cur_cov = IntronCoverage(ref, intron, cov, True)
            # else:
            #     cur_cov = IntronCoverage(ref, intron)

        ## 保存该基因的最大 intron 和所有 intron
        gene_to_max_introns[ gene_name ] = max_intron
        gene_to_measurable_introns[ gene_name ] = measurable_introns

    bam_handle.close()

    return (gene_to_max_introns, gene_to_measurable_introns)

def print_measurable(gene2trans, gene2intersect, gene2max, gene2measurable, handle):
    """
    把每个基因的所有 intron 的测量结果（coverage + junction support）打印成一张表。
    gene    reference   start   end   coverage   support_start   support_end   transcripts
    ENSG000001   chr1   1000   1200   0.03   5   3   ENST00001,ENST00002
    """
    print >> handle, 'gene\treference\tstart\tend\tcoverage\tsupport_start\tsupport_end\ttranscripts'
    for gene in gene2max:
        for measurable in gene2measurable[gene]:
            cur_line = []
            cur_line.append(gene)
            cur_line.append( measurable.ref )
            cur_line.append( str(measurable.coords[0]) )
            cur_line.append( str(measurable.coords[1]) )
            cur_line.append( str(measurable.coverage) )
            cur_line.append( str(measurable.support[0]) )
            cur_line.append( str(measurable.support[1]) )
            cur_line.append( ','.join(trans.transcript_id
                for trans in gene2trans[gene]) )
            print >> handle, '\t'.join( cur_line )
