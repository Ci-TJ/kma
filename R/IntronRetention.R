# Keep Me Around: Intron Retention Detection
# Copyright (C) 2015  Harold Pimentel <haroldpimentel@gmail.com>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# TODO: update docs below...

#' Compute intron retention
#'
#' Create an IntronRetention object.
#'
#' @param targetExpression a matrix of nothing but expression values and one
#' column containing 'target_id's
#' @param intronToUnion a table with columns (all character vectors):
#' \describe{
#'  \item{intron}{the identifier of the intron}
#'  \item{target_id}{the target_id of transcripts compatible with the intron}
#'  \item{gene}{the gene name that this intron belongs to}
#'  \item{intron_extension}{the actual coordinates of the intron quantified
#'  (including the region that overlaps the intronic region)}
#' }
#' @param groups a vector with the grouping
#' @param psi if TRUE, compute the psi value. otherwise, compute a rate
#' @return an IntronRetention object
#' @export
newIntronRetention <- function(targetExpression,
    intronToUnion,
    groups,
    unique_counts = NULL,
    psi = TRUE)
{
    # targetExpression：每一行是一个 target_id（转录本或 intron_extension），列是样本表达量 # intronToUnion：每一行是 intron 与兼容的 target_id 的对应关系 
    # groups：每个样本对应的 condition # unique_counts：可选，来自 eXpress 的 unique fragment 数 
    # psi：是否把 intron_extension 当作一个“伪转录本”加入表达矩阵，用于计算 PSI
    # TODO: verify all 'introns' are in targetExpression and all target_ids in
    # targetExpression
    labs <- setdiff(colnames(targetExpression), 'target_id') ## labs = 所有样本名（target_id 除外）

    if (length(groups) != length(labs)) {
        stop("length(groups) must be the same as the number of experiments included (and also in the same order)")
    }

    # targetExpression <- data.table(targetExpression)
    # intronToUnion <- data.table(intronToUnion)
    targetExpression <- targetExpression %>%
        arrange(target_id)

    intronToUnion <- intronToUnion %>%
        arrange(target_id)

    ## 如果 psi=TRUE，则把 intron_extension 当作一个新的 target_id 加入 # 这是 KMA 的关键：把 intron 本身当作一个“转录本”来量化
    if (psi) {
        repIntrons <- intronToUnion %>%
            select(intron, gene, intron_extension) %>%
            distinct() %>%
            mutate(target_id = intron_extension) ## 每个 intron 只保留一条 mutate(target_id = intron_extension) # intron_extension 作为 target_id
        intronToUnion <- data.table(dplyr::bind_rows(intronToUnion, repIntrons))
    }
    
    unique_counts_tbl <- NULL

    ## 如果提供 unique_counts，则需要把它 melt 成 long format，并映射到 intron
    if (!is.null(unique_counts)) {
        # TODO: verify column names are exactly the same as in targ_expression
        cat("'melting' unique counts\n")
        ## intron_extension → target_id，用于和 unique_counts 对齐
        intron_targ_tbl <- intronToUnion %>%
            select(intron, intron_extension) %>%
            distinct() %>%
            rename(target_id = intron_extension)
        ## unique_counts 原本是宽表（每列一个样本），这里 melt 成 long
        unique_counts_tbl <- melt(unique_counts, id.vars = "target_id",
            variable.name = "sample",
            value.name = "unique_counts")

        # return(list(unique_counts = unique_counts_tbl, intron_targ = intron_targ_tbl))
        ## 把 unique_counts 映射到 intron
        unique_counts_tbl <- data.table(unique_counts_tbl) %>%
            inner_join(data.table(intron_targ_tbl), by = c("target_id")) %>%
            select(-c(target_id)) # %>%
            # mutate(sample = as.character(sample))
    }


    cat("computing denominator\n")
    intronToUnion <- data.table(intronToUnion)
    targetExpression <- data.table(targetExpression)
    ## denomExp：每个 intron 的 denominator（所有兼容转录本表达量之和）
    denomExp <- left_join(intronToUnion, targetExpression, by = "target_id") %>%
        group_by(intron) %>%
        select(-(target_id)) %>%
        summarise_each(funs(sum), -matches("gene"),
            -matches("intron_extension")) %>%
        arrange(intron) %>%
        left_join(
            select(intronToUnion, intron, intron_extension) %>%
                distinct(),
            by = c("intron")) ## target_id 不再需要# gene 和 intron_extension 不参与求和# 把 intron_extension 加回来

    ## tmp_targExpression：把 target_id 改名为 intron_extension，用于 numerator 计算
    tmp_targExpression <- targetExpression %>%
        rename("intron_extension" = target_id) %>%
        data.table()

    denomExp <- data.table(denomExp)
    cat("computing numerator\n")
    ## numerator：每个 intron_extension 对应的表达量（即 intron 本身的表达）
    numExp <- select(denomExp, intron, intron_extension) %>%
        # inner_join(targetExpression, by = c("intron_extension" = "target_id")) %>%
        inner_join(tmp_targExpression, by = c("intron_extension")) %>%
        arrange(intron_extension)
    rm(tmp_targExpression)

    ## denomExp 和 numExp 转为矩阵形式，行名为 intron
    denomExp <- as.data.frame(denomExp) %>%
        arrange(intron)
    rownames(denomExp) <- denomExp$intron
    denomExp <- select(denomExp, -c(intron, intron_extension))

    numExp <- as.data.frame(numExp) %>%
        arrange(intron)
    rownames(numExp) <- numExp$intron
    numExp <- select(numExp, -c(intron, intron_extension))

    cat("computing retention\n") 
    retentionExp <- numExp / denomExp ## retention = numerator / denominator

    ## targetExpression 也去掉 target_id 列，行名设为 target_id
    rownames(targetExpression) <- targetExpression$target_id
    targetExpression$target_id <- NULL

    # TODO: include intron_extension in flat
    cat("'melting' expression\n")
    ## flat：长格式表，每行是 intron × sample
    flat <- melt_retention(retentionExp, numExp, denomExp, groups)

    flat <- data.table(flat)

    ## 如果有 unique_counts，则加入 flat
    if (!is.null(unique_counts)) {
        cat("joining unique_counts and retention data\n")
        flat <- flat %>%
            inner_join(unique_counts_tbl, by = c("intron", "sample"))
    }

    cat("sorting and grouping by (intron, condition)\n")
    flat <- flat %>%
        arrange(intron, condition) %>%
        group_by(intron, condition)

    ## intron_to_ext：每个 intron 的 intron_extension 及其长度
    intron_to_ext <- intronToUnion %>%
        select(intron, intron_extension) %>%
        distinct() %>%
        mutate(extension_len = intron_length(intron_extension))
    ## intron_to_gene：每个 intron 的 gene
    intron_to_gene <- intronToUnion %>%
        select(intron, gene) %>%
        distinct() %>%
        data.table()
    ## 合并 gene 信息
    intron_to_ext <- intron_to_ext %>%
        left_join(intron_to_gene, by = c("intron"))

    # TODO: add a list "filters" which keeps track of all the filters and their
    # calls
    ## 最终把所有对象转为 data.frame 并打包成 IntronRetention 对象
    retentionExp <- as.data.frame(retentionExp, stringsAsFactors = FALSE)
    numExp <- as.data.frame(numExp, stringsAsFactors = FALSE)
    denomExp <- as.data.frame(denomExp, stringsAsFactors = FALSE)
    targetExpression <- as.data.frame(targetExpression, stringsAsFactors = FALSE)
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)
    intron_to_ext <- as.data.frame(intron_to_ext, stringsAsFactors = FALSE)
    intronToUnion <- as.data.frame(intronToUnion, stringsAsFactors = FALSE)

    structure(list(retention = retentionExp,         # intron × sample 的 IR/PSI
            numerator = numExp,                      # numerator
            denominator = denomExp,                  # denominator
            labels = labs,                           # 样本名
            groups = groups,                         # condition
            features = targetExpression,             # 原始表达矩阵
            flat = flat,                             # 长格式表
            unique_counts = unique_counts,           # unique fragment 数
            intron_to_extension = intron_to_ext,     # intron → intron_extension
            intron_to_union = intronToUnion          # intron → target_id 映射
            ),
        class = "IntronRetention")
}

#' @export
melt_retention <- function(ret, num, denom, groupings)
{
    samp_to_condition <- data.frame(sample = colnames(ret),
        condition = groupings, stringsAsFactors = FALSE)

    ret <- ret %>% mutate(intron = rownames(ret))

    ret <- reshape2::melt(ret, id.vars = "intron",
        variable.name = "sample",
        value.name = "retention") %>%
        mutate(sample = as.character(sample))

    denom <- denom %>% mutate(intron = rownames(denom))
    denom <- reshape2::melt(denom, id.vars = "intron",
        variable.name = "sample",
        value.name = "denominator") %>%
        mutate(sample = as.character(sample))

    num <- num %>% mutate(intron = rownames(num))
    num <- melt(num, id.vars = "intron",
        variable.name = "sample",
        value.name = "numerator") %>%
        mutate(sample = as.character(sample))

    num <- data.table(num)
    denom <- data.table(denom)
    ret <- data.table(ret)

    m_res <- inner_join(num, denom, by = c("intron", "sample"))
    m_res <- m_res %>%
        inner_join(ret, by = c("intron", "sample"))

    samp_to_condition <- data.table(samp_to_condition)
    m_res <- data.table(m_res)
    left_join(m_res, samp_to_condition, by = "sample")
}


#' @export
retentionTestSingleCond <- function(retentionMat, level = 0.0, offset = 0.00)
{
    stopifnot(ncol(retentionMat) > 1)
    m <- apply(retentionMat, 1, mean)
    v <- apply(retentionMat, 1, var)
    v <- v + offset
    testStat <- (m - level) / sqrt(v)
    list(avg = m, variance = v, testStat = testStat)
}

#' Get intron lengths from an identifier
#'
#' Transform an intron identifier into lengths.
#'
#' @param intron_names the names of the introns in format 'chrom:start-stop'
#' @return an integer vector of the lengths
#' @export
intron_length <- function(intron_names)
{
    unlist(lapply(strsplit(intron_names, ":"), function(x)
        {
            coords <- as.integer(strsplit(x[2], '-')[[1]])
            coords[2] - coords[1]
        }))
}

#' Generate the null distribution
#'
#'
#' @param flat_grouped
#' @param n_samp number of samples
#' @param test_stat test statistic to use (function)
#' @return a numeric with means
#' @export
intron_null_dist <- function(flat_grouped, n_samp = 10000, test_stat = mean)
{
    all_dat <- dcast(flat_grouped, intron ~ sample, value.var = "retention")
    all_dat <- select(all_dat, -c(intron))
    all_dat <- all_dat[complete.cases(all_dat),]

    samps <- as.matrix(as.data.frame(lapply(all_dat, sample,
                n_samp,
                replace = T)))

    data <- apply(samps, 1, test_stat)

    list(data = data, ecdf = ecdf(data))
}

#' @export
intron_pval <- function(mean_val, null_ecdf)
{
    1 - null_ecdf(mean_val)
}

#' Print an IntronRetention object
#'
#' Print an IntronRetention object
#' @param ir an IntronRetention object
#' @return the unchanged object \code{ir} after printing a summary to the
#' terminal
#' @export
print.IntronRetention <- function(ir)
{
    cat(sprintf("IntronRetention object (%d introns)\n",
            nrow(ir$retention)))
    cat("----------------------------------------------\n")
    cat("Samples:\t", paste(ir$labels, collapse = " "), "\n", sep = "")
    cat("Conditions:\t", paste(ir$groups, collapse = " "), "\n", sep = "")
    invisible(ir)
}

#' @export
add_gene_names <- function(dat, ir)
{
    # TODO: ensure "intron" exists in dat

    nrow_before <- nrow(dat)
    dat <- data.table(dat) %>%
        left_join(data.table(ir$intron_to_ext), by = "intron")

    nrow_after <- nrow(dat)
    stopifnot(nrow_before == nrow_after)

    dat
}
