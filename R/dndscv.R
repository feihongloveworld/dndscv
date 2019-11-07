#' dNdScv
#'
#' Analyses of selection using the dNdScv and dNdSloc models. Default parameters typically increase the performance of the method on cancer genomic studies. Reference files are currently only available for the GRCh37/hg19 version of the human genome.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param mutations Table of mutations (5 columns: sampleID, chr, pos, ref, alt). Only list independent events as mutations.
#' @param gene_list List of genes to restrict the analysis (use for targeted sequencing studies)
#' @param refdb Reference database (path to .rda file)
#' @param sm Substitution model (precomputed models are available in the data directory)
#' @param kc List of a-priori known cancer genes (to be excluded from the indel background model)
#' @param cv Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates]
#' @param max_muts_per_gene_per_sample If n<Inf, arbitrarily the first n mutations by chr position will be kept
#' @param max_coding_muts_per_sample Hypermutator samples often reduce power to detect selection
#' @param use_indel_sites Use unique indel sites instead of the total number of indels (it tends to be more robust)
#' @param min_indels Minimum number of indels required to run the indel recurrence module
#' @param maxcovs Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)
#' @param constrain_wnon_wspl This constrains wnon==wspl (this typically leads to higher power to detect selection)
#' @param outp Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate. Note that the same genetic code must be used in the dndscv and buildref functions.
#' @param outmats Output the internal N and L matrices (default = F)
#'
#' @return 'dndscv' returns a list of objects:
#' @return - globaldnds: Global dN/dS estimates across all genes.
#' @return - sel_cv: Gene-wise selection results using dNdScv.
#' @return - sel_loc: Gene-wise selection results using dNdSloc.
#' @return - annotmuts: Annotated coding mutations.
#' @return - genemuts: Observed and expected numbers of mutations per gene.
#' @return - mle_submodel: MLEs of the substitution model.
#' @return - exclsamples: Samples excluded from the analysis.
#' @return - exclmuts: Coding mutations excluded from the analysis.
#' @return - nbreg: Negative binomial regression model for substitutions.
#' @return - nbregind: Negative binomial regression model for indels.
#' @return - poissmodel: Poisson regression model used to fit the substitution model and the global dNdS values.
#' @return - wrongmuts: Table of input mutations with a wrong annotation of the reference base (if any). 
#' 
#' @export

dndscv = function(mutations, 
                  gene_list = NULL, 
                  refdb = "hg19", 
                  sm = "192r_3w", 
                  kc = "cgc81", 
                  cv = "hg19", 
                  max_muts_per_gene_per_sample = 3, 
                  max_coding_muts_per_sample = 3000, 
                  use_indel_sites = T, 
                  min_indels = 5, 
                  maxcovs = 20, 
                  constrain_wnon_wspl = T, 
                  outp = 3, 
                  numcode = 1, 
                  outmats = F) {

    ## 1. Environment
    message("[1] Loading the environment...")
    ## 突变提取前五列
    mutations = mutations[,1:5] # Restricting input matrix to first 5 columns
    ## 转化为字符串
#sampleID                        chr  pos       ref  mut
#F16061220030-CLN                17   7577058   C    A
#F16061220030-CLN                3    47155465  C    T
#ct-C1512289427-S-CLN            17   7578441   G    T
#ct-C16042816301-S-CLN           17   7578441   G    T
#ct-C16060119301-S-CLN           17   7578441   G    T
#ct-C16060119301-S-YH095-CLN-RE  17   7578441   G    T
#ct-C16092028993-S-CLN           2    47601106  TG   CA
#ct-C16092028993-S-CLN           17   7578206   T    C
#F16053119231-CLN                17   7577538   C    T

    mutations[,c(1,2,3,4,5)] = lapply(mutations[,c(1,2,3,4,5)], as.character) # Factors to character
    ## 变异位点转化为数值型 
    mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
    ## 提取发生真实碱基突变的位点
    mutations[[3]] = as.numeric(mutations[[3]]) # Removing mutations with identical reference and mutant base
    ## 添加表头
    colnames(mutations) = c("sampleID","chr","pos","ref","mut")
    ## 提取存在NA的行，返回index
    # Removing NA entries from the input mutation table
    indna = which(is.na(mutations),arr.ind=T)
    ## 整行去除含NA的突变
    if (nrow(indna)>0) {
        mutations = mutations[-unique(indna[,1]),] # Removing entries with an NA in any row
        warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(indna[,1]))))
    }
    ## 载入参考数据库
    # [Input] Reference database
    if (refdb == "hg19") {
        data("refcds_hg19", package="dndscv")
        if (any(gene_list=="CDKN2A")) { # Replace CDKN2A in the input gene list with two isoforms  为啥
            gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
        }
    } else {
        load(refdb)
    }
    ## 输入基因列表，默认全部基因
    # [Input] Gene list (The user can input a gene list as a character vector)
    ## 
    if (is.null(gene_list)) {
        gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
    } else { # Using only genes in the input gene list
        ## 首先提取全部基因
        allg = sapply(RefCDS,function(x) x$gene_name)
        ## 基因必须存在于全部基因中
        nonex = gene_list[!(gene_list %in% allg)]
        ## 打印出不在列的基因
        if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
        RefCDS = RefCDS[allg %in% gene_list] # Only input genes
        ## gr_genes
        gr_genes = [gr_genes$names %in% gene_list] # Only input genes
    }
    ## 协变量矩阵
    # [Input] Covariates (The user can input a custom set of covariates as a matrix)
    if (is.character(cv)) {
        data(list=sprintf("covariates_%s",cv), package="dndscv")
    } else {
        covs = cv
    }
    ## 已知癌症基因 默认使用CGC81数据库
    # [Input] Known cancer genes (The user can input a gene list as a character vector)
    if (kc[1] %in% c("cgc81")) {
        data(list=sprintf("cancergenes_%s",kc), package="dndscv")
    } else {
        known_cancergenes = kc
    }
    ## 突变模型  默认 192r_3w 用户可以提供一个个性化的突变模型作为举证
    # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
    if (length(sm)==1) {
        data(list=sprintf("submod_%s",sm), package="dndscv")
    } else {
        substmodel = sm
    }
    ## 扩展的参考序列 是为了更快的访问，这部分是高级数据结构
    # Expanding the reference sequences [for faster access]
    for (j in 1:length(RefCDS)) {
        RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
        RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
        RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
        if (!is.null(RefCDS[[j]]$seq_splice)) {
            RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
            RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
            RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
        }
    }
    
    ## 突变注释
    ## 2. Mutation annotation
    message("[2] Annotating the mutations...")
    ## 三碱基构建
    nt = c("A","C","G","T")
    trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
    trinucinds = setNames(1:64, trinucs)
    trinucsubs = NULL
    for (j in 1:length(trinucs)) {
        trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
    }
    trinucsubsind = setNames(1:192, trinucsubs)
    ## 
    ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))  ## 每个基因名对应再RefCDS上的位置编号记录下来
    ## 表达的基因
    gr_genes_ind = ind[gr_genes$names]   ## 查看gr基因基因名对饮的ind   
    ## 警告：可能存在多碱基突变注释失败的情况,删除complex突变
    # Warning about possible unannotated dinucleotide substitutions
    if (any(diff(mutations$pos)==1)) {
        warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
    }
    ## 警告：在不同样本中同一个突变不同实例
    # Warning about multiple instances of the same mutation in different sampleIDs
    if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
        warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
    }
    ## 突变的起始点
    # Start and end position of each mutation
    mutations$end = mutations$start = mutations$pos
    l = nchar(mutations$ref)-1 # Deletions of multiple bases
    mutations$end = mutations$end + l
    ind = substr(mutations$ref,1,1)==substr(mutations$mut,1,1) & nchar(mutations$ref)>nchar(mutations$mut) # Position correction for deletions annotated in the previous base (e.g. CA>C)
    mutations$start = mutations$start + ind
    ## 突变对应到基因上
    # Mapping mutations to genes
    gr_muts = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$start,mutations$end))
    ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))
                                            
    mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene  ## 如果对应到多个基因，删除第二个#
    mutations$geneind = gr_genes_ind[ol[,2]]  ## 添加上对应基因再RefDCS上的位置
    mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind] ## 添加上基因名
    mutations = unique(mutations)
    ## 排除超突变样本 max_coding_muts_per_sample
    # Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
    nsampl = sort(table(mutations$sampleID))
    exclsamples = NULL
    if (any(nsampl>max_coding_muts_per_sample)) {
        message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). %0.0f samples left after filtering.',sum(nsampl>max_coding_muts_per_sample),sum(nsampl<=max_coding_muts_per_sample)))
        exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
        mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
    }
    ## 排除超突变基因
    # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
    mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
    exclmuts = NULL
    if (any(mutrank>max_muts_per_gene_per_sample)) {
        message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample (see the max_muts_per_gene_per_sample argument in dndscv)',sum(mutrank>max_muts_per_gene_per_sample)))
        exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
        mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
    }
    ## 对突变的额外的注释
    # Additional annotation of substitutions
    ## 转录本的链信息
    mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
    ##  snv
    snv = (mutations$ref %in% nt & mutations$mut %in% nt)
    if (!any(snv)) { stop("Zero coding substitutions found in this dataset. Unable to run dndscv. Common causes for this error are inputting only indels or using chromosome names different to those in the reference database (e.g. chr1 vs 1)") }
    indels = mutations[!snv,]
    ## 取出snv
    mutations = mutations[snv,]
    mutations$ref_cod = mutations$ref
    mutations$mut_cod = mutations$mut
    compnt = setNames(rev(nt), nt)
    isminus = (mutations$strand==-1)
    mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
    mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
    ##############################################################################################
    for (j in 1:length(RefCDS)) {
        RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
        ## 再 N列表中添加192 行 + 4列的矩阵
    }
    ## 子功能：获得编码位置 在给定的外显子区间
    # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
    
    chr2cds = function(pos,cds_int,strand) {
        if (strand==1) {
            return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))  ## 返回当前位置再外显子上的坐标
        } else if (strand==-1) {
            return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
        }
    }
    ## 
    # Annotating the functional impact of each substitution and populating the N matrices
    ## 每个突变简历array
    ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = codonsub = array(NA, nrow(mutations))
    
    for (j in 1:nrow(mutations)) {
    
        geneind = mutations$geneind[j]  ## 基因编号
        pos = mutations$pos[j]          ## 对应的位置
        ## 注释潜在的剪切位点
        if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution
        
            impact[j] = "Essential_Splice"; impind = 4
            pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)  ## 当前位置所在的剪切位点的编号
            cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]    ## 剪切位点的碱基
             ## ref三碱基               ## 剪切位点的上游一个碱基            当前碱基                    下游一个碱基
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            ##   -1碱基 突变碱基  +1碱基
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            aachange[j] = ntchange[j] = codonsub[j] = "."

        } else { # Coding substitution
        
            pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
            cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]   ## cds上的碱基
           ##    cds前一个碱基  ref碱基  后一个碱基
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
           ## cds前一个碱基  突变碱基 后一个碱基
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
           ## 密码子的位置
            codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
           ## 取得ref的密码子碱基
            old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
           ## 再cds上的位置
            pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
           ## 突变后的密码子
            new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
           ## 突变前的氨基酸
            old_aa = seqinr::translate(old_codon, numcode = numcode)
           ## 突变后的氨基酸
            new_aa = seqinr::translate(new_codon, numcode = numcode)
           ## ref氨基酸  氨基酸位置  老氨基酸
            aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
           ## 碱基变化  突变碱基  突变碱基
            ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
           ## G>C
            codonsub[j] = sprintf('%s>%s',paste(old_codon,collapse=""),paste(new_codon,collapse=""))
           ## 如果新氨基酸和旧氨基酸，同义突变
            # Annotating the impact of the mutation
            if (new_aa == old_aa){ 
                impact[j] = "Synonymous"; impind = 1
            # 新氨基酸为*  沉默突变
            } else if (new_aa == "*"){
                impact[j] = "Nonsense"; impind = 3
            # 老氨基酸不为* 错译突变
            } else if (old_aa != "*"){
                impact[j] = "Missense"; impind = 2
            # 老氨基酸为* 终止子丢失
            } else if (old_aa=="*") {
                impact[j] = "Stop_loss"; impind = NA
            }
        }
        
        if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
            wrong_ref[j] = 1  ## 如果ref碱基和cds碱基不一致，位点存放到防止到wrong_ref 中
        } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
          ## 影响不为空
            ## 三碱基替换 CAA>CAT  得到三碱基编号
            trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
            ##                   三碱基编号  影响编号           +1 
            RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
        }
        ### 打印提示信息
        if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g%% ...', round(j/nrow(mutations),2)*100)) }
    }
    
    mutations$ref3_cod = ref3_cod
    mutations$mut3_cod = mut3_cod
    mutations$aachange = aachange
    mutations$ntchange = ntchange
    mutations$codonsub = codonsub  ## 碱基替换
    mutations$impact = impact
    mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]  ## 添加蛋白编号
    ## 
    错误注释位点
    如果少于10%，
    
    if (any(!is.na(wrong_ref))) {
        if (mean(!is.na(wrong_ref)) < 0.1) { # If fewer than 10% of mutations have a wrong reference base, we warn the user
            warning(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base (see the affected mutations in dndsout$wrongmuts). Please identify the causes and rerun dNdScv.', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
        } else { # If more than 10% of mutations have a wrong reference base, we stop the execution (likely wrong assembly or a serious problem with the data)
            stop(sprintf('%0.0f (%0.2g%%) mutations have a wrong reference base. Please confirm that you are not running data from a different assembly or species.', sum(!is.na(wrong_ref)), 100*mean(!is.na(wrong_ref))))
        }
    ## 筛选正确注释的位点
        wrong_refbase = mutations[!is.na(wrong_ref), 1:5]
        mutations = mutations[is.na(wrong_ref),]   
    }
    ## 看是否存在indel信息；
    if (any(nrow(indels))) { # If there are indels we concatenate the tables of subs and indels
        indels = cbind(indels, data.frame(ref_cod=".", mut_cod=".", ref3_cod=".", mut3_cod=".", aachange=".", ntchange=".", codonsub=".", impact="no-SNV", pid=sapply(RefCDS,function(x) x$protein_id)[indels$geneind]))
        
        # Annotation of indels
        ## 判断ins 还是 del
        ins = nchar(gsub("-","",indels$ref))<nchar(gsub("-","",indels$mut))
        del = nchar(gsub("-","",indels$ref))>nchar(gsub("-","",indels$mut))
        ## 双核苷酸变异
        multisub = nchar(gsub("-","",indels$ref))==nchar(gsub("-","",indels$mut)) # Including dinucleotides
        ## ref 和 mut变化的碱基说
        l = nchar(gsub("-","",indels$ref))-nchar(gsub("-","",indels$mut))
        ## 构建每个indels的字符
        indelstr = rep(NA,nrow(indels))
        for (j in 1:nrow(indels)) {
            ## indel 基因编号
            geneind = indels$geneind[j]
            ## indels 的起始位置，连续的
            pos = indels$start[j]:indels$end[j]
            ## 如果时ins，pos-1 pos
            if (ins[j]) { pos = c(pos-1,pos) } # Adding the upstream base for insertions
            ## pos的位置 cds的位置
            pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
            if (length(pos_ind)>0) {
                ## 框内
                inframe = (length(pos_ind) %% 3) == 0
                ## ins frshift infram
                ## del frshift infram
                ## 双核苷酸改变 MNV
                if (ins[j]) { # Insertion
                    indelstr[j] = sprintf("%0.0f-%0.0f-ins%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
                } else if (del[j]) { # Deletion
                    indelstr[j] = sprintf("%0.0f-%0.0f-del%s",min(pos_ind),max(pos_ind),c("frshift","inframe")[inframe+1])
                } else { # Dinucleotide and multinucleotide changes (MNVs)
                    indelstr[j] = sprintf("%0.0f-%0.0f-mnv",min(pos_ind),max(pos_ind))
                }
            }
        }
        ## 标记indel的突变类型 添加到碱基变化中
        indels$ntchange = indelstr
        ## 合并snv 和 indels
        annot = rbind(mutations, indels)
    } else {
        annot = mutations
    }
    ## 样本id  chr  pos 排序
    annot = annot[order(annot$sampleID, annot$chr, annot$pos),]
    
    ## 评估 全局率 和 选择参数
    ## 3. Estimation of the global rate and selection parameters
    message("[3] Estimating global rates...")
    ## 
    Lall = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS)))
    ## 
    Nall = array(sapply(RefCDS, function(x) x$N), dim=c(192,4,length(RefCDS)))
    ## 计算三碱基变异的不同突变类型的所有碱基的综合，这个数据的来源时啥？
    L = apply(Lall, c(1,2), sum)  ## 二维矩阵
    ## 计算当前数据集合中的对应的信息
    N = apply(Nall, c(1,2), sum)
    ## 子功能  拟合突变模型
    # Subfunction: fitting substitution model
    
    fit_substmodel = function(N, L, substmodel) {
        ## 将矩阵拉伸成向量
        l = c(L); 
        n = c(N); 
        r = c(substmodel)
        n = n[l!=0]; 
        r = r[l!=0]; 
        l = l[l!=0];
        ## 
        params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
      
        indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
        colnames(indmat) = params
        for (j in 1:length(r)) {
            indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
        }
        ## 广义线性回归模型    
        ##                            offset 偏置量  以log(l)作为偏置量
        model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
        ## exp，自然对数e为底指数函数，全称Exponential(指数曲线)。
        ## 提取系数
        mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
        ## Wald置信区间 ？？？
        ci = exp(confint.default(model)) # Wald confidence intervals
        par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
        ## 返回参数和模型
        return(list(par=par, model=model))
    }
    ## 拟合所有突变率 和 3个全局选择参数
    # Fitting all mutation rates and the 3 global selection parameters
    ## 原始的突变模型，泊松分布
    poissout = fit_substmodel(N, L, substmodel) # Original substitution model
    ## 分布参数
    par = poissout$par
    ## 分布模型
    poissmodel = poissout$model
    ## 参数的第一列为第二列的表头
    parmle =  setNames(par[,2], par[,1])
    ## 最大似然子模型
    mle_submodel = par
    rownames(mle_submodel) = NULL
    ## 拟合模型 用1 和 2 全局选择参数
    # Fitting models with 1 and 2 global selection parameters
    
    s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
    par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
    s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
    par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
    globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
    sel_loc = sel_cv = NULL
                        
                        
                        
                        
                        
                        
                        
                        
                        
    ## 可变的 dnds 模型率 ： 基因突变率 从 单个基因的同义突变中推测出来
    ## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from synonymous subs in the gene only)
    ## 只从当前基因的同义突变中推测基因的突变率
    
    genemuts = data.frame(gene_name = sapply(RefCDS, 
                                             function(x) x$gene_name), 
                                             n_syn=NA, 
                                             n_mis=NA, 
                                             n_non=NA, 
                                             n_spl=NA, 
                                             exp_syn=NA, 
                                             exp_mis=NA, 
                                             exp_non=NA, 
                                             exp_spl=NA, 
                                             stringsAsFactors=F)
    genemuts[,2:5] = t(sapply(RefCDS, function(x) colSums(x$N)))
    mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
    genemuts[,6:9] = t(sapply(RefCDS, function(x) colSums(x$L*mutrates)))
    numrates = length(mutrates)
    
    if (outp > 1) {
        message("[4] Running dNdSloc...")
        
        selfun_loc = function(j) {
            y = as.numeric(genemuts[j,-1])
            x = RefCDS[[j]]
            
            # a. Neutral model: wmis==1, wnon==1, wspl==1
            mrfold = sum(y[1:4])/sum(y[5:8]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) # loglik null model
            
            # b. Missense model: wmis==1, free wnon, free wspl
            mrfold = max(1e-10, sum(y[c(1,2)])/sum(y[c(5,6)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
            llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) # loglik free wmis
            
            # c. free wmis, wnon and wspl
            mrfold = max(1e-10, y[1]/y[5]) # Correction factor of "t"
            w = y[2:4]/y[6:8]/mrfold; w[y[2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
            llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))), log=T)) # loglik free wmis, wnon, wspl
            w[w>1e4] = 1e4
            
            p = 1-pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
            return(c(w,p))
        }
        
        sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), selfun_loc)))
        colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","pall_loc")
        sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method="BH")
        sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method="BH")
        sel_loc = cbind(genemuts[,1:5],sel_loc)
        sel_loc = sel_loc[order(sel_loc$pall_loc,sel_loc$pmis_loc,-sel_loc$wmis_loc),]
    }
    
    ## 阴性 的二项式回归 +局部的同义突变
    ## 5. dNdScv: Negative binomial regression (with or without covariates) + local synonymous mutations
    
    nbreg = nbregind = NULL
    if (outp > 2) {
        
        message("[5] Running dNdScv...")
        ## 协变量
        # Covariates
        if (is.null(cv)) {
            nbrdf = genemuts[,c("n_syn","exp_syn")]
            model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
            message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))
        } else {
            covs = covs[genemuts$gene_name,]
            if (ncol(covs) > maxcovs) {
                warning(sprintf("More than %s input covariates. Only the first %s will be considered.", maxcovs, maxcovs))
                covs = covs[,1:maxcovs]
            }
            nbrdf = cbind(genemuts[,c("n_syn","exp_syn")], covs)
            ## 阴性二项式回归
            # Negative binomial regression for substitutions
            if (nrow(genemuts)<500) { # If there are <500 genes, we run the regression without covariates
                model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
            } else {
                model = tryCatch({
                    MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + . , data = nbrdf) # We try running the model with covariates
                }, warning = function(w){
                    MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf) # If there are warnings or errors we run the model without covariates
                }, error = function(e){
                    MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf) # If there are warnings or errors we run the model without covariates
                })
            }
            message(sprintf("    Regression model for substitutions (theta = %0.3g).", model$theta))
        }
        if (all(model$y==genemuts$n_syn)) {
            genemuts$exp_syn_cv = model$fitted.values
        }
        theta = model$theta
        nbreg = model
        ## 子功能：
        # Subfunction: Analytical opt_t using only neutral subs
        mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
            tml = (n_neutral+shape-1)/(exp_rel_neutral+(1/scale))
            if (shape<=1) { # i.e. when theta<=1
                tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1, since otherwise it takes meaningless values
            }
            return(tml)
        }
        ## 子功能
        # Subfunction: dNdScv per gene
        selfun_cv = function(j) {
            y = as.numeric(genemuts[j,-1])
            x = RefCDS[[j]]
            exp_rel = y[5:8]/y[5]
            ## gama
            # Gamma
            shape = theta
            scale = y[9]/theta
            ## 中性模型
            # a. Neutral model
            indneut = 1:4 # vector of neutral mutation types under this model (1=synonymous, 2=missense, 3=nonsense, 4=essential_splice)
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/y[5]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik null model
            ## 错义突变模型
            # b. Missense model: wmis==1, free wnon, free wspl
            indneut = 1:2
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
            llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
            ## 截断突变模型 
            # c. Truncating muts model: free wmis, wnon==wspl==1
            indneut = c(1,3,4)
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[2]/y[6]/mrfold; wfree[y[2]==0] = 0
            lltrunc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
            ## 自由选择模型
            # d. Free selection model: free wmis, free wnon, free wspl
            indneut = 1
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[2:4]/y[6:8]/mrfold; wfree[y[2:4]==0] = 0
            llall_unc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
          
            if (constrain_wnon_wspl == 0) {
            
                p = 1-pchisq(2*(llall_unc-c(llmis,lltrunc,ll0)),df=c(1,2,3))
                return(c(wfree,p))
                
            } else { # d2. Free selection model: free wmis, free wnon==wspl
            
                wmisfree = y[2]/y[6]/mrfold; wmisfree[y[2]==0] = 0
                wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
                llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wmisfree,wtruncfree,wtruncfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis, free wnon==wspl
                p = 1-pchisq(2*c(llall_unc-llmis,llall-c(lltrunc,ll0)),df=c(1,1,2))
                return(c(wmisfree,wtruncfree,wtruncfree,p))
            }
        }
        
        sel_cv = as.data.frame(t(sapply(1:nrow(genemuts), selfun_cv)))
        colnames(sel_cv) = c("wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv")
        sel_cv$qmis_cv = p.adjust(sel_cv$pmis_cv, method="BH")
        sel_cv$qtrunc_cv = p.adjust(sel_cv$ptrunc_cv, method="BH")
        sel_cv$qallsubs_cv = p.adjust(sel_cv$pallsubs_cv, method="BH")
        sel_cv = cbind(genemuts[,1:5],sel_cv)
        sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
        ## indel 频发：基于负二项分布回归
        ## Indel recurrence: based on a negative binomial regression (ideally fitted excluding major known driver genes)
        
        if (nrow(indels) >= min_indels) {
        
            geneindels = as.data.frame(array(0,dim=c(length(RefCDS),8)))
            colnames(geneindels) = c("gene_name","n_ind","n_induniq","n_indused","cds_length","excl","exp_unif","exp_indcv")
            geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
            geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,2]),2] = 0
            geneindels$n_induniq = as.numeric(table(unique(indels[,-1])$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,3]),3] = 0
            geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
            
            if (use_indel_sites) {
                geneindels$n_indused = geneindels[,3]
            } else {
                geneindels$n_indused = geneindels[,2]
            }
            ## 剔除 已知的癌症基因
            # Excluding known cancer genes (first using the input or the default list, but if this fails, we use a shorter data-driven list)
            geneindels$excl = (geneindels[,1] %in% known_cancergenes)
            min_bkg_genes = 50
            if (sum(!geneindels$excl)<min_bkg_genes | sum(geneindels[!geneindels$excl,"n_indused"]) == 0) { # If the exclusion list is too restrictive (e.g. targeted sequencing of cancer genes), then identify a shorter list of selected genes using the substitutions.
                newkc = as.vector(sel_cv$gene_name[sel_cv$qallsubs_cv<0.01])
                geneindels$excl = (geneindels[,1] %in% newkc)
                if (sum(!geneindels$excl)<min_bkg_genes | sum(geneindels[!geneindels$excl,"n_indused"]) == 0) { # If the new exclusion list is still too restrictive, then do not apply a restriction.
                    geneindels$excl = F
                    message("    No gene was excluded from the background indel model.")
                } else {
                    warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", paste(newkc, collapse=", ")))
                }
            }
            
            geneindels$exp_unif = sum(geneindels[!geneindels$excl,"n_indused"]) / sum(geneindels[!geneindels$excl,"cds_length"]) * geneindels$cds_length
            ## 对indels进行负二项回归
            # Negative binomial regression for indels
            
            if (is.null(cv)) {
                nbrdf = geneindels[,c("n_indused","exp_unif")][!geneindels[,6],] # We exclude known drivers from the fit
                model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
                nbrdf_all = geneindels[,c("n_indused","exp_unif")]
            } else {
                nbrdf = cbind(geneindels[,c("n_indused","exp_unif")], covs)[!geneindels[,6],] # We exclude known drivers from the fit
                if (sum(!geneindels$excl)<500) { # If there are <500 genes, we run the regression without covariates
                    model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
                } else {
                    model = tryCatch({
                        MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf) # We try running the model with covariates
                    }, warning = function(w){
                        MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf) # If there are warnings or errors we run the model without covariates
                    }, error = function(e){
                        MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf) # If there are warnings or errors we run the model without covariates
                    })
                }
                nbrdf_all = cbind(geneindels[,c("n_indused","exp_unif")], covs)
            }
            message(sprintf("    Regression model for indels (theta = %0.3g)", model$theta))
            
            theta_indels = model$theta
            nbregind = model
            geneindels$exp_indcv = exp(predict(model,nbrdf_all))
            geneindels$wind = geneindels$n_indused / geneindels$exp_indcv
            ## 对每个基因的indel频发进行统计学检验
            # Statistical testing for indel recurrence per gene
        
            geneindels$pind = pnbinom(q=geneindels$n_indused-1, mu=geneindels$exp_indcv, size=theta_indels, lower.tail=F)
            geneindels$qind = p.adjust(geneindels$pind, method="BH")
            ## fisher 联合 p-values 
            # Fisher combined p-values (substitutions and indels)
            
            sel_cv = merge(sel_cv, geneindels, by="gene_name")[,c("gene_name","n_syn","n_mis","n_non","n_spl","n_indused","wmis_cv","wnon_cv","wspl_cv","wind","pmis_cv","ptrunc_cv","pallsubs_cv","pind","qmis_cv","qtrunc_cv","qallsubs_cv")]
            colnames(sel_cv) = c("gene_name","n_syn","n_mis","n_non","n_spl","n_ind","wmis_cv","wnon_cv","wspl_cv","wind_cv","pmis_cv","ptrunc_cv","pallsubs_cv","pind_cv","qmis_cv","qtrunc_cv","qallsubs_cv")
            sel_cv$pglobal_cv = 1 - pchisq(-2 * (log(sel_cv$pallsubs_cv) + log(sel_cv$pind_cv)), df = 4)
            sel_cv$qglobal_cv = p.adjust(sel_cv$pglobal, method="BH")
        
            sel_cv = sel_cv[order(sel_cv$pglobal_cv, sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
        }
    }
    
    if (!any(!is.na(wrong_ref))) {
        wrong_refbase = NULL # Output value if there were no wrong bases
    }
    
    annot = annot[,setdiff(colnames(annot),c("start","end","geneind"))]
    
    if (outmats) {
        dndscvout = list(globaldnds = globaldnds, sel_cv = sel_cv, sel_loc = sel_loc, annotmuts = annot, genemuts = genemuts, mle_submodel = mle_submodel, exclsamples = exclsamples, exclmuts = exclmuts, nbreg = nbreg, nbregind = nbregind, poissmodel = poissmodel, wrongmuts = wrong_refbase, N = Nall, L = Lall)
    } else {
        dndscvout = list(globaldnds = globaldnds, sel_cv = sel_cv, sel_loc = sel_loc, annotmuts = annot, genemuts = genemuts, mle_submodel = mle_submodel, exclsamples = exclsamples, exclmuts = exclmuts, nbreg = nbreg, nbregind = nbregind, poissmodel = poissmodel, wrongmuts = wrong_refbase)
    }
    
} # EOF
