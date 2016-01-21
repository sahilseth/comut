## this comut plot is based on complex heatmap


if(FALSE){
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
}


reshape_ann.old <- function(x){
  require(reshape2);require(dplyr);
  
  
  ## --- remove rows which are . !
  xx = filter(x, !exonicfunc %in% c(".", "unknown"))
  
  resolve_multi_mut <- function(){
    
    multi_mut = group_by(xx, gene.refgene, tumor_name) %>%
      summarise(count = n()) %>% filter(count > 1)
    
    for(i in 1:nrow(multi_mut)){
      keys = subset(xx, gene.refgene == multi_mut$gene.refgene[i] &
                      tumor_name == multi_mut$tumor_name[i]) %>%
        select(key) %>% unlist()
      rm_keys = keys[-1]
      ## Keep only the first mutation
      xx = subset(xx, !(tumor_name == multi_mut$tumor_name[i] & key %in% rm_keys))
    }
  }
  
  ## for each tumor, each gene
  
  yy = dcast(xx, gene.refgene ~ tumor_name, 
             value.var = "exonicfunc.refgene")
  
  ## if a sample has more than one is a gene, pick one
  #table(with(yy, paste(gene.refgene, exonicfunc.refgene)))
  
  ## Top variants per gene:
  require(DT)
  require(knitr)
  tp = group_by(yy, gene.refgene) %>% summarise(count = n()) %>% arrange(-count)
  kable(tp)
  
}


reshape_ann <- function(funcwd, 
                        sample_col = "sample_plat",
                        func_col = "exonicfunc.refgene.x",
                        gene_col = "gene.refgene.x",
                        comb_type = function(x) paste(unique(x), collapse = ";")){
  
  
  # get a clean small dataset
  x2 = select_(funcwd, 
               samp = sample_col, 
               func = func_col,
               gene = gene_col)
  
  # reshape
  xwd = reshape2::dcast(x2, gene ~ samp, value.var = "func", fun.aggregate = comb_type)
  
  row.names(xwd) = xwd$gene
  
  xwd
}




# decide how each mutation shows
get_rect_funcs <- function(mat, colors_mut = colors_mut_func(), get_type = function(x) strsplit(x, ";")[[1]]){
  
  col_funcs = lapply(colors_mut, function(col){
    force(col)
    function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col, col = NA, alpha = 0.8))
    }
  })
  
  myfuncs = c(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey82", col = NA, alpha = 0.1))
    },
    AMP = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "white", col = NA))
    },
    DEL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
    },
    col_funcs
  )
  
  # only return background and the types in MAT
  all_type = unique(unlist(lapply(mat, get_type)))
  
  myfuncs[c("background", all_type)]
  
}

# accepts mutsig input, a maf file and does wonders
# in future a VCF may be supported as well

#' Create a comut plot
#'
#' @param x see details below
#' @param title 
#' @param colors_mut 
#' @param get_type 
#' @param ... 
#'
#' @export
#' 
#' @details 
#' x needs to be a data.frame, where first column is gene names or features, with subsequent 
#' columns - one for each sample.
#' 
#' 
#'
#' @examples \dontrun{
#' p <- comut(mat.mutsig, bottom_annotation = ba, show_column_names = TRUE, title = "MutSig p<0.01 in NCG");p
#' 
#' }
comut <- function(x, 
                  title = "coMut",
                  colors_mut = colors_mut_func(),
                  get_type = function(x) strsplit(x, ";")[[1]], ...
){
  
  mat = as.matrix(x[, -1])
  rownames(mat) = as.character(unlist(x[, 1]))
  
  alter_fun_list = get_rect_funcs(mat, colors_mut, get_type)
  all_type = unique(unlist(lapply(mat, get_type)))
  colors_mut = colors_mut[all_type]
  
  ComplexHeatmap::oncoPrint(mat, 
                            # function which will open up types!
                            get_type = function(x) strsplit(x, ";")[[1]],
                            alter_fun_list = alter_fun_list, 
                            col = colors_mut, 
                            column_title = title, remove_empty_columns = FALSE, ...)
}


#' Given a mat, and clinical columns; add them
#'
#' @param mat 
#'
get_bottom_anno <- function(mat, 
                            clins = c("platform", "gender", "final_pat_iaslc", "smoking_type")){
  
  sampinfo = data.frame(
    sample_id = gsub("(.*)\\.tumor_f.*", "\\1", colnames(mat)[-1]), 
    platform = gsub(".*tumor_f_(.*)", "\\1", colnames(mat)[-1]), stringsAsFactors = FALSE) %>% tbl_df()
  
  # all clin info
  clin = readRDS("~/projects2/sarco/20151110_cleanclinical/clin.rds")
  sampinfo = left_join(sampinfo, clin) %>% apply(2, as.factor) %>% data.frame(stringsAsFactors = FALSE)
  
  sampinfo <<- sampinfo
  
  library(RColorBrewer)
  set1_cols = RColorBrewer::brewer.pal(n = 9, "Set1")
  
  # clinical colors
  clin_cols = list(
    platform = c(wex = "red", ccp = "blue", rna = "green"),
    gender = c(Male = "white", Female = "black"),
    final_pat_iaslc = structure(brewer.pal(7, "Blues"), names = c("IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IV")),
    smoking_type = c(NEVER = "white", FORMER = "grey", CURRENT = "black")
  )
  
  clin_cols <- clin_cols[clins]
  
  
  str(sampinfo)
  summary(sampinfo)
  
  ann_hts = rep(unit(5, "mm"), length(clin_cols))
  ba =  dplyr::select(sampinfo, one_of(names(clin_cols))) %>% 
    ch$HeatmapAnnotation(col = clin_cols, 
                         text = ch$anno_text(sampinfo$sample_id, rot = 90, offset = 1, just = 'right'),
                         annotation_height = unit.c(ann_hts, ch$max_text_width(sampinfo$sample_id)))
}



