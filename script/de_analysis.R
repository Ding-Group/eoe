# ==================================================================
# you can create a Seurat object by using the data from single-cell portal
# 
rm(list = ls())

load('~/work/eoe/result/eoe.rda')

# =====
gene.all = rownames(eoe@assays$RNA@counts)

mt.gene  = grep('^MT-', gene.all, value = TRUE)
ribosome.gene = grep('^RPL|^RPS|^MRPS|^MRPL', gene.all, value = TRUE)
problematic.gene = c('MALAT1', 'EEF1A1', 'TPT1')

gene.to.remove = c(mt.gene, ribosome.gene, problematic.gene)

features = setdiff(gene.all, c(mt.gene, ribosome.gene, problematic.gene))

eoe@meta.data[, 'nFeature_RNA_log2'] = log2(eoe@meta.data$nFeature_RNA)

table(eoe@meta.data$patient, eoe@meta.data$disease)

# =====
meta = read_xlsx('~/work/eoe/supp_data/supp_data1_EoE_metadata.xlsx') %>% 
  as.data.frame()

rownames(meta) = meta[, 1]

eoe@meta.data[, 'treat'] = meta[ident[as.character(eoe@meta.data$orig.ident)], 'Steroid?']


cluster.marker.lr.log2 = sapply(cell.type.lev, function(z) {
  cat(z, sep='\n')
  
  ident.1 = (eoe@active.ident  %in% z) & 
    (eoe@meta.data$disease == 'Active')
  ident.1 = names(eoe@active.ident)[ident.1]
  
  ident.2 = (eoe@active.ident  %in% z) & 
    (eoe@meta.data$disease %in% c('Ctrl'))
  ident.2 = names(eoe@active.ident)[ident.2]
  
  
  if (length(ident.1) < 10 | length(ident.2) < 10) {
    return (NULL)
  }
  
  cluster.markers = FindMarkers(eoe, ident.1, ident.2, test.use = 'LR', 
                                latent.vars = c('nFeature_RNA_log2', 'version', 'loc', 'treat'), 
                                features = features)
  
  return (cluster.markers)
}, simplify = FALSE)



FilterAmbient = function(case, gene, cell.type) {
  # 
  gene = intersect(gene, rownames(eoe@assays$RNA@counts))
  
  if (length(gene) <= 0) {
    return(NULL)
  }
  
  id = eoe@meta.data$orig.ident == case
  cell = names(eoe@active.ident)[id]
  cell = intersect(cell, colnames(eoe@assays$RNA@counts))
  
  eoe.sub = subset(eoe, cells = cell)
  
  if (sum(eoe.sub@active.ident == cell.type) <= 10) {
    return (NULL)
  }
  
  res = FindMarkers(eoe.sub, cell.type, test.use = 'LR',  
                    latent.vars = c('nFeature_RNA_log2'), 
                    setdiff(unique(eoe.sub@active.ident), cell.type), 
                    features = gene)
  
  res
}


SelectDEGene = function(cell.type, de.list.up, de.list.down, disease='Active', ctrl='Ctrl') {
  id = eoe@meta.data$disease == disease
  case = unique(eoe@meta.data$orig.ident[id])
  
  gene = de.list.up$gene
  
  out = sapply(case, function(z) {
    cat(z, sep='\n')
    FilterAmbient(z, gene = gene, cell.type = cell.type)
  }, simplify = FALSE)
  
  
  res = sapply(out, function(z) {
    id = z$avg_log2FC < 0 
    neg = rownames(z)[id]
    
    id = z$avg_log2FC > 0 & z$p_val_adj < 0.01 
    pos = rownames(z)[id]
    
    res = list(pos=pos, neg=neg)
    
    res
  }, simplify = FALSE)
  
  pos = table(unlist(sapply(res, function(z) z$pos)))
  neg = table(unlist(sapply(res, function(z) z$neg)))
  
  gene.add = names(pos)
  g = intersect(names(pos), names(neg))
  
  gene.up.filt = setdiff(gene.add, g[which(pos[g] - neg[g] < 0)])
  
  
  # =====
  id = eoe@meta.data$disease %in% ctrl
  case = unique(eoe@meta.data$orig.ident[id])
  
  gene = de.list.down$gene
  
  out = sapply(case, function(z) {
    cat(z, sep='\n')
    FilterAmbient(z, gene = gene, cell.type = cell.type)
  }, simplify = FALSE)
  
  
  res = sapply(out, function(z) {
    id = z$avg_log2FC < 0 
    neg = rownames(z)[id]
    
    id = z$avg_log2FC > 0 & z$p_val_adj < 0.01 
    pos = rownames(z)[id]
    
    res = list(pos=pos, neg=neg)
    
    res
  }, simplify = FALSE)
  
  pos = table(unlist(sapply(res, function(z) z$pos)))
  neg = table(unlist(sapply(res, function(z) z$neg)))
  
  gene.add = names(pos)
  g = intersect(names(pos), names(neg))
  
  gene.down.filt = setdiff(gene.add, g[which(pos[g] - neg[g] < 0)])
  
  res = list(up=gene.up.filt, down=gene.down.filt)
  
  res
}


# =====
marker.up.filt = sapply(cluster.marker.lr.log2, function(z) {
  cat(0, '\n')
  
  if (length(z) <= 0) {
    return(NULL) 
  }
  z = z %>% tibble::rownames_to_column(var='gene') %>% 
    filter(!(gene %in% gene.to.remove)) %>% 
    filter(avg_log2FC > log2(1.5)) %>% 
    filter(pct.1 > 0.2) %>% 
    filter(p_val_adj < 0.01)
}, simplify = FALSE)


marker.down.filt = sapply(cluster.marker.lr.log2, function(z) {
  if (length(z) <= 0) {
    return(NULL) 
  }
  z = z %>% tibble::rownames_to_column(var='gene') %>% 
    filter(!(gene %in% gene.to.remove)) %>% 
    filter(avg_log2FC < -log2(1.5)) %>% 
    filter(pct.2 > 0.2) %>% 
    filter(p_val_adj < 0.01)
}, simplify = FALSE)


cell.type.lev = levels(eoe@active.ident)

marker.filt.ambient = sapply(cell.type.lev, function(z) {
  res = SelectDEGene(z, 
                     de.list.up   = marker.up.filt[[z]], 
                     de.list.down = marker.down.filt[[z]])
  
  res
}, simplify = FALSE)




