# ==================================================================
# you can create a Seurat object by using the data from single-cell portal
# 

rm(list = ls())

load('~/work/eoe/result/eoe.rda')
load('~/work/eoe/data/net.rda')

gene.all = unique(c(net$`Ligand.1`, net$`Ligand.2`,
                    net$`Receptor.1`, net$`Receptor.2`, net$`Receptor.3`))

gene = intersect(gene.all, rownames(eoe@assays$RNA@counts))


# =====
db.pair = list()

for (i in seq(nrow(net))) {
  gene.1 = as.character(net[i, 1:2])
  gene.1 = gene.1[!is.na(gene.1)]
  
  gene.1 = intersect(gene.1, gene)
  
  # =====
  gene.2 = as.character(net[i, 3:5])
  gene.2 = gene.2[!is.na(gene.2)]
  
  gene.2 = intersect(gene.2, gene)
  
  db.pair[[i]] = list(gene.1, gene.2)
}


# =====
id = eoe@meta.data$disease == 'Active'

x = eoe@assays$RNA@counts[gene, id]

cell.type = as.character(eoe@active.ident)[id]


# ==============================================================================
ligand.receptor.count = list()

cell.type.lev = levels(eoe@active.ident)

for (cell.type.i in cell.type.lev) {
  cat(cell.type.i, '\n')
  
  id = which(cell.type.i == cell.type)
  
  y = x[, id, drop=FALSE]
  
  out = sapply(seq(length(db.pair)), function(z) {
    gene = db.pair[[z]][[1]]

    if (length(gene) <= 1) {
      min.count = 1.0
      local.norm = (sum(y[gene, ] > 0) + min.count) / ncol(y)
      local.norm = min(local.norm, 1.0)
      
      res = sum(y[gene, , drop=FALSE]) / sum(x[gene, ]) * local.norm 
    } else {
      res = sapply(gene, function(zz) {
        min.count = 1.0
        local.norm = (sum(y[zz, ] > 0) + min.count) / ncol(y)
        local.norm = min(max(local.norm, 0.25), 1.0)
        
        sum(y[zz, , drop=FALSE]) / sum(x[zz, ]) * local.norm
      })
      res = (prod(res))^(1/length(res)) * (length(res))
    }
    res[is.na(res)] = 0
    
    res
  })
  
  out1 = sapply(seq(length(db.pair)), function(z) {
    gene = db.pair[[z]][[2]]
    
    if (length(gene) <= 1) {
      min.count = 1.0 
      local.norm = (sum(y[gene, ] > 0) + min.count) / ncol(y)
      local.norm = min(local.norm, 1.0)
      
      res = sum(y[gene, , drop=FALSE]) / sum(x[gene, ]) * local.norm 
    } else {
      res = sapply(gene, function(zz) {
        min.count = 1.0
        local.norm = (sum(y[zz, ] > 0) + min.count) / ncol(y)
        local.norm = min(max(local.norm, 0.25), 1.0)
        
        sum(y[zz, , drop=FALSE]) / sum(x[zz, ]) * local.norm
      })
      res = (prod(res))^(1/length(res)) * (length(res))
    }
    res[is.na(res)] = 0
    
    res
  })
  
  ligand.receptor.count[[cell.type.i]] = list(out, out1)
}



# ==============================================================================
net.add = net %>% mutate(name = paste(.[[1]], .[[2]], .[[3]], .[[4]], .[[5]], 
                                        sep=':'))
names(db.pair) = as.data.frame(net.add$name)[,1]

out = sapply(cell.type.lev, function(z1) {
  res = sapply(cell.type.lev, function(z2) {
    x2 = ligand.receptor.count[[z2]][[2]] * ligand.receptor.count[[z1]][[1]] +
      ligand.receptor.count[[z2]][[1]] * ligand.receptor.count[[z1]][[2]]
    
    names(x2) = seq(length(db.pair))
    
    x2 = sort(x2, decreasing = TRUE)
    
    # 
    num.int = 10
    id = as.numeric(names(x2))
    
    ligand.all = vector('list', num.int)
    recepter.all = vector('list', num.int)
    
    # 
    GetLigandReceptor = function(idx) {
      ligand = net[idx, 1:2]
      ligand = as.character(ligand)
      ligand = setdiff(ligand, '')
      ligand = unique(ligand)
      
      receptor = net[idx, 3:5]
      receptor = as.character(receptor)
      receptor = setdiff(receptor, '')
      receptor = unique(receptor)
      
      list(ligand=ligand, receptor=receptor)
    }
    
    idx = 1
    lig.rec = GetLigandReceptor(id[idx])
    
    ligand.all[[1]] = lig.rec$ligand
    recepter.all[[1]] = lig.rec$receptor
    
    # 
    CheckOverlap = function(lig.rec, gene.list) {
      sat = sapply(gene.list, function(z) {
        ll = length(z)
        
        if (ll <= 0) {
          return (FALSE)
        }
        
        lig = lig.rec$ligand
        rec = lig.rec$receptor
        
        if (sum(lig %in% z) == length(lig)) {
          return (TRUE)
        }
        
        if (sum(rec %in% z) == length(rec)) {
          return (TRUE)
        }
        
        if (sum(z %in% lig) == ll) {
          return (TRUE)
        }
        
        if (sum(z %in% rec) == ll) {
          return (TRUE)
        }
        
        return (FALSE)
      })
      
      if (sum(sat) > 0) {
        return (TRUE)
      }
      
      return (FALSE)
    }
    
    idx.all = rep(1, num.int)
    
    for (i in 2:length(id)) {
      
      lig.rec = GetLigandReceptor(id[i])
      
      sat = CheckOverlap(lig.rec = lig.rec, gene.list = ligand.all)
      if (sat) {
        next
      }
      sat = CheckOverlap(lig.rec = lig.rec, gene.list = recepter.all)
      if (sat) {
        next
      }
      
      idx = idx + 1
      ligand.all[[idx]] = lig.rec$ligand
      recepter.all[[idx]] = lig.rec$receptor
      
      idx.all[idx] = i
      
      if (idx >= num.int) {
        break
      }
    }
    
    
    net[as.numeric(names(x2)[idx.all]), ]
    
    x2 = sum(x2[idx.all])
    
    x2
  })
  
  th = sort(res, decreasing = TRUE)[4]
  res[res <= th] = 0
  res
}, simplify = TRUE)

num.int = rowSums(out > 0)

aa = range(num.int)
aa = c(aa[1], aa[1]+round(diff(aa)/2), aa[2])
col_fun_int = circlize::colorRamp2(aa, c(distinct.col[2], 'white', distinct.col[1]))

id = rowSums(out) > 0
id.col = colSums(out) > 0

library(ComplexHeatmap)

row_ha = rowAnnotation(
  '# interactions'=num.int[id], 
  col = list(# 'Sum of rows' = col_fun_row, 
    '# interactions' = col_fun_int
  ))


aa = out[id, id.col]

pheatmap(aa, scale='column', cluster_cols = TRUE, cluster_rows=TRUE, 
        name='Column (Z-score)', left_annotation = row_ha, 
        clustering_method='ward.D2')

