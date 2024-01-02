# ==================================================================
# you can create a Seurat object by using the data from single-cell portal
# 
load('~/work/eoe/result/eoe.rda')

x = melt(table(eoe@meta.data$orig.ident, eoe@active.ident))
colnames(x) = c('patient', 'cell_type', 'count')

patient = eoe@meta.data$orig.ident
disease = eoe@meta.data$disease

id = !duplicated(patient)
disease = disease[id]
names(disease) = patient[id]

# =====
ver = eoe@meta.data$version
ver = ver[id]
names(ver) = patient[id]

# =====
loc = eoe@meta.data$loc
loc = loc[id]
names(loc) = patient[id]


x = data.frame(x,
               loc=loc[as.character(x$patient)], 
               disease=disease[as.character(x$patient)],
               version=ver[as.character(x$patient)])


detach(package:plyr)

x = x %>% dplyr::group_by(patient) %>% 
  dplyr::mutate(total=sum(count))

x$disease = factor(x$disease, levels = c('Ctrl', 'Remission', 'Active'))



# ======
meta = read_xlsx('~/work/eoe/supp_data/supp_data1_EoE_metadata.xlsx') %>% 
  as.data.frame()

rownames(meta) = meta[, 1]

x = data.frame(x, steroid=meta[ident[as.character(x$patient)], 'Steroid?'])


out = setNames(vector('list', length(unique(x$cell_type))), unique(x$cell_type))

for (z in levels(x$cell_type)) {
  model.nb = NULL
  
  cat(z, sep='\n')
  
  id = x$cell_type == z
  xx = x[id, ]
  
  if (sum(xx$count) == 0) {
    next 
  }
  
  model.nb = glm.nb(count ~  version + disease + loc + steroid + offset(log(total)), 
                    data=xx, maxit=1000)

  res = coef(summary(model.nb))
  
  out[[z]] = res
}


pvalue.active = sapply(out, function(z) z['diseaseActive', 4])

pvalue.remission = sapply(out, function(z) z['diseaseRemission', 4])

