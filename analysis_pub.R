#####
#data step
#####
species = read.csv("species.csv", row.names = 1)
tax = read.csv("tax.csv", row.names = 1)
met = read.csv("met.csv", row.names = 1)
library(phyloseq)
pd = merge_phyloseq(otu_table(as.matrix(species), taxa_are_rows = T), sample_data(met), tax_table(as.matrix(tax)))

#####
#differential abundance testing
#####
library(ANCOMBC)

out = ancombc(phyloseq = pd, formula = "group", lib_cut = 1)
res = out$res
sig = c()
for(i in 1:4){
  sig = c(sig, row.names(res$diff_abn)[which(res$diff_abn[,i] == TRUE)])
  
}
con_ref = res$q_val[sig,]

ord = unique(met$group)
sample_data(pd)$group = ordered(sample_data(pd)$group, levels = ord[c(1,2,3,4,5)])
out = ancombc(phyloseq = pd, formula = "group", lib_cut = 1)
res = out$res
sig = c()
for(i in 1:3){
  sig = c(sig, row.names(res$diff_abn)[which(res$diff_abn[,i] == TRUE)])
  
}
ec_ref = res$q_val[sig,1:3]
names(ec_ref)[1:3] = ord[c(2,3,4)]

sample_data(pd)$group = ordered(sample_data(pd)$group, levels = ord[c(2,3,4,5,1)])
out = ancombc(phyloseq = pd, formula = "group", lib_cut = 1)
res = out$res
sig = c()
for(i in 1:2){
  sig = c(sig, row.names(res$diff_abn)[which(res$diff_abn[,i] == TRUE)])
  
}
ei_ref = res$q_val[sig,1:2]
names(ei_ref)[1:2] = ord[c(3,4)]

sample_data(pd)$group = ordered(sample_data(pd)$group, levels = ord[c(3,4,5,1,2)])
out = ancombc(phyloseq = pd, formula = "group", lib_cut = 1)
res = out$res
sig = row.names(res$diff_abn)[which(res$diff_abn[,1] == TRUE)]
wc_ref = res$q_val[sig,1,drop = FALSE]
names(wc_ref)[1] = ord[4]

con_ref$z = row.names(con_ref)
con_ref = melt(con_ref, id.vars = "z")
con_ref$variable = paste("Control_v_", con_ref$variable, sep = "")
con_ref = con_ref[con_ref$value < 0.05,]

ec_ref$z = row.names(ec_ref)
ec_ref = melt(ec_ref, id.vars = "z")
ec_ref$variable = paste("Egg+Cellulose_v_", ec_ref$variable, sep = "")
ec_ref = ec_ref[ec_ref$value < 0.05,]

ei_ref$z = row.names(ei_ref)
ei_ref = melt(ei_ref, id.vars = "z")
ei_ref$variable = paste("Egg+Inulin_v_", ei_ref$variable, sep = "")
ei_ref = ei_ref[ei_ref$value < 0.05,]

wc_ref$z = row.names(wc_ref)
wc_ref = melt(wc_ref, id.vars = "z")
wc_ref$variable = paste("Whey+Cellulose_v_", wc_ref$variable, sep = "")

pw = rbind(con_ref, ec_ref, ei_ref, wc_ref)
pw$z = sub("\\..*", "", pw$z)
pw$species = tax$Species[match(pw$z, row.names(tax))]

hold = pw[!duplicated(pw),]

#to change the p-value cutoff for plotting
# pw = pw[pw$value < 0.05,]

color = colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(25)

temp = data.frame(otu_table(pd))
temp = data.frame(t(temp))
temp = temp/rowSums(temp)
temp = temp[,names(temp) %in% pw$z]
pw$species = sub("_subsp.*", "", pw$species)
newnames <- lapply(pw$species[match(names(temp), pw$z)],
  function(x) bquote(italic(.(x))))

met$ab = met$sample
met$ab = ifelse(met$group == "Control", paste("CON", met$ab, sep = ""), met$ab)
met$ab = ifelse(met$group == "Egg+Cellulose", paste("EC", met$ab, sep = ""), met$ab)
met$ab = ifelse(met$group == "Egg+Inulin", paste("EI", met$ab, sep = ""), met$ab)
met$ab = ifelse(met$group == "Whey+Cellulose", paste("WC", met$ab, sep = ""), met$ab)
met$ab = ifelse(met$group == "Whey+Inulin", paste("WI", met$ab, sep = ""), met$ab)

row.names(temp) = met$ab[match(row.names(temp), row.names(met))]

library(pheatmap)
# png('heatmap.png', res = 600, width = 9, height = 18, units = 'in')
# print(pheatmap(as.matrix(t(temp)), color = color,  labels_row = as.expression(newnames)), fontsize_col = 12, fontsize_row = 12)
# dev.off()

temp2 = temp[,which(colMeans(temp) > 0.001)]
newnames <- lapply(pw$species[match(names(temp2), pw$z)],
                   function(x) bquote(italic(.(x))))

# png('heatmap2.png', res = 600, width = 9, height = 10, units = 'in')
# print(pheatmap(as.matrix(t(temp2)), color = color,  labels_row = as.expression(newnames)),fontsize_col = 12, fontsize_row = 12)
# dev.off()


con = colMeans(temp[grep("CON", row.names(temp)),])
con = con * 100
pw$'Control' = con[match(pw$z, names(con))]

ec = colMeans(temp[grep("EC", row.names(temp)),])
ec = ec * 100
pw$'Egg+Cellulose' = ec[match(pw$z, names(ec))]

ei = colMeans(temp[grep("EI", row.names(temp)),])
ei = ei * 100
pw$'Egg+Inulin' = ei[match(pw$z, names(ei))]

wc = colMeans(temp[grep("WC", row.names(temp)),])
wc = wc * 100
pw$'Whey+Cellulose' = wc[match(pw$z, names(wc))]

wi = colMeans(temp[grep("WI", row.names(temp)),])
wi = wi * 100
pw$'Whey+Inulin' = wi[match(pw$z, names(wi))]

pw = pw[,-1]
names(pw)[1:2] = c("comparison", "q_value")

#the percentage of diet groups microbiome comprised of DA species
da = data.frame(sig_sum_abun = rowSums(temp))
da$group = met$group[match(row.names(da), met$ab)]
aggregate(da$sig_sum_abun, list(da$group), mean)

temp = pw[!duplicated(pw),]

t2 = tax
t2$Species = sub("_subsp.*", "", t2$Species)

temp$Genus = t2$Genus[match(temp$species, t2$Species)]
temp$Family = t2$Family[match(temp$species, t2$Species)]
temp$Order = t2$Order[match(temp$species, t2$Species)]
temp$Class = t2$Class[match(temp$species, t2$Species)]
temp$Phylum = t2$Phylum[match(temp$species, t2$Species)]
temp = temp[,c(13,12,11,10,9,3,1,2,4,5,6,7,8)]
temp = temp[order(temp$q_value),]

pw = pw[!duplicated(pw),]
data.frame(table(pw$comparison))
dd = data.frame(table(pw$comparison))
dd = dd[order(dd$Freq, decreasing = TRUE),]
dd = dd[grep("Control", dd$Var1),]

tt = data.frame(t(otu_table(pd)))
tt[tt > 0] = 1
tt$group = met$group[match(row.names(tt), row.names(met))]
spc_pr_group = c()
for(i in 1:length(unique(tt$group))){
  temp = tt[tt$group == unique(tt$group)[i],-ncol(tt)]
  temp = data.frame(t(colSums(temp)))
  row.names(temp) = unique(tt$group)[i]
  temp[temp > 0] = 1
  spc_pr_group = rbind(spc_pr_group, temp)
}
rowSums(spc_pr_group)
dd$total_species = c(154, 218, 258, 272)
dd$Distance = dd$Freq/dd$total_species
dd$Var1 = sub(".*group", "", dd$Var1)
dd$Var1 = ordered(dd$Var1, levels = dd$Var1)
p = ggplot(dd, aes(Var1, Distance)) + geom_point(size = 5, color = "grey60") + geom_point(size = 2) + ylim(0,1) + labs(x = "", y = "Distance from Control\n(Proportion of Species DA)") + theme_bw(base_size = 15)
# ggsave("DA_from_control.png", p, height = 4, width = 6, dpi = 300, units = "in")

#####
#barplot
#####
spec = data.frame(otu_table(pd))
s = c(tax_table(pd)[,7])
row.names(spec) = s
spec = data.frame(t(spec))
spec = spec/rowSums(spec)

who = names(sort(colMeans(spec), decreasing = TRUE))[1:19]
f = spec[,names(spec) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = cbind(met, f)
m = melt(dd, id.vars = c("sample", "ab", "group"), measure.vars = who)
col = scale_fill_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(19 - 11), brewer_pal(palette = "Spectral")(11))),"#CCCCCC"), name = "Species")
cbb = scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "thistle4", "yellowgreen", "violetred1"))

p = ggplot(m, aes(ab, fill = variable)) + geom_bar(aes(weight = value), position = position_stack(reverse = TRUE))  +
  theme_bw(base_size = 15) + facet_wrap(~ group, scales = "free") + 
  col  + theme(axis.text.x = element_blank()) +
  xlab("Sample") + ylab("Relative Abundance") +  theme(legend.justification=c(1,0), legend.position=c(1.008,0)) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(reverse=TRUE,ncol = 2)) + scale_y_continuous(labels = percent_format())
# ggsave("bar.png", p, width = 16.75, height = 8, units = "in", dpi = 600)
p2 = p
#####
#beta
#####
v = vegdist(t(otu_table(pd)))

a = adonis(v ~ group, data = met)
a = adonis(v ~ protein * fiber, data = met) #note the summed R2 is the sames as just for group, but here you can talk about interactions
holda = a
# # write.csv(data.frame(a[[1]]), "adonis.csv")
#####
#read in metabolomics
#####
chem = read.csv("variable_measures_exp19.csv", row.names = 1)
chem = chem[,-1]
look = c("560_576", "583_581")
temp = chem[row.names(chem) %in% look,] #no clue why 'c' is making data frames
temp = rbind(temp, temp)
row.names(temp) = c("560","583", "576", "581")
chem = chem[!row.names(chem) %in% look,] 
chem = rbind(chem, temp)
row.names(chem) = paste("X", row.names(chem), sep = "")
rm(look, temp)

#####
#correlate metabolomics with species
#####

t1 = data.frame(t(chem))
t2 = data.frame(t(data.frame(otu_table(pd))))
t2 = t2/rowSums(t2)
t2 = data.frame(t(t2))
t2 = t2[,names(t2) %in% names(t1)]
keep = c()
for(i in unique(met$group)){
  hold = row.names(met)[met$group == i]
  hold = t2[,names(t2) %in% hold]
  hold = row.names(t2)[which(rowMeans(hold) >= 0.02)]
  keep = c(keep, hold)
}
keep = unique(keep)
t2 = t2[row.names(t2) %in% keep,]
t1 = rbind(t1, t2)
t1 = t(t1)
t1 = cor(t1, method = "spearman")
t1 = t1[,colnames(t1) %in% row.names(t2)]
t1 = t1[!rownames(t1) %in% row.names(t2),]
temp = data.frame(tax_table(pd))
annotation_col = data.frame(temp[,2][match(colnames(t1), row.names(temp))])
colnames(t1) = temp[,7][match(colnames(t1), row.names(temp))]
m = melt(t1)

m = m[order(m$value, decreasing = TRUE),]
length(which(m$value > 0.25))
library(pheatmap)

colnames(t1) = sub("_", " ", colnames(t1))
newnames <- lapply(colnames(t1),function(x) bquote(italic(.(x))))

row.names(annotation_col) = colnames(t1)
colnames(annotation_col) = "Phylum"

cclass = read.csv("chemical_classes_v2.csv")
row.names(cclass) = cclass$chemical
cclass = cclass[,1,drop = FALSE]

# png('species_metabalomics_heatmap_phylum.png', res = 600, width = 10, height = 17, units = 'in')
# print(pheatmap((t1),fontsize_col = 10, fontsize_row = 10, display_numbers = FALSE, labels_col  = as.expression(newnames), annotation_col = annotation_col, annotation_row = cclass))
# dev.off()

#####
#read in metagenomic function
#####
f = read.delim("pathway_prediction.txt", row.names = 1)
dc = f[,(ncol(f) - 2):ncol(f)]
f = f[, -((ncol(f) - 2):ncol(f))]
k = read.delim("key", header = F)
k[,1] = paste("X", k[,1], sep = "")
names(f) = k[,1][match(names(f), k[,2])]
f = f[,names(f) %in% row.names(chem)]
f = f[rowSums(f) > 0,]
f = data.frame(t(f))
all(row.names(f) %in% row.names(chem))
all(row.names(chem) %in% row.names(f))

#####
#test group effect on function
#####
temp = vegdist(f)
t2 = met[row.names(met) %in% row.names(f),]
a = adonis(temp ~ protein * fiber, dat = t2)
1 - 0.66869
# write.csv(data.frame(a[[1]]), "function_adonis.csv")

both = data.frame(set = c(rep("Composition", 4), rep("Function", 4)), variance = rep(row.names(a[[1]])[1:4],2), r2 = c(holda[[1]][1:4,5], a[[1]][1:4,5]))
both$variance = sub("protein:fiber", "interaction", both$variance)
both$variance = sub("Residuals", "residual", both$variance)

p = ggplot(both, aes(set,fill = variance) )+ geom_bar(aes(weight = r2), position = "fill") +
  theme_bw(base_size = 15)  + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("") + ylab("Relative Abundance") + scale_y_continuous(labels = percent_format())
p$data$variance = factor(p$data$variance, ordered = TRUE, levels = c("protein", "fiber", "interaction", "residual"))
p + guides(fill = guide_legend(reverse=TRUE))

p = ggplot(both, aes(set, fill = variance)) + geom_bar(aes(weight = r2), position = position_stack(reverse = TRUE))  +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("") + ylab("Total Variance") +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(reverse=TRUE)) + scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("orchid4", "#E69F00", "#56B4E9", "#999999"),
                    name="Variance\nExplained",
                    breaks=c("protein", "fiber", "interaction", "residual"))
p$data$variance = factor(p$data$variance, ordered = TRUE, levels = c("protein", "fiber", "interaction", "residual"))
# ggsave("adonis_bar.png", p, width = 5, height = 6, units = "in", dpi = 600)


#####
#correlate metabolomics and function
#####

t1 = data.frame(t(chem))
t2 = data.frame(t(f))
t1 = rbind(t1, t2)
t1 = t(t1)
t2 = t1
t1 = cor(t1, method = "spearman")
t1 = t1[,colnames(t1) %in% names(f)]
t1 = t1[!rownames(t1) %in% names(f),]
m = melt(t1)

m = m[order(m$value, decreasing = TRUE),]
length(which(m$value > 0.25))

dc = dc[row.names(dc) %in% colnames(t1),]
dc$level2 = ordered(dc$level2 , levels = unique(dc$level2))
annotation_row = data.frame(dc$level2)
row.names(annotation_row) = row.names(dc)
# png('noclustmbfunction_metabalomics_heatmap.png', res = 600, width = 15, height = 22, units = 'in')
# print(pheatmap(t(t1), #color=colorRampPalette(c("steelblue3", "white", "red"))(50),
#                fontsize_col = 10, fontsize_row = 10, display_numbers = FALSE, annotation_row = annotation_row, cluster_rows = FALSE))#,
# #treeheight_row = 0, treeheight_col = 0, fontsize_number = 12))
# dev.off()

# png('clustmbfunction_metabalomics_heatmap.png', res = 600, width = 15, height = 22, units = 'in')
# print(pheatmap(t(t1), #color=colorRampPalette(c("steelblue3", "white", "red"))(50), fontsize_col = 10, fontsize_row = 10, display_numbers = FALSE, annotation_row = annotation_row))#,
# treeheight_row = 0, treeheight_col = 0, fontsize_number = 12))
# dev.off()

dc$level3 = ordered(dc$level3 , levels = unique(dc$level3))
annotation_row = data.frame(dc$level3)
row.names(annotation_row) = row.names(dc)
# png('l3_clustmbfunction_metabalomics_heatmap.png', res = 600, width = 15, height = 22, units = 'in')
# print(pheatmap(t(t1), fontsize_col = 10, fontsize_row = 10, display_numbers = FALSE, annotation_row = annotation_row))
# dev.off()

#####
#cluster the explanatory variables to make composite variables
#####

pom = function(var, add.one = FALSE){
  var = var - min(var, na.rm = TRUE)
  var = var/max(var, na.rm = TRUE)
  if(add.one == TRUE){
    var = var + 1
  }
  return(var)
}
for(i in 1:ncol(f)){
  f[,i] = pom(f[,i])
}

#positive
b = seq(0.1,.69,0.01)
pc = vector("list", length(b))
for(i in 1:length(b)){
  keep = as.character(m$Var2[which((m$value) >= b[i])])
  f_sub = f[,names(f) %in% keep]
  pc[[i]] = capscale(f_sub ~ 1, distance = "bray")$CA$u
  print(i)
}
for(i in 1:60){print(dim(pc[[i]]))}
possig = c()
for(i in 1:length(pc)){
  temp = data.frame(threshold = b[i], p.value =  anova(capscale(chem ~ ., data.frame(pc[[i]]), distance = "bray"))$`Pr(>F)`[1])
  possig = rbind(possig, temp)
  print(i)
}

#negative
pc = vector("list", length(b))
for(i in 1:length(b)){
  keep = as.character(m$Var2[which((m$value) <= -b[i])])
  if(!identical(keep, character(0))){
    f_sub = f[,names(f) %in% keep]
    pc[[i]] = capscale(f_sub ~ 1, distance = "bray")$CA$u
  }
}

drop = c()
for(i in 1:length(pc)){
  drop = c(drop, !is.null(pc[[i]]))
}
pc = pc[drop]
negsig = c()
for(i in 1:length(pc)){
  temp = data.frame(threshold = -b[i], p.value =  anova(capscale(chem ~ ., data.frame(pc[[i]]), distance = "bray"))$`Pr(>F)`[1])
  negsig = rbind(negsig, temp)
  print(i)
}

possig$set = "positive" 
negsig$set = "negative"
all = rbind(possig, negsig)

p = ggplot(all, aes(threshold, p.value)) + geom_point(size = 3, alpha = 0.4) + theme_bw(base_size = 15) +
  labs(x = "Spearman Correlation Threshold", y = "p-value") + geom_hline(yintercept = 0.05, color = "grey60", linetype = "dashed") 
# ggsave("rda_subsets.png", p, height = 4, width = 6, dpi = 300, units = "in")
pxx = p
keep = unique(m$Var2[which(m$value >= 0.61)])
f_sub = f[,names(f) %in% keep]
x = capscale(f_sub ~ 1, distance = "bray")
dbrda = capscale(chem ~ ., data.frame(x$CA$u), distance = "bray")
anova(dbrda)
summary(dbrda)

comb = c()
for(i in 1:length(keep)){
  for(j in 1:ncol(chem)){
    temp = cor.test(f[,keep[i]], chem[,j], use = "pairwise.complete.obs", method = "spearman")
    temp = data.frame(kegg = keep[i], chemical = names(chem[j]), cor = temp$estimate, pval = temp$p.value)
    comb = rbind(comb, temp)
  }
  print(i)
}
comb = comb[comb$pval <= 0.05,]
comb$level1 = dc$level1[match(comb$kegg, row.names(dc))]
comb$level2 = dc$level2[match(comb$kegg, row.names(dc))]
comb$level3 = dc$level3[match(comb$kegg, row.names(dc))]
comb = comb[order(comb$chemical),]
row.names(comb) = NULL
comb = comb[,c(3,4,2,1,5,6,7)]

hope = dcast(formula = chemical ~ level1, data = comb, value.var = "cor")
hope[is.na(hope)] = 0
row.names(hope) = hope[,1]
hope = hope[,-1]
# png('sig_metabalomics_mbfunction_heatmap.png', res = 600, width = 12, height = 7.5, units = 'in')
# print(pheatmap(t(hope), treeheight_row = 0, treeheight_col = 0, fontsize_number = 12))
# dev.off()

# sig = dc$level1[row.names(dc) %in% names(f_sub)]
# comb$included.in.final.dbrda = ifelse(comb$level1 %in% sig, "YES", "NO")

# write.csv(comb, "sig_chem_mbfunction_cors.csv", row.names = FALSE)

test2 = data.frame(table(comb$chemical))
test2$binom = NA
for(i in 1:nrow(test2)){
  test2$binom[i] = binom.test(test2$Freq[i], nrow(comb), p = 1/133, alternative = "g")$p.value
}
test2 = test2[order(test2$binom),] #these are the terms that were over-represented in the list of significant (0.05) correlations
names(test2) = c("chemical", "significant.associations", "binomial.p.value")
row.names(test2) = NULL
# write.csv(test2, "chemical_sig_association_enrichment.csv", row.names = FALSE)

set = test2[1:12,]
set$chemical = c("Butyric acid", "PC36:0AA", "Citric acid", "C5MDC", "C6:1", "Creatine", "Phenylalanine", "Putrescine", "Propionic acid", "Alpha aminoadipic acid", "LYSOC18.0", "Pyruvic Acid")
set$log = -log10(set$binomial.p.value)
set$chemical = ordered(set$chemical, levels = rev(set$chemical))
p = ggplot(set, aes(chemical)) + geom_bar(aes(weight = log), alpha = 0.6) + coord_flip() + theme_bw(base_size = 15) + 
  labs(y = "-log10(p-value)", x = "")
# ggsave("plasma_log10.png", p, height = 4, width = 6, dpi = 300, units = "in")






