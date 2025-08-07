#/home/huangkeyun/hechang/miniconda/envs/r.4.4.1/bin/R --max-ppsize=500000

library(readr)
library(dplyr)
library(mlr3)
library(mlr3filters)
library(mlr3fselect)
library(mlr3learners)
library(mlr3viz)
library(mlr3pipelines)
library(matrixTests)
library(future)

set.seed(123)


output_mat = '/BioII/lulab_b/huangkeyun/hechang/fragmentome/Degradation/matrix/gut/'
ref_genes = read.csv("/BioII/lulab_b/huangkeyun/hechang/fragmentome/Degradation/ref/long_RNA.gencode.bed",header=F,sep='\t')

kl_file_list = list.files("/BioII/lulab_b/baiyilan/hechang/results/gut/") %>% grep("kl",.,value=T)

kl_list = lapply(kl_file_list, function(x) {
    if (file.exists(paste0('/BioII/lulab_b/baiyilan/hechang/results/gut/',x))) {
        one = read_delim(paste0('/BioII/lulab_b/baiyilan/hechang/results/gut/',x))$kl
        return(one)
    } else {
        return(NA)
    }
    return(one)
}) %>% do.call(cbind,.)
samples_id = lapply(kl_file_list,function(x){
  strsplit(x,'\\.')[[1]][1]
})%>% unlist()
colnames(kl_list) = samples_id
rownames(kl_list) = ref_genes[, 4]

#write.csv(kl_list, paste0(output_mat,"/kl_all.csv"), quote = F)

expr = read.csv('/data2/lulab1/exRNA/published/exOmics_RNA/GSE133684/level3.0_Expression/GSE133684_count_matrix_TPM.txt',sep='\t')
metadata = read.csv('/data2/lulab1/exRNA/published/exOmics_RNA/GSE133684/metadata/sample_classes.txt',sep='\t')
kl = kl_list[,rownames(metadata)]
expr = expr[,rownames(metadata)]


feature_select = function(matrix_without_genes,group,nf=1000){
    result = row_oneway_equalvar(matrix_without_genes, group)
    ordered_genes = rownames(result[order(result$pvalue),])[1:nf]
    return(ordered_genes)
}

base_learner_A1 = po("select", id = 'A1',selector = selector_grep("^A_")) %>>%
po("learner", lrn("classif.ranger", predict_type = "prob")) %>% as_learner()

base_learner_B1 = po("select", id = 'B1',selector = selector_grep("^B_")) %>>%
po("learner", lrn("classif.ranger", predict_type = "prob")) %>% as_learner()


# 使用mlr3pipelines的安全Stacking（交叉验证生成元特征）
stacking_graph = ppl(
"stacking",
base_learners = list(base_learner_A1,base_learner_B1), 
super_learner = lrn("classif.glmnet", predict_type = "prob", alpha = 1,lambda = 0.01),  # 元学习器
A1.classif.ranger.resampling.folds = 5,
B1.classif.ranger.resampling.folds = 5
)
ensemble_all = as_learner(stacking_graph)

ML_flow_2 = function(ensemble,matrix1,matrix2, group, id,fold=5){
    genes_selected_A = feature_select(matrix1,group)
    genes_selected_B = feature_select(matrix2,group)

    mat1 = data.frame(t(matrix1[genes_selected_A,]))
    mat2 = data.frame(t(matrix2[genes_selected_B,]))

    print('anova done.')

    colnames(mat1) = paste0("A_", 1:dim(mat1)[2])
    colnames(mat2) = paste0("B_", 1:dim(mat2)[2])

    comb_mat = cbind(mat1, mat2)
    comb_mat$target = group

    task = as_task_classif(comb_mat, target = "target", id = id)
    ensemble$id = id

    resampling = rsmp("cv", folds = fold)
    
    plan(multisession)
    print('multisession on, start training...')
    rr = resample(task, ensemble, resampling)
    plan(sequential)

    ms = msrs(c("classif.auc","classif.fbeta", "classif.sensitivity"))
    return(data.frame(rf = rr$aggregate(ms)))
}

train_models_cv = function(task, cvfold=10){

    learner_rf = lrn("classif.ranger", predict_type = "prob")

    cv = rsmp("cv", folds = cvfold)

    res_rf = resample(task, learner_rf, cv)

    ms = msrs(c("classif.auc","classif.fbeta", "classif.sensitivity"))

    return(data.frame(rf = res_rf$aggregate(ms)))
}

ML_flow_self = function(raw_matrix, group, method, dataset,fold){
    genes_selected = feature_select(raw_matrix,group)
    mat = data.frame(t(raw_matrix[genes_selected,]), class = group)
    task = TaskClassif$new(id = 'ML', backend = mat, target = "class")
    res = train_models_cv(task,fold)
    return(cbind(res, methods = rep(method,3), datasets = rep(dataset,3)))
}

ML_flow_cross = function(ensemble,matrix1kl,matrix1expr,matrix1ks,matrix2kl,matrix2expr,matrix2ks,group1,group2,id,fold=5){
    #matrix1x: train; matrix2x:test
    matrixkl = cbind(matrix1kl[,-1], matrix2kl[,-1])
    matrixexpr = cbind(matrix1expr[,-1], matrix2expr[,-1])
    matrixks = cbind(matrix1ks[,-1],matrix2ks[,-1])

    genes_selected_A = feature_select(matrix1kl[,-1],group1)
    genes_selected_B = feature_select(matrix1expr[,-1],group1)
    genes_selected_C = feature_select(matrix1ks[,-1],group1)

    mat1 = data.frame(t(matrixkl[genes_selected_A,]))
    mat2 = data.frame(t(matrixexpr[genes_selected_B,]))
    mat3 = data.frame(t(matrixks[genes_selected_C,]))

    print('anova done.')

    colnames(mat1) = paste0("A_", 1:1000)
    colnames(mat2) = paste0("B_", 1:1000)
    colnames(mat3) = paste0("C_", 1:1000)

    comb_mat = cbind(mat1, mat2, mat3)
    comb_mat$target = c(group1, group2)

    task = as_task_classif(comb_mat, target = "target", id = id)
    ensemble$id = id

    ensemble$train(task, row_ids = 1:(dim(matrix1kl)[2]-1))
    pred = ensemble$predict(task, row_ids = dim(matrix1kl)[2]:dim(matrixkl)[2])

    ms = msrs(c("classif.auc","classif.fbeta", "classif.sensitivity"))
    return(data.frame(rf = pred$score(ms)))
}


###TEP2015
kl2015 = read_delim(paste0(MATRIX, 'TEP2015/kl_all.csv'))
ks2015 = read_delim(paste0(MATRIX, 'TEP2015/ks_all.csv'))

md_2015 = read.csv(paste0(METADATA, '/TEP2015_sample_classes.txt'),header=T,sep='\t')
group_2015 = md_2015
group_2015[which(group_2015!='HC'),1] = 'C'
group_2015 = factor(group_2015[match(colnames(kl2015)[-1], rownames(group_2015)),1])

expr_2015 = read_delim(paste0(EXPR_MATRIX,"TEP2015.txt")) %>% z_score_normalize()

res_2015 = ML_flow_self(ensemble_all,kl2015,expr_2015,ks2015,group_2015, 'TEP2015')

###TEP2017
kl2017 = read_delim(paste0(MATRIX, 'TEP2017/kl_all.csv'))
ks2017 = read_delim(paste0(MATRIX, 'TEP2017/ks_all.csv'))

md_2017 = read.csv(paste0(METADATA, '/TEP2017_sample_classes.txt'),header=T,sep='\t')
group_2017 = md_2017
group_2017[which(group_2017!='Healthy Control'),1] = 'C'
group_2017[which(group_2017=='Healthy Control'),1] = 'HC'
group_2017 = factor(group_2017[match(colnames(kl2017)[-1], rownames(group_2017)),1])


expr_2017 = read_delim(paste0(EXPR_MATRIX,"TEP2017.txt")) %>% z_score_normalize()

res_2017 = ML_flow_self(ensemble_all,kl2017, expr_2017,ks2017,group_2017, 'TEP2017')

###TEP2022
kl2022 = read_delim(paste0(MATRIX, 'TEP2022/kl_all.csv'))
ks2022 = read_delim(paste0(MATRIX, 'TEP2022/ks_all.csv'))
ks2022 = ks2022[, -match('SRR15781541',colnames(ks2022))]

md_2022 = read.csv(paste0(METADATA, '/sample_cross_new_TEP_2022.txt'),header=T,sep='\t')
group_2022 = md_2022[match(colnames(kl2022)[-1], md_2022$sample),]$group
group_2022[which(group_2022!='Cancer')] = 'HC'
group_2022[which(group_2022=='Cancer')] = 'C'
group_2022 = factor(group_2022)

expr_2022 = read_delim(paste0(EXPR_MATRIX,"TEP2022.txt"))
expr_2022 = expr_2022[,-match(c('SRR15781676', 'SRR15781678', 'SRR15781815', 'SRR15782862','SRR15781541'), colnames(expr_2022))]
md_2022 = read.csv(paste0(METADATA, '/sample_cross_new_TEP_2022.txt'),header=T,sep='\t')
expr_2022 = expr_2022[,-which(is.na(match(colnames(expr_2022), md_2022$sample)))[2]]
expr_2022 = z_score_normalize(expr_2022)

res_2022 = ML_flow_self(ensemble_all,kl2022, expr_2022,ks2022, group_2022, 'TEP2022')


res_stacking_single = rbind(res_2015, res_2017, res_2022)
res_stacking_single = cbind(msrs=rep(c('AUC','F-beta','Sensitivity'), times=3), 
    res_stacking_single, 
    methods = rep('Comb', 9),
    datasets = rep(c('2015','2017','2022'),each=3))
write.csv(res_stacking_single, paste0(MATRIX, 'stacking_single_all.csv'),row.names=F)


res_2017_2015 = ML_flow_cross(ensemble_all,kl2017, expr_2017,ks2017, kl2015, expr_2015,ks2015, group_2017, group_2015, '2017test2015')
res_2017_2022 = ML_flow_cross(ensemble_all,kl2017, expr_2017, ks2017,kl2022, expr_2022, ks2022,group_2017, group_2022, '2017test2022')
res_2022_2015 = ML_flow_cross(ensemble_all,kl2022, expr_2022,ks2022, kl2015, expr_2015,ks2015, group_2022, group_2015, '2022test2015')
res_2022_2017 = ML_flow_cross(ensemble_all,kl2022, expr_2022,ks2022, kl2017, expr_2017, ks2017,group_2022, group_2017, '2022test2017')

res_stacking_cross = rbind(res_2017_2015, res_2017_2022, res_2022_2015, res_2022_2017)
res_stacking_cross = cbind(msrs=rep(c('AUC','F-beta','Sensitivity'), times=4), 
    res_stacking_cross, 
    methods = rep('Comb', 12),
    datasets = rep(c('2017 test 2015','2017 test 2022','2022 test 2015','2022 test 2017'),each=3))
write.csv(res_stacking_cross, paste0(MATRIX, 'stacking_cross_all.csv'),row.names=F)
