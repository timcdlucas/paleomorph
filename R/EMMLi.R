#################################################################
##  EMMLi 
##
##   Inputs: corr = lower triangle correlation matrix, 
##       N_sample = the number of samples used to calculate the correlation matrix
##  					mod = landmark classification, 
##				 saveAs = specify where to save the results.
##			  
##  Output: A .csv file, resembling the AIC worksheet. 
##
##  Prabu (p.siva@ucl.ac.uk)
#################################################################

##############
# The MIT License (MIT)

# Copyright (c) 2015 Pratheesan Sivasubramaniam

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#############


#' Evaluating modularity with maximum likelihood
#'
#' Calculates the AICc values of different models of modularity.
#'
#'@param corr Lower triangle or full correlation matrix. n x n square matrix for n landmarks.
#'@param N_sample The number of samples used
#'@param mod A data frame defining the models. The first column should contain the landmark names. Subsequent columns should define which landmarks are contained with which module. If a landmark should be ignored for a specific model, the element should be NA.
#'@param saveAs A character string defining the filename and path for where to save output.
#'
#'@export
#'@return NULL. The output is saved to the file defined by the saveAs argument.


EMMLi = function(corr, N_sample, mod, saveAs){
  
  # Check inputs
  stopifnot(is.numeric(corr), dim(corr)[1] == dim(corr)[2], is.numeric(N_sample), N_sample > 0, is.character(saveAs))
  # Check that models are given as integers
  if (!all(sapply(mod[, -1], function(i) i %% 1 == 0) | is.na(mod[, -1]))) stop('mod should contain a column of names and then columns of integers defining models')
  if (dim(corr)[1] != dim(corr)[2]) stop('corr should be a square matrix.')
  
  # Create null model
  mod$No.modules = 1

  # Create varlist variable
  varlist = paste0('mod$', names(mod)[-1])

  # make the upper triangle of corr NA, we only use the lower triangle.
  corr[upper.tri(corr, diag = T)] = NA

  # array of coefficient matrix, NAs are removed.
  corr_list = (as.array(corr[!is.na(corr)])) 
  symmet = corr

  # a symmetric matrix formed from the coefficient matrix, only used to find intermodular coefficients.
  
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)] 
  all_modules = list()
  for(colnam in varlist){
    col = array(eval(parse(text = colnam)))
    
    # na.omit() used to remove unclassified landmarks.
    modNF = stats::na.omit(cbind(mod[, 1], col)) 
    w = unique(modNF[, 2])
    
    modules = list()
    for(i in seq(length(w))){
      # identify landmarks within a class
      fg = modNF[modNF[, 2] == w[i], ] 
      
      # coefficients between identified landmarks.
      l = corr[fg[, 1], fg[, 1]] 
      modules[[i]] = (as.array(l[!is.na(l)]))
    }
    names(modules) = paste("Module", w)
    
    between_mod = list()
    betweenModules = list()
    withinModules = list()
    unintegrated = list()
    betweenFloat = list()
    
    if (length(w) > 1){ #check that the num. of modules is greater than 1

      # all possible combination of modules
      cb = utils::combn(w, 2) 
      for (i in seq(dim(cb)[2])){
        fg1 = modNF[modNF[, 2] == cb[1, i], ]
        fg2 = modNF[modNF[, 2] == cb[2, i], ]
        
        # setdiff(A,B) - present in A but not in B
        between = symmet[as.integer(setdiff(fg2[, 1], fg1[, 1])), as.integer(setdiff(fg1[, 1], fg2[, 1]))] 
        between_mod[[i]] = between[!is.na(between)]
      }
      names(between_mod) = paste(cb[1, ], "to", cb[2, ])
      
      betweenModules['betweenModules'] = list(as.vector(rle(unlist(between_mod))$values))
      withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
    }
    
    unintegrated_list = setdiff(corr_list,unlist(c(modules,between_mod)))
    unintegrated['unintegrated'] = list(as.vector(rle(unlist(unintegrated_list))$values))
    # unintegrated['unintegrated+between'] = list(as.vector(rle(unlist(c(unintegrated_list,betweenModules)))$values))
    
    if (length(unintegrated_list) != 0){
      betweenFloat['betweenFloat'] = list(as.vector(rle(unlist(c(betweenModules['betweenModules'], unintegrated['unintegrated'])))$values))
      all_modules[[colnam]] = c(modules, between_mod, withinModules, betweenModules, betweenFloat, unintegrated)
    } else {
      all_modules[[colnam]] = c(modules, between_mod, withinModules, betweenModules, unintegrated)
    }
    
  }
  
  LogL = function(z_r, z_p) {-0.5 * log(var) - ((z_r - z_p)^2) / (2 * var)}
  
  output = list()
  maxlogL = list()
  logp = list()
  for (m in seq(length(all_modules))){
    # table of max. likelihood and rho
    maxres = matrix(, nrow = 2, ncol = (length(all_modules[[m]]))) 
    dimnames(maxres) = list(c("MaxL", "MaxL_p"), c(names(all_modules[[m]][seq((length(all_modules[[m]])))])))
    
    for(g in seq((length(all_modules[[m]])))){
      r = unlist(unname(all_modules[[m]][g]))
      n_value = length(r)
      z_r = 0.5 * log((1 + abs(r)) / (1 - abs(r)))
      n = N_sample 
      var = 1 / (n - 3)
      
      #rho
      p = seq(0, 0.99, 0.01) 

      #zeta
      z_p = 0.5 * log((1 + p)/(1 - p)) 
      
      LogL_table = outer(z_r, z_p, LogL)
      Likelihoods = colSums(LogL_table)
      MaxIndex = which.max(Likelihoods)
      
      MaxL = Likelihoods[MaxIndex]
      MaxL_p = p[MaxIndex]
      maxres[1, g] = MaxL
      maxres[2, g] = MaxL_p
    }
    
    # list of maximum likelihood for all modules.
    output[[names(all_modules)[m]]] = maxres 
    
    #calculate sum of modular likelihood and k 
    if (dim(output[[m]])[2] == 2){ 
      
      maxlogL[[names(all_modules)[m]]][['default']] = c(output[[m]][1], 2)
      logp[[names(all_modules)[m]]][['default']] = output[[m]][, 1, drop = FALSE]
      
    } else if (dim(output[[m]])[2] == 6){   
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_between = output[[m]]['MaxL', grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      mod_between_p = output[[m]][,grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      if (mod_between['unintegrated'] == 0){
        K = length(mod_between) - 1
      } else {
        K = length(mod_between)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      
      within_between = output[[m]]['MaxL', c('withinModules','betweenModules','unintegrated')]
      within_between_p = output[[m]][, c('withinModules','betweenModules','unintegrated')]
      if (within_between['unintegrated'] == 0){
        K = length(within_between) - 1
      } else {
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K+1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
    } else {      
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_to = output[[m]]['MaxL', grep('Module |to |unintegrated', names(output[[m]][1, ]))]
      mod_to_p = output[[m]][, grep('Module |to |unintegrated', names(output[[m]][1,]))]
      if (mod_to['unintegrated'] == 0){
        K = length(mod_to)-1
      } else {
        K = length(mod_to)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + sep.between']] = c(sum(mod_to), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + sep.between']] = mod_to_p
      
      within_between = output[[m]]['MaxL', c('withinModules', 'betweenModules', 'unintegrated')]
      within_between_p = output[[m]][, c('withinModules', 'betweenModules', 'unintegrated')]
      if (within_between['unintegrated'] == 0){
        K = length(within_between) - 1
      } else {
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K + 1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
      
      mod_between = output[[m]]['MaxL', grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      mod_between_p = output[[m]][, grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      if (mod_between['unintegrated'] == 0){
        K = length(mod_between) - 1 
      } else {
        K = length(mod_between) 
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      to_within = output[[m]]['MaxL', grep('to |withinModules|unintegrated', names(output[[m]][1, ]))]
      to_within_p = output[[m]][, grep('to |withinModules|unintegrated', names(output[[m]][1, ]))]
      if (to_within['unintegrated'] == 0){
        K = length(to_within) - 1
      } else {
        K = length(to_within) 
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + sep.between']] = c(sum(to_within), K + 1)
      logp[[names(all_modules)[m]]][['same.Mod + sep.between']] = to_within_p
      
      if (output[[m]]['MaxL', grep('unintegrated', names(output[[m]][1, ]))] != 0){
        
        sepmod_samebetweenunintegrated = output[[m]][, grep('Module |betweenFloat',names(output[[m]][1, ]))]
        K = length(sepmod_samebetweenunintegrated['MaxL', ])
        maxlogL[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = c(sum(sepmod_samebetweenunintegrated['MaxL', ]), K + 1)
        logp[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = sepmod_samebetweenunintegrated
        
        samemod_samebetweenunintegrated = output[[m]][, grep('withinModules|betweenFloat', names(output[[m]][1, ]))]
        K = length(samemod_samebetweenunintegrated['MaxL', ])
        maxlogL[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = c(sum(samemod_samebetweenunintegrated['MaxL', ]), K + 1)
        logp[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = samemod_samebetweenunintegrated
        
      }
    }
  }
  
  results = matrix(unlist(maxlogL), ncol = 2, byrow = TRUE)
  a = names(unlist(maxlogL))[seq(1, dim(results)[1] * dim(results)[2], dim(results)[2])]
  dimnames(results) = list(a, c('MaxL', 'K'))
  
  n = length(which(!is.na(corr) == TRUE))
  AICc = apply(results, 1, function(x) -2 * x['MaxL'] + 2 * x['K'] + (2 * x['K'] * (x['K'] + 1)) / (n - x['K'] - 1))
  dAICc = AICc - min(AICc)
  Model_L = exp(-0.5 * dAICc)
  Post_Pob = Model_L / sum(Model_L)
  
  results = cbind(results, n, AICc, dAICc, Model_L, Post_Pob)
  
  n = names(all_modules)
  nm = unlist(strsplit(n, split = 'mod\\$'))[seq(2, 2 * length(n), 2)]
  
  a = names(unlist(maxlogL))[seq(1, dim(results)[1] * 2, 2)]
  b = unlist(strsplit(a, split = 'mod\\$'))
  b = unlist(strsplit(b, split = '1$'))
  rownames(results) = b

  
  s = 1
  i = length(nm)
  o = order(results[grep(nm[i], rownames(results)), 2])
  t = results[grep(nm[i], rownames(results)), , drop = FALSE][o, , drop = FALSE]
  s = s + length(o)
  
  for(i in 1:(length(nm) - 1)){
    o = order(results[grep(paste(nm[i], '.s', sep = ""), rownames(results)), 2])
    t = rbind(t, results[grep(paste(nm[i], '.s', sep = ""), rownames(results)), , drop = FALSE][o, , drop = FALSE])
    s = s + length(o)
  }
  
  results = t
  
  rholist = list()
  h = 1

  for (i in 1:(length(logp))){
    #print(i)
    for (j in 1:length(logp[[i]])){
      
      rholist[h] = logp[[i]][j]
      
      h = h + 1
    }
    
  } 
  
  rho_output = rholist[which(Post_Pob > 0.01)]
  
 
  rholist_name = names(which(Post_Pob > 0.01))
  rholist_name = unlist(strsplit(rholist_name, split = 'mod\\$'))
  rholist_name = unlist(strsplit(rholist_name, split = '1$'))
  
  utils::write.table(results, file = saveAs, row.names = TRUE, col.names = NA, sep = ",")
  cat("\n\n", file = saveAs, append = TRUE)
  for(q in 1:length(rho_output)){
    cat(rholist_name[q], "\n", file = saveAs, append = TRUE)
    write(paste(c('', colnames(rho_output[[q]])), collapse = ','), saveAs, append = TRUE)
    write.table(rho_output[q], saveAs, row.names = TRUE, col.names = FALSE, sep = ",", append = TRUE)
    cat("\n", file = saveAs, append = TRUE)
  }
}
