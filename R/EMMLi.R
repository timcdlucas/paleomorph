#################################################################
##  EMMLi 
##
##   Inputs: corr = lower triangle coeffecient matrix, 
##  					mod = landmark classification, 
##				varlist = list of classfications (mod columns), 
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

EMMLi = function(corr, N_sample, mod, varlist,saveAs){
  
  varlist[length(varlist)+1] = "mod$No.modules"
  mod$No.modules = 1
  
  corr[upper.tri(corr, diag=T)] = NA
  corr_list = (as.array(corr[!is.na(corr)])) # array of coefficient matrix, NAs are removed.
  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)] #a symmetric matrix formed from the coefficient matrix, only used to find intermodular coefficients.
  
  all_modules = list()
  for(colnam in varlist){
    col = array(eval(parse(text = colnam)))
    
    modNF = na.omit(cbind(mod[,1],col)) # na.omit() used to remove unclassified landmarks.
    w = unique(modNF[,2])
    
    modules = list()
    for(i in seq(length(w))){
      fg = modNF[modNF[,2] == w[i],] # identify landmarks within a class
      
      l = corr[fg[,1],fg[,1]] # coefficients between identified landmarks.
      modules[[i]] = (as.array(l[!is.na(l)]))
    }
    names(modules) = paste("Module",w)
    
    between_mod = list()
    betweenModules = list()
    withinModules = list()
    unintegrated = list()
    betweenFloat = list()
    
    if (length(w) > 1){ #check that the num. of modules is greater than 1
      cb = combn(w,2) # all possible combination of modules
      for (i in seq(dim(cb)[2])){
        fg1 = modNF[modNF[,2] == cb[1,i],]
        fg2 = modNF[modNF[,2] == cb[2,i],]
        
        between = symmet[as.integer(setdiff(fg2[,1], fg1[,1])), as.integer(setdiff(fg1[,1], fg2[,1]))] #setdiff(A,B) - present in A but not in B
        between_mod[[i]] = between[!is.na(between)]
      }
      names(between_mod) = paste(cb[1,],"to",cb[2,])
      
      betweenModules['betweenModules'] = list(as.vector(rle(unlist(between_mod))$values))
      withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
    }
    
    unintegrated_list = setdiff(corr_list,unlist(c(modules,between_mod)))
    unintegrated['unintegrated'] = list(as.vector(rle(unlist(unintegrated_list))$values))
    # unintegrated['unintegrated+between'] = list(as.vector(rle(unlist(c(unintegrated_list,betweenModules)))$values))
    
    if(length(unintegrated_list) != 0){
      betweenFloat['betweenFloat'] = list(as.vector(rle(unlist(c(betweenModules['betweenModules'],unintegrated['unintegrated'])))$values))
      all_modules[[colnam]] = c(modules,between_mod,withinModules,betweenModules,betweenFloat,unintegrated)
    }else{
      all_modules[[colnam]] = c(modules,between_mod,withinModules,betweenModules,unintegrated)
      
    }
    
  }
  
  LogL = function(z_r,z_p) {-0.5*log(var) - ((z_r - z_p)^2) / (2*var)}
  
  output = list()
  maxlogL = list()
  logp = list()
  for (m in seq(length(all_modules))){
    maxres = matrix(,nrow=2,ncol=(length(all_modules[[m]]))) # table of max. likelihood and rho
    dimnames(maxres) = list(c("MaxL", "MaxL_p"),c(names(all_modules[[m]][seq((length(all_modules[[m]])))])))
    
    for(g in seq((length(all_modules[[m]])))){
      r = unlist(unname(all_modules[[m]][g]))
      n_value = length(r)
      z_r = 0.5*log((1+abs(r))/(1-abs(r)))
      n = N_sample #24 # or 1830
      var = 1/(n-3)
      
      p = seq(0, 0.99, 0.01) #rho
      z_p = 0.5*log((1+p)/(1-p)) #zeta
      
      LogL_table = outer(z_r,z_p,LogL)
      Likelihoods = colSums(LogL_table)
      MaxIndex = which.max(Likelihoods)
      
      MaxL = Likelihoods[MaxIndex]
      MaxL_p = p[MaxIndex]
      maxres[1,g] = MaxL
      maxres[2,g] = MaxL_p
    }
    
    output[[names(all_modules)[m]]] = maxres # list of maximum likelihood for all modules.
    
    if (dim(output[[m]])[2] == 2){ #calculate sum of modular likelihood and k 
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated',drop=F]}
      
      maxlogL[[names(all_modules)[m]]][['default']] = c(output[[m]][1], 2)
      logp[[names(all_modules)[m]]][['default']] = output[[m]][,1,drop=F]
      
    } else if(dim(output[[m]])[2] == 6){   
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_between = output[[m]]['MaxL',grep('Module |betweenModules|unintegrated',names(output[[m]][1,]))]
      mod_between_p = output[[m]][,grep('Module |betweenModules|unintegrated',names(output[[m]][1,]))]
      if(mod_between['unintegrated'] == 0){
        K = length(mod_between)-1
      }else{
        K = length(mod_between)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K+1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      
      within_between = output[[m]]['MaxL',c('withinModules','betweenModules','unintegrated')]
      within_between_p = output[[m]][,c('withinModules','betweenModules','unintegrated')]
      if(within_between['unintegrated'] == 0){
        K = length(within_between)-1
      }else{
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K+1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
    } else{      
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_to = output[[m]]['MaxL',grep('Module |to |unintegrated',names(output[[m]][1,]))]
      mod_to_p = output[[m]][,grep('Module |to |unintegrated',names(output[[m]][1,]))]
      if(mod_to['unintegrated'] == 0){
        K = length(mod_to)-1
      }else{
        K = length(mod_to)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + sep.between']] = c(sum(mod_to), K+1)
      logp[[names(all_modules)[m]]][['sep.Mod + sep.between']] = mod_to_p
      
      within_between = output[[m]]['MaxL',c('withinModules','betweenModules','unintegrated')]
      within_between_p = output[[m]][,c('withinModules','betweenModules','unintegrated')]
      if(within_between['unintegrated'] == 0){
        K = length(within_between)-1
      }else{
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K+1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
      
      mod_between = output[[m]]['MaxL',grep('Module |betweenModules|unintegrated',names(output[[m]][1,]))]
      mod_between_p = output[[m]][,grep('Module |betweenModules|unintegrated',names(output[[m]][1,]))]
      if(mod_between['unintegrated'] == 0){
        K = length(mod_between)-1 #length(output[[m]]['MaxL',grep('Module ',names(output[[m]][1,]))])*2+1
      }else{
        K = length(mod_between) #length(output[[m]]['MaxL',grep('Module ',names(output[[m]][1,]))])*2+2
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K+1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      to_within = output[[m]]['MaxL',grep('to |withinModules|unintegrated',names(output[[m]][1,]))]
      to_within_p = output[[m]][,grep('to |withinModules|unintegrated',names(output[[m]][1,]))]
      if(to_within['unintegrated'] == 0){
        K = length(to_within)-1 #length(output[[m]]['MaxL',grep('Module |to',names(output[[m]][1,]))])+1
      }else{
        K = length(to_within) #length(output[[m]]['MaxL',grep('Module |to',names(output[[m]][1,]))])+2
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + sep.between']] = c(sum(to_within), K+1)
      logp[[names(all_modules)[m]]][['same.Mod + sep.between']] = to_within_p
      
      if(output[[m]]['MaxL',grep('unintegrated',names(output[[m]][1,]))] != 0){
        
        sepmod_samebetweenunintegrated = output[[m]][,grep('Module |betweenFloat',names(output[[m]][1,]))]
        K = length(sepmod_samebetweenunintegrated['MaxL',])
        maxlogL[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = c(sum(sepmod_samebetweenunintegrated['MaxL',]), K+1)
        logp[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = sepmod_samebetweenunintegrated
        
        samemod_samebetweenunintegrated = output[[m]][,grep('withinModules|betweenFloat',names(output[[m]][1,]))]
        K = length(samemod_samebetweenunintegrated['MaxL',])
        maxlogL[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = c(sum(samemod_samebetweenunintegrated['MaxL',]), K+1)
        logp[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = samemod_samebetweenunintegrated
        
      }
    }
  }
  
  results = matrix(unlist(maxlogL),ncol=2,byrow=TRUE)
  a = names(unlist(maxlogL))[seq(1,dim(results)[1]*dim(results)[2],dim(results)[2])]
  dimnames(results) = list(a,c('MaxL', 'K'))
  
  n = length(which(!is.na(corr) == TRUE))
  AICc = apply(results, 1, function(x) -2*x['MaxL']+2*x['K'] + (2*x['K']*(x['K']+1))/(n-x['K']-1))
  dAICc = AICc-min(AICc)
  Model_L = exp(-0.5*dAICc)
  Post_Pob = Model_L/ sum(Model_L)
  
  results = cbind(results,n,AICc,dAICc,Model_L,Post_Pob)
  
  n = names(all_modules)
  nm = unlist(strsplit(n, split='mod\\$'))[seq(2,2*length(n),2)]
  
  a = names(unlist(maxlogL))[seq(1,dim(results)[1]*2,2)]
  b = unlist(strsplit(a, split='mod\\$'))
  b = unlist(strsplit(b, split='1'))
  rownames(results) = b
  
  #t = matrix(,nrow=8,ncol=2)
  
  s=1
  i = length(nm)
  o = order(results[grep(nm[i],rownames(results)),2])
  t = results[grep(nm[i],rownames(results)),,drop=F][o,,drop=F]
  s= s+length(o)
  
  for(i in 1:(length(nm)-1)){
    o = order(results[grep(paste(nm[i],'.s',sep=""),rownames(results)),2])
    t = rbind(t,results[grep(paste(nm[i],'.s',sep=""),rownames(results)),,drop=F][o,,drop=F])
    #rownames(t[s:(s+length(o)-1),]) <- rownames(results[grep(nm[i],rownames(results)),,drop=F][o,,drop=F])
    s= s+length(o)
  }
  
  results = t
  
  rholist = list()
  h=1
  #rholist[h] = logp[1]
  #h=1
  for (i in 1:(length(logp))){
    #print(i)
    for (j in 1:length(logp[[i]])){
      
      rholist[h] = logp[[i]][j]
      
      h = h+1
    }
    
  } 
  
  rho_output = rholist[which(Post_Pob > 0.01)]
  
 
  rholist_name = names(which(Post_Pob > 0.01))
  rholist_name = unlist(strsplit(rholist_name, split='mod\\$'))
  rholist_name = unlist(strsplit(rholist_name, split='1'))
  
  write.table(results, file=saveAs, row.names=TRUE, col.names=NA, sep=",")
  cat("\n\n", file = saveAs, append = TRUE)
  for(q in 1:length(rho_output)){
    cat(rholist_name[q],"\n", file = saveAs, append = TRUE)
    write.table(rho_output[q], saveAs, row.names=T,col.names=NA, sep=",", append=TRUE)
    cat("\n", file = saveAs, append = TRUE)
  }
}
