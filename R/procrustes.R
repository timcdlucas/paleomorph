setwd("A:/CBER/CBER Programmer/ITEM/cvm")

a = unlist(read.csv("input.csv", sep=",", header=F))
ar <- array(a, dim=c(3,55, 18))
#ar2 = array(ar, dim=c(18,55, 3))

ar <- aperm(ar, c(3,2,1))

ar[,1,1]

scorea[a_, m_, n_] := Module[ {i, j, k, v, score},
                              score = 0.0;
                              For[i = 1, i <= m, i++,
                                  For[j = i+1, j <= m, j++,
                                      For[k = 1, k <= n, k++,
                                          If[StringQ[a[[i,k]]], Continue[]];
                                          If[StringQ[a[[j,k]]], Continue[]];
                                          v = a[[i,k]] - a[[j,k]];
                                          score += v.v;
                                          ];
                                      ];
                                  ];
                              Return[score];
                              ];
mat = ar
m = 18
n = 55

scorea = function(mat,m,n){
  sumw = 0
  for(i in seq(m)){
    for(j in seq(i,m)){
      for(k in seq(n)){
        nzi = length(mat[i,k,][mat[i,k,] == 0] == TRUE)
        nzj = length(mat[j,k,][mat[j,k,] == 0] == TRUE)
        nzsum = nzi + nzj
        if(nzsum == 0){
          w = mat[i,k,] - mat[j,k,]
          sumw = sumw + w %*% w
        }
      }
    }
  }
  return(sumw)
}

scorea(ar,m,n)
#ar

#dist[v1_,v2_] := Sqrt[ (v1-v2).(v1-v2) ];
dist = function(v1,v2){
  return(sqrt((v1-v2)%*%(v1-v2)))
}

#fnction deltaa[a1_, a2_, m_, n_]
deltaa = function(a,na,m,n){
  delta = 0
  for(i in seq(m)){
    for(j in seq(n)){
      if (sum(a[i,j,]) != 0 && sum(na[i,j,]) != 0){
        delta = delta + dist(a[i,j,],na[i,j,])
      }
    }
  }
  return(delta)
}

deltaa(a,na,m,n)

rpdecompose[M_] 
rpdecompose = function(m){
  
  
}

pcrstep[a_,m_,n_] 
pcrstep = function(na,m,n){
  for(count in seq(1000)){
    na2 = na
    for(i in seq(2,m)){
      ta = matrix(0,nrow=n,ncol=3)
      for(j in seq(m)){
        if(i != j){
          ta = ta + na[j,,] 
        }
        c = t(na[i,,]) %*% ta
        r = rpdecompose
      }
    }
    
  }
  
}


pctstep[a_,m_,n_] := Module[ {i, j, k, c, bx, by, bz, v, na},
                             
(** The linear algebra here is a little tricky. * We have M equations in M unknowns, but one equation is redundant* and the solution space is invariant under a common translation.* For now, we deal with this by throwing away the first equation* and forcing the first unknown to equal 0.* A more robust solution might be: retain all M equations and solve by* least squares, and force the sum of the unknowns to 0.*)
                             
                             c  = Table[0.0, {i,1,m-1}, {j,1,m-1}];
                             bx = Table[0.0, {i,1,m-1}];
                             by = Table[0.0, {i,1,m-1}];
                             bz = Table[0.0, {i,1,m-1}];
                             
                             For[i = 2, i <= m, i++,
                                 For[j = 1, j <= m, j++,
                                     If[i == j, Continue[]];
                                     For[k = 1, k <= n, k++,
                                         If[StringQ[a[[i,k]]] || StringQ[a[[j,k]]], Continue[]];
                                         c[[i-1,i-1]] += 1.0;
                                         If[j > 1, c[[i-1,j-1]] -= 1.0];
                                         bx[[i-1]] += a[[i,k,1]] - a[[j,k,1]];
                                         by[[i-1]] += a[[i,k,2]] - a[[j,k,2]];
                                         bz[[i-1]] += a[[i,k,3]] - a[[j,k,3]];
                                         ];
                                     ];
                                 ];
                             
                             (* FIXME: emit warning here if matrix is badly conditioned *)
                             
                             (*
                                * FIXME: this could be made more efficient by using
                              * LUDecomposition[] and LUBackSubstitution[] instead of LinearSolve[]
                              *)
                             v = Transpose[ {LinearSolve[c, bx],
                                             LinearSolve[c, by],
                                             LinearSolve[c, bz]
                             }];
                             
                             na = Table[
                    If[i==1, a[[1]], lshift[a[[i]], -v[[i-1]]]],
{i,1,m}
];
Print["tstep: score=", scorea[na,m,n], " delta=", deltaa[a,na,m,n]];
Return[na];
];

function pctstep[ar,18,55]

mat = ar
m = 18
n =55

c = matrix(0,nrow=m-1,ncol=m-1)
bx = matrix(0,nrow=1,ncol=m-1)
by = matrix(0, nrow=1,ncol=m-1)
bz = matrix(0, nrow=1,ncol=m-1)

for(i in seq(2,m)){
  for(j in seq(m)){
    if (i != j) {
      for(k in seq(n)){
        if (sum(a[i,k,]) != 0 && sum(a[j,k,]) != 0) {
          c[i-1,i-1] = c[i-1,i-1] + 1
          if (j > 1) {
            c[i-1,j-1] = c[i-1,j-1] - 1
          }
      #}
          bx[i-1] = bx[i-1] +  a[i,k,1] - a[j,k,1]
          by[i-1] = by[i-1] +  a[i,k,2] - a[j,k,2]
          bz[i-1] = bz[i-1] +  a[i,k,3] - a[j,k,3]
      }
          
        }
     }
    }
}

v = cbind(solve(c,t(bx)),solve(c,t(by)),solve(c,t(bz)))

lshift[L_,v_] := Module[ {i},
                         Table[If[StringQ[L[[i]]], L[[i]], L[[i]]+v], {i,1,Length[L]}]
                         ];

na = array(NA, dim=c(55, 3,18))
na[,,1] = a[1,,]

for(i in seq(2,18)){
  na[,,i] = t(t(a[i,,])-v[i-1,])  
}

na[na < 0] = 0
na <- aperm(na, c(3,1,2))
h = scorea(na,m,n)


