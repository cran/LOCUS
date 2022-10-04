Ltrinv<-function(x,V,d=TRUE){ Y = matrix(0,ncol = V,nrow = V);
Y[upper.tri(Y,d)]=x;return(Y + t(Y) -d*diag(diag(Y)))  }
