# Loading Pima-Surui Data
data(geno.ps.x)
data(geno.ps.y)
Geno <- geno.ps.x

# Infer population structure with first 30 markers
a = PSMix( K=2,Geno[,1:30],eps=1e-8,seed=as.integer(1+abs(rcauchy(1,scale=180))), MarkerVar=FALSE )
t( round(a[[1]],3) )

# Infer population structure with all 100 markers
a = PSMix( K=2,Geno,eps=1e-8,seed=as.integer(1+abs(rcauchy(1,scale=180))), MarkerVar=FALSE )
t( round(a[[1]],3) )
