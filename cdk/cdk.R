pkgs <- installed.packages()
if (! 'rcdk' %in% pkgs[,1]){
  install.packages('rcdk')
}
cdk_fingerprint <- function(smi, tps = c('pubchem', 'kr')){
  library(rcdk)
  mol <- parse.smiles(smi)[[1]]
  fps <- lapply(tps, function(tp) get.fingerprint(mol, type=tp))
  vec <- lapply(fps, function(fp) {
    v <- rep(0, fp@nbit)
    v[fp@bits] <- 1
    v
  })
  vec <- unlist(vec)
  return(vec)
}