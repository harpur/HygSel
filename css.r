#source('css.r')                                            # path for Folder

css=function(input.data, tests=c('fst','daf','xpehh'), ihs=F)
{
test.stat=input.data
gt=tests; gt1=length(gt)  ## default
if(ihs==T){test.stat$ihs=abs(test.stat$ihs)}
gt1=length(gt)
test.rank = test.stat
for (i in 1:length(gt)){
test.rank[,gt[i]]    <- rank(test.rank[,gt[i]],       na.last = "keep")
}
for (i in 1:length(gt)){
test.rank[,gt[i]]    <- test.rank[,gt[i]]    / (max(test.rank[,gt[i]],    na.rm = T) + 1)
}

test.z = test.stat
for (i in 1:length(gt)){
test.z[,gt[i]]    <- qnorm(test.rank[,gt[i]])
}
test.mz <- apply(test.z[,gt], 1, mean)
p.test.mz <- pnorm(test.mz, 0, 1/sqrt(gt1), lower.tail = F) #  zbar ~ N(0, 1/sqrt(3))
### add results to the input file
p.test.mz=as.data.frame(p.test.mz); test.mz=as.data.frame(test.mz);
results= cbind(test.mz, p.test.mz, -log10(p.test.mz)); names(results)=c('mz','p_mz','css');
ouput.data=cbind(input.data, results)
ouput.data=as.data.frame(ouput.data)
print("MeanZ, P of meanZ, -log10 of P of meanZ has been calculated and added to the output file as: 'mz', 'p_mz', 'css'")
ouput.data

  }  # End... CSS




##################################################################################