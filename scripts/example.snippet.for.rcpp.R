


var_d = function(d,coor1,coor2){
  i = match(coor1,PCA_DT[,grid_id])
  j = match(coor2,PCA_DT[,grid_id])
  charvec = paste0("PC",1:d)
  sq_diff_sum=0
  for(sum_ind in 1:d){
    sq_diff_sum = sq_diff_sum +(PCA_DT[,eval(parse(text = charvec[sum_ind]))][i] - 
                                  PCA_DT[,eval(parse(text = charvec[sum_ind]))][j]  )^2
  }
  return(sq_diff_sum)
}

# var_sc_prereq is a data table with some million rows, and a column for each number in dvec:

dummy_fct = function(d){
  return( var_d(d,var_sc_prereq[,grid_id1],var_sc_prereq[,grid_id2]))
}

dummy_list = mclapply(dvec,dummy_fct,mc.cores = 10)