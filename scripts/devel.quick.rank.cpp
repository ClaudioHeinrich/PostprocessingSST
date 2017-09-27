NumericVector QuickRank(NumericMatrix M,int i_Obs)
{
  int i;
  std::vector<pair<double,int>> a;
  for(i = 0; i < i_obs; ++i)
    {
      a[i] = make_pair(M(0,i),1);
    }
  for(i = i_obs; i < M.ncol(); ++i)
    {
      a[i] = make_pair(M(0,i),0);
    }

  sort(a.begin(), a.end(),ordering());

  NumericVector b = 
  for(i = 0; i < i_obs;++i)b[i] = 1.0;
  NumericVector
			     
}

struct ordering {
    bool operator ()(pair<double, int> const& a, 
                     pair<double, int> const& b) {
        return a.first < b.first;
    }
};
