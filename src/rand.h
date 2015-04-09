double sTruncNorm(double bd, double mu, double var, int lower);
double TruncNorm(double lb, double ub, double mu, double var, int invcdf);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);


