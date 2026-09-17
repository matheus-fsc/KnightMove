#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
typedef unsigned long long u64;
int n,V,deg[64],adjv[64][8],LMAX;
static void dfs(int start,int cur,int depth,u64 vis,u64*c){
  int d=deg[cur];
  for(int i=0;i<d;i++){
    int w=adjv[cur][i];
    if(w==start){ if(depth>=3) c[depth]++; }
    else if(w>start && !((vis>>w)&1ULL) && depth<LMAX)
      dfs(start,w,depth+1,vis|(1ULL<<w),c);
  }
}
int main(int argc,char**argv){
  n=atoi(argv[1]); LMAX=atoi(argv[2]); V=n*n;
  int dr[8]={-2,-2,-1,-1,1,1,2,2}, dc[8]={-1,1,-2,2,-2,2,-1,1};
  for(int r=0;r<n;r++)for(int c=0;c<n;c++){int u=r*n+c;
    for(int k=0;k<8;k++){int rr=r+dr[k],cc=c+dc[k];
      if(rr>=0&&rr<n&&cc>=0&&cc<n) adjv[u][deg[u]++]=rr*n+cc;}}
  /* tarefas = (start, primeiro vizinho > start) */
  int ts[4096],tw[4096],T=0;
  for(int s=0;s<V;s++) for(int i=0;i<deg[s];i++) if(adjv[s][i]>s){ ts[T]=s; tw[T]=adjv[s][i]; T++; }
  u64 tot[64]={0};
  double t0=omp_get_wtime();
  #pragma omp parallel
  { u64 loc[64]={0};
    #pragma omp for schedule(dynamic,1)
    for(int t=0;t<T;t++){ int s=ts[t],w=tw[t];
      dfs(s,w,2,(1ULL<<s)|(1ULL<<w),loc); }
    #pragma omp critical
    for(int L=0;L<64;L++) tot[L]+=loc[L];
  }
  double el=omp_get_wtime()-t0;
  u64 g=0;
  for(int L=3;L<=LMAX;L++) if(tot[L]){ printf("L=%2d: %llu\n",L,tot[L]); g+=tot[L]; }
  printf("TOTAL = %llu   [%.1fs, %d threads]\n",g,el,omp_get_max_threads());
  return 0;
}
