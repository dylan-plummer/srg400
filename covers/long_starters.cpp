// Re-embed the K5 edge labels while preserving the 30 prescribed translation differences.
#include <array>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using U=uint64_t;
using namespace std;
int old[10]={15,41,29,0,48,34,6,19,11,23};
int edge[10], lab[10],delta[50][50],hs[50][50],sh[10];
U wanted=0;long long visits=0,hits=0;int limit=1000;
int dif(int a,int b){return 5*((b/5-a/5+5)%5)+(b%5-a%5+5)%5;}
bool recurse(unsigned assigned,U used,U occupied) {
    if(assigned==1023){
        if(used!=wanted)throw 1;
        int h=0;for(int i=0;i<10;i++)for(int j=0;j<i;j++)h+=hs[lab[i]][lab[j]];
        cout<<h;for(int x:lab)cout<<' '<<x;cout<<'\n';
        return ++hits==limit;
    }
    visits++;
    int best=-1,score=-1;
    for(int i=0;i<10;i++)if(!(assigned>>i&1)){
        int w=0;for(int j=0;j<10;j++)if((assigned>>j&1)&&(edge[i]&edge[j]))w++;
        if(w>score){best=i;score=w;}
    }
    for(int x=25*sh[best];x<25*(sh[best]+1);x++){
        if((occupied>>x)&1)continue;
        U add=0;bool ok=true;
        for(int j=0;j<10;j++)if((assigned>>j&1)&&(edge[best]&edge[j])){
            int t=delta[x][lab[j]];if(t<0){ok=false;break;}
            U bit=U(1)<<t;
            if(!(wanted&bit) || ((used|add)&bit)){ok=false;break;}add|=bit;
        }
        if(!ok)continue;lab[best]=x;
        if(recurse(assigned|(1<<best),used|add,occupied|(U(1)<<x)))return true;
    }
    return false;
}
int main(int argc,char**argv){
    if(argc>1)limit=atoi(argv[1]);
    int k=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)edge[k++]=1<<i|1<<j;
    int ids[3][25];for(auto& l:ids)fill(begin(l),end(l),-1);int id=0;
    for(int i=0;i<50;i++)for(int j=i+1;j<50;j++){
        int t=dif(i%25,j%25),sheet=i/25+j/25;
        if(sheet!=1)t=min(t,dif(j%25,i%25));
        int dx=t/5,dy=t%5;bool adjacent=sheet==1?(dy+3*dx*dx)%5==0:!dx && (dy==(sheet?2:1)||dy==(sheet?3:4));
        hs[i][j]=hs[j][i]=adjacent;
        if(!adjacent&&ids[sheet][t]<0)ids[sheet][t]=id++;
        delta[j][i]=delta[i][j]=ids[sheet][t];
    }
    if(id!=42)return 1;
    for(int j=0;j<10;j++){
        sh[j]=old[j]/25;
        for(int i=0;i<j;i++)if(edge[i]&edge[j])wanted|=U(1)<<delta[old[i]][old[j]];
    }
    if(popcount(wanted)!=30)return 2;
    lab[0]=0;recurse(1,0,1);
    cerr<<"visited "<<visits<<" partial embeddings, emitted "<<hits<<" long starters\n";
}
