// Enumerate long starters using the remaining difference classes after two short orbits.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
using namespace std;
using U=uint64_t;
int cls[50][50],edges[10],label[10],sheet[10],order[10],prior[10][10],nprior[10],pat;
U target;
long long nodes=0,hits=0;
void visit(int depth,U used,U taken){
    if(depth==10){
        if(used!=target)throw 1;
        cout<<target<<' '<<pat;for(auto x:label)cout<<' '<<x;cout<<'\n';hits++;
        return;
    }
    nodes++;int i=order[depth],base=25*sheet[i];
    for(int x=base;x<base+(depth==0?1:25);x++){
        if((taken>>x)&1)continue;
        U occupied=used;
        bool ok=true;
        for(int k=0;k<nprior[depth];k++){
            int c=cls[x][label[prior[depth][k]]];
            if(c<0||!(target>>c&1)||((occupied>>c)&1)){ok=false;break;}
            occupied|=U(1)<<c;
        }
        if(ok){label[i]=x;visit(depth+1,occupied,taken|(U(1)<<x));}
    }
}
int main(){
    int k=0;for(int a=0;a<5;a++)for(int b=a+1;b<5;b++)edges[k++]=1<<a|1<<b;
    for(auto& row:cls)for(int& c:row)cin>>c;
    int count=0;
    while(cin>>target>>pat){
        for(int i=0;i<10;i++)sheet[i]=(pat>>i)&1;
        int availability[3]={};
        for(int a=0;a<50;a++)for(int b=0;b<a;b++)if(cls[a][b]>=0&&(target>>cls[a][b]&1))availability[a/25+b/25]++;
        int selected=0;
        for(int pos=0;pos<10;pos++){
            int best=-1,bestScore=-1;
            for(int i=0;i<10;i++)if(!(selected>>i&1)){
                int score=0;
                for(int j=0;j<10;j++)if((edges[i]&edges[j])&&(selected>>j&1))
                    score+=30000/(1+availability[sheet[i]+sheet[j]]);
                if(score>bestScore){best=i;bestScore=score;}
            }
            order[pos]=best;nprior[pos]=0;
            for(int j=0;j<10;j++)if((selected>>j&1)&&(edges[j]&edges[best]))prior[pos][nprior[pos]++]=j;
            selected|=1<<best;
        }
        visit(0,0,0);count++;
    }
    cerr<<count<<" complete problems, "<<nodes<<" partial labelings, "<<hits<<" normalized embeddings\n";
}
