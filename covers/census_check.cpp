// Independently fill each triangular graph, without quotienting by host symmetries.
// Each recursion step assigns a K5 edge (hence a vertex of its line graph).
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using Mask = std::uint64_t;
int difference_id[50][50], endpoints[10], image_of[10], order_of[10];
int profile, problems;
Mask target;
std::vector<int> earlier[10];
long long answers;
void assign(int position, Mask available, Mask unused[2]) {
    if(position == 10) {
        if(available) throw 1;
        std::cout << target << ' ' << profile;
        for(int point : image_of) std::cout << ' ' << point;
        std::cout << '\n'; ++answers;
        return;
    }
    const int vertex=order_of[position], side=(profile >> vertex)&1;
    Mask candidates=position ? unused[side] : Mask(1) << (25*side);
    while(candidates) {
        int p=std::countr_zero(candidates);
        Mask bit=Mask(1)<<p;
        candidates ^= bit;
        Mask left=available;
        bool permitted=true;
        for(int u : earlier[position]) {
            int j=difference_id[p][image_of[u]];
            if(j<0 || !(left & (Mask(1)<<j))) { permitted=false; break; }
            left ^= Mask(1)<<j;
        }
        if(permitted) {
            image_of[vertex]=p; unused[side]^=bit;
            assign(position+1,left,unused);
            unused[side]^=bit;
        }
    }
}
int main() {
    for(auto& row : difference_id) for(int& cell : row) std::cin >> cell;
    int j=0;
    for(int a=0;a<5;a++) for(int b=a+1;b<5;b++) endpoints[j++]=(1<<a)|(1<<b);
    while(std::cin >> target >> profile) {
        int weight[3]={}, used=0;
        for(int side=0;side<3;side++) {
            int origin=side==2?25:0, offset=side==0?0:25;
            for(int x=offset;x<offset+25;x++) {
                int d=difference_id[origin][x];
                if(d>=0 && (target & (Mask(1)<<d))) weight[side]++;
            }
            weight[side]=500/(weight[side]+1);
        }
        for(int pos=0;pos<10;pos++) {
            int best=-1,best_score=-1;
            for(int i=9;i>=0;i--) if(!(used&(1<<i))) {
                int score=0;
                for(int k=0;k<10;k++) if(i!=k && (endpoints[k]&endpoints[i])) {
                    int pairtype=((profile>>i)&1)+((profile>>k)&1);
                    score+= weight[pairtype]*((used>>k&1) ? 100 : 1);
                }
                if(score>best_score) {best=i;best_score=score;}
            }
            order_of[pos]=best; earlier[pos].clear();
            for(int k=0;k<10;k++) if((used&(1<<k)) && (endpoints[k]&endpoints[best])) earlier[pos].push_back(k);
            used |= 1<<best;
        }
        Mask first=(Mask(1)<<25)-1;
        Mask unused[2]={first, first<<25};
        assign(0,target,unused); ++problems;
    }
    std::cerr << problems << " unreduced problems, " << answers << " normalized fillings\n";
}
