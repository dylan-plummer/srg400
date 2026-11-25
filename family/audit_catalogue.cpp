// Rebuild every entry from the original graph6 file, without the starter or Python code.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <openssl/sha.h>
extern "C" {
#include <nauty.h>
}
using Row = std::array<uint64_t,7>;
using Matrix = std::array<Row,400>;
struct Block {int p,q,c,d,negative;};

void set_edge(Matrix &a,int i,int j) {
    a[i][j/64] |= uint64_t(1)<<(j%64);
    a[j][i/64] |= uint64_t(1)<<(i%64);
}
int has_edge(const Matrix &a,int i,int j) {
    return (a[i][j/64]>>(j%64))&1;
}
Matrix read_seed(const std::string &name) {
    std::ifstream in(name);
    std::string line;
    std::getline(in,line);
    assert(line.size()==13304 && line.substr(0,4)=="~?EO");
    Matrix a{};
    int t=0;
    for(int j=1;j<400;j++) for(int i=0;i<j;i++,t++)
        if(((line[4+t/6]-63)>>(5-t%6))&1) set_edge(a,i,j);
    return a;
}
std::vector<Block> read_blocks(const Matrix &seed) {
    std::vector<Block> blocks;
    for(int p=0;p<50;p++) for(int q=p+1;q<50;q++) {
        int base=has_edge(seed,8*p,8*q), c=0,d=0,weight=0;
        for(int x=0;x<8;x++) weight+=has_edge(seed,8*p+x,8*q);
        if(weight==0) continue;
        assert(weight==4);
        for(int j=0;j<3;j++) {
            c |= (base^has_edge(seed,8*p+(1<<j),8*q))<<j;
            d |= (base^has_edge(seed,8*p,8*q+(1<<j)))<<j;
        }
        assert(c&&d);
        blocks.push_back({p,q,c,d,1-base});
    }
    assert(blocks.size()==1050);
    return blocks;
}
Matrix rebuild(const std::vector<Block> &blocks,const std::string &chars,int signs) {
    Matrix a{};
    for(auto v:blocks) {
        int first=7*(v.p/25)+v.c-1, second=7*(v.q/25)+v.d-1;
        int sign=v.negative^((signs>>first)&1)^((signs>>second)&1);
        int c=chars[first]-'0',d=chars[second]-'0';
        for(int x=0;x<8;x++) for(int y=0;y<8;y++)
            if((__builtin_parity(x&c)^__builtin_parity(y&d))==sign)
                set_edge(a,8*v.p+x,8*v.q+y);
    }
    return a;
}
void verify(const Matrix &a) {
    for(int i=0;i<400;i++) for(int j=i;j<400;j++) {
        int n=0;
        for(int w=0;w<7;w++) n+=__builtin_popcountll(a[i][w]&a[j][w]);
        assert(n==(i==j?168:(has_edge(a,i,j)?68:72)));
    }
}
std::string canonical_digest(const Matrix &a,int expected_group) {
    constexpr int n=400, words=SETWORDSNEEDED(n);
    std::array<setword,n*words> g{},can{};
    std::array<int,n> lab{},partition{},orbits{};
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon=TRUE;
    statsblk stats{};
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
        if(has_edge(a,i,j)) ADDONEEDGE(g.data(),i,j,words);
    densenauty(g.data(),lab.data(),partition.data(),orbits.data(),&options,&stats,words,n,can.data());
    assert(stats.errstatus==0 && stats.grpsize1==expected_group && stats.grpsize2==0);
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(can.data()),sizeof(can),digest);
    std::ostringstream out;
    for(unsigned char v:digest) out<<std::hex<<std::setfill('0')<<std::setw(2)<<int(v);
    return out.str();
}
int main(int argc,char **argv) {
    if(argc!=5) {std::cerr<<"usage: audit_catalogue seed.g6 catalogue.tsv first past_last\n"; return 2;}
    Matrix seed=read_seed(argv[1]);
    auto blocks=read_blocks(seed);
    int first=std::stoi(argv[3]),end=std::stoi(argv[4]);
    std::ifstream in(argv[2]);std::string line;std::getline(in,line);
    int id,order,mask,count=0;std::string p,q;
    while(in>>id>>order>>p>>q>>mask) if(id>=first&&id<end) {
        Matrix a=rebuild(blocks,p+q,mask);
        if(id==0) assert(a==seed);
        verify(a);
        std::cout<<id<<' '<<order<<' '<<canonical_digest(a,order)<<'\n';
        count++;
    }
    assert(count==end-first);
}
