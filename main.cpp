#include <iostream>
#include <fstream>
#include <deque>
#include <algorithm>
#include <vector>
#include <chrono>

template<typename T>
typename T::size_type LevenshteinDistance(const T &source, const T &target) {
    if (source.size() > target.size()) {
        return LevenshteinDistance(target, source);
    }

    using TSizeType = typename T::size_type;
    const TSizeType min_size = source.size(), max_size = target.size();
    std::vector<TSizeType> lev_dist(min_size + 1);

    for (TSizeType i = 0; i <= min_size; ++i) {
        lev_dist[i] = i;
    }

    for (TSizeType j = 1; j <= max_size; ++j) {
        TSizeType previous_diagonal = lev_dist[0], previous_diagonal_save;
        ++lev_dist[0];

        for (TSizeType i = 1; i <= min_size; ++i) {
            previous_diagonal_save = lev_dist[i];
            if (source[i - 1] == target[j - 1]) {
                lev_dist[i] = previous_diagonal;
            } else {
                lev_dist[i] = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;
            }
            previous_diagonal = previous_diagonal_save;
        }
    }

    return lev_dist[min_size];
}

class Sequence{
public:
    std::string title;
    std::string source;
    std::int32_t score;

    Sequence(std::int32_t score = 0){
        this->score = score;
        title.reserve(500);
        source.reserve(1'000);
    }
    friend std::istream& operator>>(std::istream& in, Sequence& seq){
        seq.source.clear();
        seq.title.clear();
        std::string buffer;
        char c = in.get();
        seq.title.resize(300);
        in.getline(seq.title.data(), 1000, '\n');
        while(in.peek() != '>' && !in.eof()){
            in >> buffer; in.get();
            seq.source += buffer;
            buffer.clear();
        }
        return in;
    }
    bool operator<(const Sequence& seq){
        return score > seq.score;
    }
};

int main(int argc, char** argv){
    if(argc != 3) exit(1);
    std::ifstream fin_source(argv[1]);
    std::ifstream fin_target(argv[2]);
    std::ofstream fout("result", std::ios::out | std::ios::trunc);
    std::vector<Sequence> top_sequences(100, Sequence(1'000));
    Sequence source, inp;
    fin_source >> source;
    std::uint32_t idx = 1;
    while(!fin_target.eof()){
        fin_target >> inp;
        inp.score = LevenshteinDistance(source.source, inp.source);
        if(inp.score < top_sequences[0].score){
            top_sequences.insert(std::lower_bound(top_sequences.begin(),top_sequences.end(), inp), inp);
            top_sequences.erase(top_sequences.begin());
        }
        std::cout << idx++ << '\n';
    }
    for(auto it = top_sequences.rbegin(); it != top_sequences.rend(); ++it){
        const auto v = *it;
        std::string_view title(v.title);
        if(auto it = std::find(v.title.begin(), v.title.end(),0); it != v.title.end())
            title = std::string_view(v.title.data(), it - v.title.begin());
        fout << '>' << title << '\n' << "score: " << v.score << '\n' << v.source << '\n';
    }
    return 0;
}
