#ifndef JUMP_ARRAY_H_
#define JUMP_ARRAY_H_

#include <vector>
#include <cstdint>

using namespace std;

class JumpArray
{
public:
    JumpArray();
    JumpArray(int64_t size);
    JumpArray(const JumpArray& rhs) = default;

    JumpArray& operator=(const JumpArray& rhs);
    void swap(JumpArray& rhs);

    void delete_index(int64_t index);
    void delete_indices(const vector<int64_t>& indices);

    int64_t size() const { return jumps.size()-1; }
    int64_t space() const { return size() - n_deleted; }

    vector<int64_t> get_indices() const;

    struct JumpIter
    {
        JumpIter(JumpArray& jarr, int64_t idx) : jarr(jarr), idx(idx) {}

        int64_t operator*() const { return idx; }

        JumpIter& operator++()
        {
            idx++;
            idx += jarr.jumps[idx];
            return *this;
        }

        friend bool operator==(const JumpIter& lhs, const JumpIter& rhs) { return lhs.idx == rhs.idx; }
        friend bool operator!=(const JumpIter& lhs, const JumpIter& rhs) { return lhs.idx != rhs.idx; }

    private:
        JumpArray& jarr;
        int64_t idx;
    };

    JumpIter begin() { return JumpIter(*this, jumps.front()); }
    JumpIter end() { return JumpIter(*this, size()); }

private:
    vector<int64_t> jumps;
    int64_t n_deleted;
};

#endif
