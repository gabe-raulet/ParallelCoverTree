#include "JumpArray.h"
#include <utility>

JumpArray::JumpArray() : JumpArray(0) {}

JumpArray::JumpArray(int64_t size) : jumps(size+1, 0), n_deleted(0) {}

JumpArray& JumpArray::operator=(const JumpArray& rhs)
{
    JumpArray tmp(rhs);
    tmp.swap(*this);
    return *this;
}

void JumpArray::swap(JumpArray& rhs)
{
    ::swap(jumps, rhs.jumps);
    ::swap(n_deleted, rhs.n_deleted);
}

void JumpArray::delete_index(int64_t index)
{
    if (index < 0 || index >= size() || !!jumps[index])
        return;

    jumps[index] += (jumps[index+1]+1);
    int64_t v = jumps[index];
    index--;


    while (index >= 0 && !!jumps[index])
    {
        jumps[index] = ++v;
        index--;
    }

    n_deleted++;
}

void JumpArray::delete_indices(const vector<int64_t>& indices)
{
    for (int64_t index : indices)
        delete_index(index);
}

vector<int64_t> JumpArray::get_indices() const
{
    vector<int64_t> indices;

    int64_t i = 0;

    while (true)
    {
        i += jumps[i];
        if (i >= size()) break;
        indices.push_back(i);
        i++;
    }

    return indices;
}
