template <float_type R, int D> requires is_pos_int<D>
R Traits<R, D>::Distance::operator()(const Point& p, const Point& q)
{
    if constexpr (D == 1) return p < q? q - p : p - q;
    else
    {
        Real dist = 0;
        for (int i = 0; i < D; ++i) dist += (p[i] - q[i])*(p[i] - q[i]);
        return sqrt(dist);
    }
}

template <float_type R, int D> requires is_pos_int<D>
void Traits<R, D>::pack_serialized_point(SerializedPoint& record, const Point& p)
{
    constexpr int dim = D;

    char *dim_dest = record.data();
    char *pt_dest = dim_dest + sizeof(dim);
    char *pt_src;

    if constexpr (D == 1) pt_src = (char*)&p;
    else pt_src = (char*)p.data();

    memcpy(dim_dest, &dim, sizeof(dim));
    memcpy(pt_dest, pt_src, sizeof(p));
}

template <float_type R, int D> requires is_pos_int<D>
void Traits<R, D>::unpack_serialized_point(const SerializedPoint& record, Point& p)
{
    int dim;

    const char *dim_src = record.data();
    const char *pt_src = record.data() + sizeof(int);
    char *pt_dest;

    if constexpr (D == 1) pt_dest = (char*)&p;
    else pt_dest = (char*)p.data();

    memcpy(&dim, dim_src, sizeof(dim));
    assert(dim == D);
    memcpy(pt_dest, pt_src, sizeof(p));
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::write_to_file(Iter first, Iter last, const char *fname)
{
    ofstream os;
    SerializedPoint record;

    os.open(fname, ios::binary | ios::out);

    while (first != last)
    {
        pack_serialized_point(record, *first++);
        os.write(record.data(), sizeof(record));
    }

    os.close();
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::write_to_file(Iter first, Iter last, const char *fname, MPIEnv::Comm comm)
{
    if (comm.is_distributed())
    {
        vector<char> mymem;
        write_to_mem(first, last, mymem);
        comm.file_write_at_all_once(fname, mymem);
    }
    else
    {
        write_to_file(first, last, fname);
    }
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::write_to_mem(Iter first, Iter last, vector<char>& mem)
{
    char *dest;
    size_t num_bytes;
    SerializedPoint record;

    num_bytes = ::distance(first, last) * sizeof(record);
    mem.resize(num_bytes);
    dest = mem.data();

    while (first != last)
    {
        pack_serialized_point(record, *first++);
        memcpy(dest, record.data(), sizeof(record));
        dest += sizeof(record);
    }
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::read_from_mem(Iter d_first, const vector<char>& mem)
{
    Point p;
    size_t n;
    const char *src;
    SerializedPoint record;

    n = mem.size() / sizeof(record);
    src = mem.data();

    for (size_t i = 0; i < n; ++i)
    {
        memcpy(record.data(), src, sizeof(record));
        unpack_serialized_point(record, p); *d_first++ = p;
        src += sizeof(record);
    }
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::read_from_file(Iter d_first, const char *fname)
{
    Point p;
    ifstream is;
    SerializedPoint record;
    size_t filesize, n;

    filesize = get_file_size(fname);
    is.open(fname, ios::binary | ios::in);

    assert(filesize % sizeof(record) == 0);
    n = filesize / sizeof(record);

    for (size_t i = 0; i < n; ++i)
    {
        is.read(record.data(), sizeof(record));
        unpack_serialized_point(record, p);
        *d_first++ = p;
    }

    is.close();
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter>
void Traits<R, D>::read_from_file(Iter d_first, const char *fname, MPIEnv::Comm comm)
{
    if (comm.is_distributed())
    {
        vector<char> mymem;
        comm.file_read_at_all_once(fname, mymem, sizeof(SerializedPoint));
        read_from_mem(d_first, mymem);
    }
    else
    {
        read_from_file(d_first, fname);
    }
}

template <float_type R, int D> requires is_pos_int<D>
template <class RandomGen, class RandomDist>
void Traits<R, D>::fill_random_point(Point& point, RandomGen& gen, RandomDist& dist)
{
    if constexpr (D == 1) point = dist(gen);
    else generate(point.begin(), point.end(), [&gen, &dist] () { return dist(gen); });
}

template <float_type R, int D> requires is_pos_int<D>
template <class Iter, class RandomGen, class RandomDist>
void Traits<R, D>::fill_random_points(Iter first, Iter last, RandomGen& gen, RandomDist& dist)
{
    for_each(first, last, [&gen, &dist](Point& point) { fill_random_point(point, gen, dist); });
}

template <float_type R, int D> requires is_pos_int<D>
string Traits<R, D>::repr(const Point& point, int precision)
{
    stringstream ss;
    ss << fixed << showpoint << showpos << setprecision(precision);

    if constexpr (D == 1) ss << point;
    else
    {
        ss << "[";
        copy(point.begin(), point.end()-1, ostream_iterator<Real>(ss, ","));
        ss << point.back() << "]";
    }

    return ss.str();
}

template <float_type R, int D> requires is_pos_int<D>
size_t Traits<R, D>::PointHash::operator()(const Point& p) const noexcept
{
    if constexpr (same_as<R, Point>) return hash<R>{}(p);
    else
    {
        size_t h = 0;
        for (const R& e : p) h ^= hash<R>{}(e) + (h << 6) + (h >> 2);
        return h;
    }
}
