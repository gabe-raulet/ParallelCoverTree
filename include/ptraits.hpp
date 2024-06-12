template <class Real, int D>
Real Traits<Real, D>::Distance::operator()(const Point& p, const Point& q)
{
    Real dist = 0;
    for (int i = 0; i < D; ++i)
    {
        dist += (p[i] - q[i]) * (p[i] - q[i]);
    }
    return sqrt(dist);
}

template <class Real, int D>
string Traits<Real, D>::name()
{
    stringstream ss;
    ss << "[D=" << D << ",Real=FP" << sizeof(Real)*8 << "]";
    return ss.str();
}

template <class Real, int D>
template <class RandomGen, class RandomDist>
void Traits<Real, D>::fill_random(Point& point, RandomGen& gen, RandomDist& dist)
{
    generate(point.begin(), point.end(), [&gen, &dist] () { return dist(gen); });
}

template <class Real, int D>
template <class RandomGen, class RandomDist>
void Traits<Real, D>::fill_random_vec(vector<Point>& points, RandomGen& gen, RandomDist& dist)
{
    for (Point& p : points) fill_random(p, gen, dist);
}

template <class Real, int D>
void Traits<Real, D>::write_to_file(const vector<Point>& points, const char *outfname)
{
    ofstream os;
    os.open(outfname, ios::binary | ios::out);

    int dim = dimension();

    for (const Point& p : points)
    {
        os.write((char*)&dim, sizeof(int));
        os.write((char*)p.data(), sizeof(Point));
    }

    os.close();
}

template <class Real, int D>
void Traits<Real, D>::read_from_file(vector<Point>& points, const char *infname)
{
    Point p;
    ifstream is;
    int dim = dimension(), d;

    points.clear();
    is.open(infname, ios::binary | ios::in);

    is.seekg(0, is.end);
    size_t filesize = is.tellg();
    is.seekg(0, is.beg);

    size_t n = filesize / (sizeof(Point) + sizeof(int));
    points.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        is.read((char*)&d, sizeof(int));
        assert(d == dim);
        is.read((char*)p.data(), sizeof(Point));
        points.push_back(p);
    }

    is.close();
}
