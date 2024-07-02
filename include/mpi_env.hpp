namespace MPIEnv
{
    int initialize(int *argc, char **argv[])
    {
        if (is_initialized()) return MPI_ERR_OTHER;
        cache.reset();
        return MPI_Init(argc, argv);
    }

    int finalize()
    {
        if (is_finalized()) return MPI_SUCCESS;
        cache.reset();
        return MPI_Finalize();
    }

    void exit(int err)
    {
        finalize();
        ::exit(err);
    }

    bool is_initialized()
    {
        int flag;
        MPI_Initialized(&flag);
        return (flag != 0);
    }

    bool is_finalized()
    {
        int flag;
        MPI_Finalized(&flag);
        return (flag != 0);
    }

    bool comms_equal(MPI_Comm lhs, MPI_Comm rhs)
    {
        int result;
        MPI_Comm_compare(lhs, rhs, &result);
        return (result == MPI_IDENT) || (result == MPI_CONGRUENT);
    }

    bool type_info_compare::operator()(const type_info *lhs, const type_info *rhs) const
    {
        return lhs->before(*rhs);
    }

    template <class T> void mpi_create_type(MPI_Datatype& dtype)
    {
        if constexpr (is_array_type<T>)
        {
            using U = typename array_info<T>::value_type;
            static constexpr int count = array_info<T>::size;
            MPI_Type_contiguous(count, mpi_type<U>(), &dtype);
        }
        else
        {
            MPI_Type_contiguous(sizeof(T), MPI_CHAR, &dtype);
        }

        MPI_Type_commit(&dtype);
    }

    template <class T> MPI_Datatype mpi_commit_type()
    {
        const type_info *t = &typeid(T);
        MPI_Datatype dtype = cache.get_type(t);

        if (dtype == MPI_DATATYPE_NULL)
        {
            mpi_create_type<T>(dtype);
            cache.set_type(t, dtype);
        }

        return dtype;
    }

    template <class T> MPI_Datatype mpi_type()
    {
        if      constexpr (is_same_v<T, char>)               return MPI_CHAR;
        else if constexpr (is_same_v<T, signed char>)        return MPI_SIGNED_CHAR;
        else if constexpr (is_same_v<T, short>)              return MPI_SHORT;
        else if constexpr (is_same_v<T, int>)                return MPI_INT;
        else if constexpr (is_same_v<T, long>)               return MPI_LONG;
        else if constexpr (is_same_v<T, long long>)          return MPI_LONG_LONG;
        else if constexpr (is_same_v<T, unsigned char>)      return MPI_UNSIGNED_CHAR;
        else if constexpr (is_same_v<T, unsigned short>)     return MPI_UNSIGNED_SHORT;
        else if constexpr (is_same_v<T, unsigned int>)       return MPI_UNSIGNED;
        else if constexpr (is_same_v<T, unsigned long>)      return MPI_UNSIGNED_LONG;
        else if constexpr (is_same_v<T, unsigned long long>) return MPI_UNSIGNED_LONG_LONG;
        else if constexpr (is_same_v<T, float>)              return MPI_FLOAT;
        else if constexpr (is_same_v<T, double>)             return MPI_DOUBLE;
        else if constexpr (is_same_v<T, long double>)        return MPI_LONG_DOUBLE;
        else if constexpr (is_same_v<T, bool>)               return MPI_CXX_BOOL;
        else                                                 return mpi_commit_type<T>();
    }

    TypeCache::TypeCache() {}
    TypeCache::~TypeCache() { reset(); }

    MPI_Datatype TypeCache::get_type(const type_info* t)
    {
        auto pos = type_map.find(t);

        if (pos == type_map.end())
            return MPI_DATATYPE_NULL;
        else
            return pos->second;
    }

    void TypeCache::set_type(const type_info *t, MPI_Datatype dtype)
    {
        type_map.insert({t, dtype});
    }

    void TypeCache::free_types()
    {
        if (!is_finalized())
            for (auto& [t, dtype] : type_map)
                MPI_Type_free(&dtype);
    }

    void TypeCache::reset() { free_types(); }

    void Comm::init(MPI_Comm comm)
    {
        MPI_Comm_dup(comm, commbuf);
        MPI_Comm_rank(commbuf[0], &myrank);
        MPI_Comm_size(commbuf[0], &nprocs);
    }

    Comm::Comm() : Comm(MPI_COMM_NULL) {}

    Comm::Comm(MPI_Comm comm) { init(comm); }

    Comm::Comm(const Comm& rhs) { init(rhs.getcomm()); }

    Comm::~Comm() { MPI_Comm_free(commbuf); }

    void Comm::log_at_root_rank(ostream& os, const char *msg, const char *func) const
    {
        if (!myrank)
        {
            char *buf;
            stringstream ss;

            asprintf(&buf, "[msg::%s] :: ", func);
            ss << buf << msg;

            os << ss.str() << endl;
            free(buf);
        }
    }

    void Comm::log_at_root_rank(ostream& os, const char *msg, const char *func, const CommTimer& timer) const
    {
        if (!myrank)
        {
            char *buf;
            stringstream ss;

            asprintf(&buf, "[maxtime=%.3f,avgtime=%.3f,msg::%s] :: ", timer.get_max_time(), timer.get_avg_time(), func);
            ss << buf << msg;

            os << ss.str() << endl;
            free(buf);
        }

    }

    void Comm::log_at_all_ranks(ostream& os, const char *msg, const char *func) const
    {
        char *buf;
        stringstream ss;

        asprintf(&buf, "[myrank=%d,msg::%s] :: ", myrank, func);
        ss << buf << msg << "\n";

        string mystr = ss.str();
        vector<char> mystrbuf(mystr.begin(), mystr.end()), strbuf;

        gatherv(mystrbuf, strbuf, 0);

        if (!myrank)
        {
            string str(strbuf.begin(), strbuf.end());
            os << str << flush;
        }

        free(buf);
    }

    void Comm::log_at_all_ranks(ostream& os, const char *msg, const char *func, const CommTimer& timer) const
    {
        char *buf;
        stringstream ss;

        asprintf(&buf, "[myrank=%d,time=%.3f,msg::%s] :: ", myrank, timer.get_my_time(), func);
        ss << buf << msg << "\n";

        string mystr = ss.str();
        vector<char> mystrbuf(mystr.begin(), mystr.end()), strbuf;

        gatherv(mystrbuf, strbuf, 0);

        if (!myrank)
        {
            free(buf);
            string str(strbuf.begin(), strbuf.end());
            asprintf(&buf, "[maxtime=%.3f,avgtime=%.3f,msg::%s]", timer.get_max_time(), timer.get_avg_time(), func);
            os << str << buf << endl;
        }

        free(buf);
    }

    void Comm::swap(Comm& rhs) noexcept
    {
        MPI_Comm tmp[1];

        ::swap(myrank, rhs.myrank);
        ::swap(nprocs, rhs.nprocs);

        memcpy(tmp, commbuf, sizeof(MPI_Comm));
        memcpy(commbuf, rhs.commbuf, sizeof(MPI_Comm));
        memcpy(rhs.commbuf, tmp, sizeof(MPI_Comm));
    }

    bool Comm::operator==(const Comm &rhs) const
    {
        return comms_equal(getcomm(), rhs.getcomm());
    }

    Comm& Comm::operator=(const Comm& rhs)
    {
        init(rhs.getcomm());
        return *this;
    }

    int Comm::barrier() const
    {
        return MPI_Barrier(getcomm());
    }

    template <class T>
    bool Comm::is_same_val(T val) const
    {
        auto comm = getcomm();
        bool is_same;

        if (nprocs == 1) return true;

        auto dtype = mpi_type<T>();

        if constexpr (is_arithmetic_v<T>)
        {
            T buf[2] = {-val, val};
            MPI_Allreduce(MPI_IN_PLACE, buf, 2, dtype, MPI_MIN, comm);
            is_same = (buf[0] == -buf[1]);
        }
        else
        {
            vector<T> buffer;
            if (myrank == 0) buffer.resize(nprocs);

            MPI_Gather(&val, 1, dtype, buffer.data(), 1, dtype, 0, comm);

            if (myrank == 0) is_same = equal(buffer.begin()+1, buffer.end(), buffer.begin());

            MPI_Bcast(&is_same, 1, mpi_type<bool>(), 0, comm);
        }

        return is_same;
    }

    template <class T>
    bool Comm::are_same_vals(const vector<T>& vals) const
    {
        auto comm = getcomm();

        if (nprocs == 1) return true;

        int count = vals.size();

        if (!is_same_val(count)) return false;

        if (count == 0) return true;

        auto dtype = mpi_type<T>();

        vector<T> buf = vals, rbuf;

        if constexpr (is_arithmetic_v<T>)
        {
            buf.resize(count<<1);
            transform(buf.begin(), buf.begin() + count, buf.begin() + count, [](T x) { return -x; });

            MPI_Allreduce(MPI_IN_PLACE, buf.data(), count<<1, dtype, MPI_MIN, comm);

            for (int i = 0; i < count; ++i)
                if (-buf[i] != buf[i+count])
                    return false;

            return true;
        }
        else
        {
            if (myrank == 0) rbuf.resize(count*nprocs);

            MPI_Gather(buf.data(), count, dtype, rbuf.data(), count, dtype, 0, comm);

            bool are_same = true;

            if (myrank == 0)
            {
                auto it = rbuf.begin() + count;

                for (int i = 1; are_same && i < nprocs; it += count, ++i)
                    if (!equal(buf.begin(), buf.end(), it))
                        are_same = false;
            }

            MPI_Bcast(&are_same, 1, mpi_type<bool>(), 0, comm);
            return are_same;
        }
    }

    template <class T>
    int Comm::reduce(const T* sendbuf, T* recvbuf, int count, int root, MPI_Op op) const
    {
        return MPI_Reduce(sendbuf, recvbuf, count, mpi_type<T>(), op, root, getcomm());
    }

    template <class T>
    int Comm::reduce(const T& sendbuf, T& recvbuf, int root, MPI_Op op) const
    {
        return reduce(&sendbuf, &recvbuf, 1, root, op);
    }

    template <class T>
    int Comm::reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, MPI_Op op) const
    {
        int count = sendbuf.size();
        if (!is_same_val(count)) return MPI_ERR_OTHER;
        if (myrank == root) recvbuf.resize(count);
        return reduce(sendbuf.data(), recvbuf.data(), count, root, op);
    }

    template <class T>
    int Comm::reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, const vector<MPI_Op>& ops) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();
        int count = sendbuf.size();
        vector<int> counts(2); counts[0] = count, counts[1] = ops.size();

        if (!are_same_vals(counts) || count != ops.size())
            return MPI_ERR_OTHER;

        if (myrank == root) recvbuf.resize(count);

        for (int i = 0; i < count; ++i)
        {
            int err = MPI_Reduce(&sendbuf[i], &recvbuf[i], 1, dtype, ops[i], root, comm);
            if (err != MPI_SUCCESS) return err;
        }

        return MPI_SUCCESS;
    }

    template <class T>
    int Comm::bcast(T* buffer, int count, int root) const
    {
        return MPI_Bcast(buffer, count, mpi_type<T>(), root, getcomm());
    }

    template <class T>
    int Comm::bcast(T& buffer, int root) const
    {
        return bcast(&buffer, 1, root);
    }

    template <class T>
    int Comm::bcast(vector<T>& buffer, int root) const
    {
        return bcast(buffer.data(), static_cast<int>(buffer.size()), root);
    }

    template <class T>
    int Comm::exscan(const T* sendbuf, T* recvbuf, int count, MPI_Op op, T identity) const
    {
        fill(recvbuf, recvbuf + count, identity);
        return MPI_Exscan(sendbuf, recvbuf, count, mpi_type<T>(), op, getcomm());
    }

    template <class T>
    int Comm::exscan(const T& sendbuf, T& recvbuf, MPI_Op op, T identity) const
    {
        return exscan(&sendbuf, &recvbuf, 1, op, identity);
    }

    template <class T> int Comm::allreduce(const T* sendbuf, T* recvbuf, int count, MPI_Op op) const
    {
        return MPI_Allreduce(sendbuf, recvbuf, count, mpi_type<T>(), op, getcomm());
    }

    template <class T> int Comm::allreduce(const T& sendbuf, T& recvbuf, MPI_Op op) const
    {
        return allreduce(&sendbuf, &recvbuf, 1, op);
    }

    template <class T> int Comm::allreduce(T* buffer, int count, MPI_Op op) const
    {
        return MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_type<T>(), op, getcomm());
    }
    template <class T> int Comm::allreduce(T& buffer, MPI_Op op) const
    {
        return allreduce(&buffer, 1, op);
    }

    template <class T> int Comm::gatherv(const vector<T>& sendbuf, vector<T>& recvbuf, int root) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int sendcount = sendbuf.size();
        vector<int> recvcounts, displs;

        if (myrank == root)
        {
            recvcounts.resize(nprocs);
            displs.resize(nprocs);
        }

        int err = MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, root, comm);
        if (err != MPI_SUCCESS) return err;

        if (myrank == root)
        {
            exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), static_cast<int>(0));
            recvbuf.resize(recvcounts.back() + displs.back());
        }

        return MPI_Gatherv(sendbuf.data(), sendcount, dtype, recvbuf.data(), recvcounts.data(), displs.data(), dtype, root, comm);
    }

    template <class T> int Comm::allgatherv(const vector<T>& sendbuf, vector<T>& recvbuf) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        vector<int> recvcounts(nprocs), displs(nprocs);

        recvcounts[myrank] = sendbuf.size();

        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
        if (err != MPI_SUCCESS) return err;

        exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), static_cast<int>(0));
        recvbuf.resize(recvcounts.back() + displs.back());

        return MPI_Allgatherv(sendbuf.data(), recvcounts[myrank], dtype, recvbuf.data(), recvcounts.data(), displs.data(), dtype, comm);
    }

    template <class T> int Comm::scatterv(const vector<T>& sendbuf, const vector<int>& sendcounts, vector<T>& recvbuf, int root) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int valid = !!(sendcounts.size() == nprocs);
        MPI_Bcast(&valid, 1, MPI_INT, root, comm);

        if (!valid) return MPI_ERR_OTHER;

        int recvcount;
        vector<int> displs;

        if (myrank == root)
        {
            displs.resize(nprocs);
            exclusive_scan(sendcounts.begin(), sendcounts.end(), displs.begin(), static_cast<int>(0));
        }

        int err = MPI_Scatter(sendcounts.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, root, comm);
        if (err != MPI_SUCCESS) return err;

        recvbuf.resize(recvcount);

        return MPI_Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), dtype, recvbuf.data(), recvcount, dtype, root, comm);
    }

    int Comm::file_write_at_all_once(const char *fname, const char *mybuf, int count, bool truncate) const
    {
        MPI_Offset mysize = count;
        MPI_Offset fileoffset;

        exscan(mysize, fileoffset, MPI_SUM, static_cast<MPI_Offset>(0));

        MPI_File fh;
        MPI_File_open(getcomm(), fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        if (truncate)
        {
            MPI_Offset filesize;
            MPI_File_get_size(fh, &filesize);
            truncate = (filesize > 0);
            bcast(truncate, 0);
            if (truncate) MPI_File_set_size(fh, 0);
        }

        MPI_File_write_at_all(fh, fileoffset, mybuf, count, MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);

        return MPI_SUCCESS;
    }

    template <string_type String>
    int Comm::file_write_at_all_once(const char *fname, const String& mybuf, bool truncate) const
    {
        return file_write_at_all_once(fname, mybuf.data(), static_cast<int>(mybuf.size()), truncate);
    }

    template <string_type String>
    int Comm::file_read_at_all_once(const char *fname, String& mybuf, int unitsize) const
    {
        MPI_File fh;
        MPI_File_open(getcomm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

        MPI_Offset filesize, fileoffset;
        MPI_File_get_size(fh, &filesize);

        assert(filesize % unitsize == 0);
        size_t numunits = filesize / unitsize;

        vector<MPI_Offset> counts(nprocs), displs(nprocs);
        get_balanced_counts(counts, numunits);
        mybuf.resize(counts[myrank]*unitsize);

        exclusive_scan(counts.begin(), counts.end(), displs.begin(), static_cast<MPI_Offset>(0));
        fileoffset = displs[myrank] * unitsize;

        MPI_File_read_at_all(fh, fileoffset, mybuf.data(), static_cast<int>(mybuf.size()), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);

        return MPI_SUCCESS;
    }

    template <class Real, class Index>
    void ArgmaxPair<Real, Index>::mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype)
    {
        ArgmaxPair *in = (ArgmaxPair*)_in;
        ArgmaxPair *inout = (ArgmaxPair*)_inout;

        for (int i = 0; i < *len; ++i)
            if (inout[i].value < in[i].value)
            {
                inout[i].value = in[i].value;
                inout[i].index = in[i].index;
            }
    }

    template <class Real, class Index>
    void ArgmaxPair<Real, Index>::create_mpi_handlers(MPI_Datatype& MPI_ARGMAX_PAIR, MPI_Op& MPI_ARGMAX)
    {
        int blklens[2] = {1,1};
        MPI_Aint disps[2] = {offsetof(ArgmaxPair, index), offsetof(ArgmaxPair, value)};
        MPI_Datatype types[2] = {mpi_type<Index>(), mpi_type<Real>()};
        MPI_Type_create_struct(2, blklens, disps, types, &MPI_ARGMAX_PAIR);
        MPI_Type_commit(&MPI_ARGMAX_PAIR);
        MPI_Op_create(&mpi_argmax, 1, &MPI_ARGMAX);
    }
}
