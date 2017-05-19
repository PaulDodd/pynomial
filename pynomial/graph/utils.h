//
// utils.h
//
#ifndef SHARED_FILE_UTILS_H
#define SHARED_FILE_UTILS_H

// #include "SharedInclude.h"
#ifdef c_plus_plus_11
#include <chrono>
#include <random>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace utils
{

using namespace std;
typedef long long int64;

#define MIN4(a, b, c, d)                 min(min(a,b), min(c,d))
#define BOUNDED(x, xmin, xmax)          (x >= xmin && x <= xmax)


#ifdef c_plus_plus_11
namespace rng_util{

template<class RNG>
inline void seed_generator(RNG& generator, const size_t& n = 100)
{
    vector<size_t> seeds;
    try {
        random_device rd;
        for(size_t i = 0; i < n; i++)
            seeds.push_back(rd());
    } catch (...) {
        cout << "random_device is not available..." << endl;
        seeds.push_back(size_t(std::chrono::system_clock::now().time_since_epoch().count()));
        seeds.push_back(size_t(getpid()));
    }
//    cout << "seeds = ";
//    for(size_t i = 0; i < n; i++)
//        cout << seeds[i] << " ";
//    cout << endl;

    seed_seq seq(seeds.begin(), seeds.end());
    generator.seed(seq);
}

}
#endif

namespace math_util{

template<class TVal>
inline vector<TVal> vector_range(TVal first = 0, TVal last = 0, TVal inc = 1)
{
    vector<TVal> v;
    for(TVal x = first; x < last; x+=inc) v.push_back(x);
    return v;
}

template<typename TVal >
inline vector<TVal> Intersection(const vector<TVal>& v1, const vector<TVal>& v2)
{
    assert(std::is_sorted(v1.begin(), v1.end()) && std::is_sorted(v2.begin(), v2.end()));
    vector<TVal> intersection(min(v1.size(), v2.size()));
    less<TVal> comp; // TODO: extract this.
    typename vector<TVal>::iterator it = std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), intersection.begin(), comp);
    intersection.resize(it-intersection.begin());
    return intersection;
}

template<class TVal>
inline vector<TVal> Union(const vector<TVal>& v1, const vector<TVal>& v2)
{
    assert(std::is_sorted(v1.begin(), v1.end()) && std::is_sorted(v2.begin(), v2.end()));
    vector<TVal> u(v1.size()+v2.size());
    typename vector<TVal>::iterator it = std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), u.begin());
    u.resize(it-u.begin());
    return u;
}


template<class ListType, class RelationType>
inline std::vector<size_t> RelationClasses(const ListType& set, const RelationType& R)
{
    // Computational Complexity: O(N^2 + N*P) N = set.size(),  P = # of partitions
    // Memory Complexity: O(N)

    std::vector<size_t> partitions(vector_range<size_t>(0, set.size())); // move constructor makes this not so bad
    for(size_t i = 0; i < set.size(); i++)
    {
        for(size_t j = i+1; j < set.size(); j++)
        {
            if(R(set[i], set[j]))
            {
                partitions[j] = partitions[i];
            }
        }
    }

    std::set<size_t> s(partitions.begin(), partitions.end());
    size_t ct = 0;
    for(std::set<size_t>::iterator iter = s.begin(); iter != s.end(); iter++)
    {
        if(ct == *iter)
        {
            ct++;
            continue;
        }

        for(size_t i = 0; i < partitions.size(); i++)
        {
            if(partitions[i] == *iter)
            {
                partitions[i] = ct;
            }
        }
        ct++;
    }

    return partitions;
};

#ifdef c_plus_plus_11
inline constexpr size_t invalid_id(){ return size_t(-1); }

inline std::vector< std::vector<size_t> > PartionSets(const std::vector<size_t>& p)
{
    // Assuming that the partition ids in p are 0 based and consecutive
    // Computational Complexity: O(N)
    // Memory Complexity: O(N)

    std::vector< std::vector<size_t> > sets;
    for(size_t i = 0; i < p.size(); i++)
    {
        if(!(p[i] < sets.size()))
            while(!(p[i] < sets.size())) sets.push_back(vector<size_t>());
        sets[p[i]].push_back(i);
    }
    return sets;
}


inline std::vector<size_t> PartitionUnion(const std::vector<size_t>& p1 , const std::vector<size_t>& p2)
{
    // Computational Complexity: O(N^2) worst case.
    // Memory Complexity: O(N)
    assert(p1.size() == p2.size());
    std::vector<size_t> partition(p1.size(), invalid_id());
    std::vector< std::vector<size_t> > sets1 = PartionSets(p1), sets2 = PartionSets(p2); // O(N)
    size_t pid = 0;

    for(size_t i = 0; i < p1.size(); i++)
    {
        std::queue<size_t> worklist1, worklist2;
        if(partition[i] == invalid_id()) // seeing this partition for the first time.
        {
            worklist1.push(i);
            worklist2.push(i);
            partition[i] = pid++;
        }

        while(!worklist1.empty() || !worklist2.empty())
        {
            if(!worklist1.empty())
            {
                size_t x = worklist1.front(); worklist1.pop();
                for(size_t y : sets1[p1[x]])
                {
                    if(partition[y] == invalid_id())
                    {
                       partition[y] = partition[i];
                       worklist2.push(y);
                    }
                    else if (partition[y] < partition[i] || partition[i] < partition[y])
                    {
                        cout << "logical error!" << endl;
                    }
                }
            }
            if(!worklist2.empty())
            {
                size_t x = worklist2.front(); worklist2.pop();
                for(size_t y : sets2[p2[x]])
                {
                    if(partition[y] == invalid_id())
                    {
                        partition[y] = partition[i];
                        worklist1.push(y);
                    }
                    else if (partition[y] < partition[i] || partition[i] < partition[y])
                    {
                        cout << "logical error!" << endl;
                    }
                }
            }
        }
    }

    return partition;
}

template<class ReverseIterator>
inline bool next_combination(ReverseIterator first, ReverseIterator last, int N)
{
    //    returns next combination.
    if(N == 0 || first == last)
        return false;

    int i = 0;
    ReverseIterator begin = first;
    while(first != last && i < N)
    {
        ReverseIterator it = first++;
        if((++(*it)) < (N - i)) {
            while(it != begin) {
                first = it--;
                *it = *first + 1;
            }
            return true;
        }
        i++;
    }
    return false;
}

#endif
}



#ifdef c_plus_plus_11
class TrueFunction
{
    public:
        template<typename... TVals>
        bool operator()(TVals...){ return true; }
};

class FalseFunction
{
    public:
        template<typename... TVals>
        bool operator()(TVals...){ return false; }
};

template<class ListType, class ReturnType > // put here for template specialization.
class access
{
    template<class SubListType, class... IndexList>
    typename enable_if<sizeof...(IndexList) == 0,  ReturnType >::type& get_ith(SubListType& list, size_t I, IndexList...) const
    {
        return list[I];
    }

    template<class SubListType, class... IndexList>
    typename enable_if< 0 < sizeof...(IndexList),  ReturnType >::type& get_ith(SubListType& list, size_t I, IndexList... Others) const
    {
        return get_ith(list[I], Others...);
    }

    public:
//        template<class... IndexList>
//        ReturnType& operator () (ListType& list, size_t I, IndexList... Others) const { return get_ith< ReturnType >(list, I, Others...); }

        template<class... IndexList>
        const ReturnType& operator () (ListType& list, size_t I, IndexList... Others) const {  return get_ith(list, I, Others...); }
};

#else
class TrueFunction
{
    public:
        template<typename TVal>
        bool operator()(TVal, TVal){ return true; }
};


#endif


/************************************************************************************************************************************************************************/
// std::complex class helper functions.
/************************************************************************************************************************************************************************/
template < class TVal>
class complex_greater_than : binary_function<std::complex<TVal>, std::complex<TVal>, bool>
{
    public:
        bool operator()(const std::complex<TVal>& x, const std::complex<TVal>& y ) { return std::norm(x) > std::norm(y); }
};

/************************************************************************************************************************************************************************/
// std::vector class helper functions.
/************************************************************************************************************************************************************************/
#ifdef c_plus_plus_11
template<class TVal>
inline void vector_copy_assign(vector< vector<TVal> >& dest, const vector< vector<TVal> >& src)
{
    dest.resize(src.size());
    for(size_t i = 0; i < src.size(); i++)
    {
        dest[i].resize(src[i].size());
        for(size_t j = 0; j < src[i].size(); j++)
        {
            const TVal& element = src[i][j];
            dest[i][j] = element;
        }
    }
}

template<class TVal>
inline double mean(const vector<TVal>& v)
{
    double sum = 0;
    for(TVal x : v)
        sum+=x;
    return v.size() ? sum/v.size() : 0;
}
template<class TVal>
inline pair<double, double> mean_std_dev(const vector<TVal>& v)
{
    double sum = 0;
    double m = mean(v), d = 0.0;
    for(TVal x : v){
        d = x - m;
        sum+=d*d;
    }
    return pair<double, double> (m, ((v.size() > 1) ? sqrt(sum/(v.size()-1)) : 0));
}

template<class TVal>
inline double std_dev(const vector<TVal>& v)
{
    return mean_std_dev(v).second;
}
#endif


template <class ItemType, class Compare = less<ItemType> >
class circular_list_equal_to
{
    Compare     m_less_than;

    template<class ListType>
    inline bool COMP(const ListType& A, const ListType& B, const size_t& n, size_t& i, size_t& j)
    {
        bool bComp = true;
        for(size_t k = 0; k < n && bComp; k++) {
            if(m_less_than(A[(i+k)%n], B[(j+k)%n])) {
                j += k+1;
                bComp = false;
            }
            else if(m_less_than(B[(j+k)%n], A[(i+k)%n])) {
                i += k+1;
                bComp = false;
            }
        }
        return bComp;
    }

    public:
        circular_list_equal_to(Compare lt = Compare()) : m_less_than(lt){ }

        template<class ListType>
        bool operator()(const ListType& A, const ListType& B, const size_t& n)
        {
            bool bEqual = false, bContinue = true;
            size_t i = 0, j = 0;
            vector<size_t> DA, DB;
            while( bContinue ) {
                bEqual = COMP( A, B, n, i, j);
                bContinue = i < n && j < n && !bEqual;
            }
            return bEqual;
        }
};

typedef circular_list_equal_to<char> circular_string_equal_to;


//inline bool CircularCompareString(const string& A, const string& B)
//{
//    bool bEqual = false, bContinue = A.length() == B.length();
//    size_t i = 0, j = 0, n = A.length();
//    vector<size_t> DA, DB;
//    while( bContinue )
//    {
//        bEqual = COMP(A, B, A.length(), i, j);
//        bContinue = i < n && j < n && !bEqual;
//    }
//    return bEqual;
//}


template<class TVal>
inline size_t argmin(const vector<TVal>& v)
{
//    if(!v.size())
//        throw runtime_error("Empty List encountered in utils::argmin.");
    TVal m = v[0];
    size_t mndx = 0;
    for(size_t i = 1; i < v.size(); i++){
        if(v[i] < m){
            m = v[i];
            mndx = i; } }
    return mndx;
}

template<class TVal>
inline TVal mininum(const vector<TVal>& v)
{
//    if(!v.size())
//        throw runtime_error("Empty List encountered in utils::minimum.");
    return v[argmin(v)];
}

template<class TVal>
inline size_t argmax(const vector<TVal>& v)
{
//    if(!v.size())
//        throw runtime_error("Empty List encountered in utils::argmax.");
    if(!v.size())
        return 0;

    TVal m = v[0];
    size_t mndx = 0;
    for(size_t i = 1; i < v.size(); i++){
        if(v[i] > m){
            m = v[i];
            mndx = i; } }
    return mndx;
}

template<class TVal>
inline TVal maximum(const vector<TVal>& v)
{
//    if(!v.size())
//        throw runtime_error("Empty List encountered in utils::maximum.");
    if(!v.size())
        return TVal();
    return v[argmax(v)];
}

template<class TVal>
inline void FromArrayToVector(vector<TVal>& vec, const TVal * pArray, const size_t& szArray)
{
    assert(pArray);
    for( int i = 0; i < szArray && pArray; i ++ )
    {
        vec.push_back(pArray[i]);
    }
}

template<class TVal>
inline void FromVectorToArray(const vector<TVal>& vec, TVal * pArray)
{
    assert(pArray);
    for( int i = 0; i < vec.size() && pArray; i ++ )
    {
        pArray[i] = vec[i];
    }
}

template<class TVal>
inline vector<TVal> VectorCat(const vector<TVal>& v1, const vector<TVal>& v2)
{
    vector<TVal> ret = v1;
    for(size_t i = 0; i < v2.size(); i++)
        ret.push_back(v2[i]);
    return ret;
}

template<class TVal>
inline void VectorAppend(vector<TVal>& v1, const vector<TVal>& v2)
{
    for(size_t i = 0; i < v2.size(); i++)
        v1.push_back(v2[i]);
}
template<class TVal>
inline vector<TVal> VectorSlice(const vector<TVal>& v, const vector<size_t> i)
{
    vector<TVal> ret;
    for(size_t s = 0; s < i.size(); s++)
        ret.push_back(v[i[s]]);
    return ret;
}

#ifdef c_plus_plus_11
template<typename TVal, typename Compare = equal_to<TVal> >
inline size_t FindInVec(const vector<TVal>& v, const TVal& val)
{
    bool bFound = false;
    size_t ndx = v.size();
    Compare IsEqual;
    for(size_t i = 0; i < v.size() && !bFound; i++)
    {
        if(IsEqual(val, v[i]))
        {
            bFound = true;
            ndx = i;
        }
    }
    return ndx;
}

template<typename TVal, typename Compare = equal_to<TVal> >
inline bool IsInVec(const vector<TVal>& v, const TVal& val)
{
    return (FindInVec<TVal, Compare>(v, val) < v.size());
}

template<typename TVal, typename Compare = equal_to<TVal> >
inline bool PushUnique(vector<TVal>& v, const TVal& elem)
{
    bool bPushedElem = false;
    if(!IsInVec<TVal, Compare>(v, elem))
    {
        v.push_back(elem);
        bPushedElem = true;
    }

    return bPushedElem;
}


#else

template< typename TVal >
inline size_t FindInVec(const vector<TVal>& v, const TVal& val)
{
    bool bFound = false;
    size_t ndx = v.size();
    equal_to<TVal> IsEqual;
    for(size_t i = 0; i < v.size() && !bFound; i++)
    {
        if(IsEqual(val, v[i]))
        {
            bFound = true;
            ndx = i;
        }
    }
    return ndx;
}

template<typename TVal >
inline bool IsInVec(const vector<TVal>& v, const TVal& val)
{
    return (FindInVec<TVal>(v, val) < v.size());
}

template<typename TVal >
inline bool PushUnique(vector<TVal>& v, const TVal& elem)
{
    bool bPushedElem = false;
    if(!IsInVec<TVal>(v, elem))
    {
        v.push_back(elem);
        bPushedElem = true;
    }

    return bPushedElem;
}

template<typename TVal >
inline vector<TVal> Intersection(const vector<TVal>& v1, const vector<TVal>& v2)
{
    vector<TVal> intersection;

    for(size_t i = 0; i < v1.size(); i++)
    {
        if(IsInVec<TVal>(v2, v1[i]))
            PushUnique(intersection, v1[i]); // no repeats allowed.
    }

    return intersection;
}

#endif

template<class TVal>
inline TVal index_sum(const vector<TVal> v, const vector<size_t>& ndx)
{
    TVal sum = TVal(0);
    for(vector<size_t>::const_iterator i = ndx.begin(); i != ndx.end(); i++)
        sum += v[*i];
    return sum;
}


template< class Key, class Value >
inline std::map<Value, std::vector<Key> > preimage(const std::map<Key, Value>& m)
{
    std::map<Value, std::vector<Key> > preim;
    typename std::map<Key,Value>::const_iterator miter = m.cbegin();
    for(; miter != m.cend(); miter++)
    {
        std::pair<typename std::map<Value, std::vector<Key> >::iterator, bool> fiter = preim.insert(pair<Value, std::vector<Key> >(miter->second, std::vector<Key>()));
        fiter.first->second.push_back(miter->first);
    }
    return preim;
}



/************************************************************************************************************************************************************************/
// Some helpful math functions
/************************************************************************************************************************************************************************/

// sign(x): returns +1 or -1 if val is greater or less that 0 respectively.
// modified this to also return 0 if fabs(val) < tol.
template<class TVal>
inline int sign(const TVal& val, const double tol = 0.0)
{
    assert(tol >= 0.0);
    return ( fabs(val) <= tol ? 0 : ((val < 0) ? -1 : 1));
}

inline int mod_dist(const int& a, const int& b, const int& n)
{
    int an = abs(a % n), bn = abs(b % n);
    int d = abs(an-bn);
    if(abs(d) > (n/2))
    {
        d = n - d;
    }
    return d;
}

template< class TVal>
inline TVal arccos(const TVal& x){
    if(x < -1)
        return acos(-1);
    else if(x > 1)
        return acos(1);
    else
        return acos(x);
}


// this is strange but we have to do it this way because floating point template parameters are not supported.
// tol = constant * 10^power;
template<class TVal, int power = -6, size_t constant = 1>
struct float_is_equal : binary_function<TVal, TVal, bool>
{
    TVal tol;
    float_is_equal() : tol(TVal(constant)*pow(TVal(10.0), TVal(power))) {}
    bool operator()(const TVal& x, const TVal& y) { return fabs(x-y) < tol;}
};

template<typename TVal, int power = -6, size_t constant = 1>
struct float_vec_is_equal : binary_function<vector<TVal>, vector<TVal>, bool>
{
    float_is_equal<TVal, power, constant> f;
    float_vec_is_equal() {}
    bool operator()(const vector<TVal>& x, const vector<TVal>& y)
    {
        bool bEqual = true;
        //cout << "x.size:"<< x.size() << " y.size():"<< y.size()<<endl;
        if(x.size() != y.size()) bEqual = false;
        for(size_t i = 0; i < x.size() && bEqual; i++) bEqual = f(x[i], y[i]);
        return bEqual;
    }
};



/************************************************************************************************************************************************************************/
// std::string class helper functions.
/************************************************************************************************************************************************************************/


template<class TVal>
inline string NumberToString(TVal num, const size_t& szFill = 0, const char& chFill = '0', const std::ios_base::fmtflags setflags = std::ios_base::dec, const std::ios_base::fmtflags unsetflags = std::ios_base::fmtflags(0)) // TODO: add formatting options.
{
    // could alternatively use sprintf.
    string str;
    stringstream stream;
    stream.unsetf(unsetflags);
    stream.setf(setflags);
    if(szFill > 0)
        stream << setfill(chFill) << setw(int(szFill));
    stream<<num;
    stream>>str;
    return str;
}

inline string StripWhiteSpace(const string& str)
{
    // strips the white space only at the begining or end of a string.
    string ret = str;
    size_t back = ret.length() > 0 ? ret.length()-1 : std::string::npos;
    while(std::isspace(ret[0]) || std::isspace(ret[back]))
    {
        if(std::isspace(ret[0]))
            ret.erase(0,1);
        if(std::isspace(ret[back]))
            ret.erase(back,1);
        back = ret.length() > 0 ? ret.length()-1 : std::string::npos;
    }
    return ret;
}

inline vector<string> SplitString(const string& str, const string& splitStr, bool bStripWhiteSpace = true, bool bAllowSuccessiveDelim = false)
{
    vector<string> strvec;
    size_t pos = 0; size_t prevPos = 0;
    string temp;
    bool bMore = (str.length() != 0);

    if(str.length() == 1)
    {
        strvec.push_back((str == splitStr) ? "" : str);
        return strvec;
    }
    else if(splitStr.length() == 0)
    {
        cout << "ERROR! Must use a valid delimiter. " << endl;
        throw(0xfffff);
    }

    while (bMore) {
        pos = str.find(splitStr, prevPos);

        if ( pos == string::npos)
        {
            pos = str.length();
            //cout << "string: "<< str << " @  ("<< prevPos << ", " << pos <<")" << endl;
            bMore = false;
        }
        if((!bAllowSuccessiveDelim) && ((pos-prevPos) == 0))
        {
            prevPos = pos + 1; // skip the extra delim.
            continue;
        }

        temp = str.substr((prevPos), (pos-prevPos));
        // cout << "sub string: "<< temp << " @  ("<< prevPos << ", " << pos <<")" << endl;

        prevPos = pos + 1;
        strvec.push_back(bStripWhiteSpace ? StripWhiteSpace(temp) : temp);
    }
    return strvec;
}

// rotate str by rot postions ( positive is forward negative is backward)
inline string InvertString(const string& str)
{
    string invert;
    size_t n = str.length();
    for(size_t i = 0; i < n; i++)
        invert += str[n-1-i];
    return invert;
}


inline std::string rgb_str(const size_t& r, const size_t& g, const size_t& b)
{
    return (utils::NumberToString(r, 2, '0', std::ios_base::hex, std::ios_base::dec)+
            utils::NumberToString(g, 2, '0', std::ios_base::hex, std::ios_base::dec)+
            utils::NumberToString(b, 2, '0', std::ios_base::hex, std::ios_base::dec));
}

/************************************************************************************************************************************************************************/
// General helper functions.
/************************************************************************************************************************************************************************/

template <class TVal, class UVal>
inline bool IsInList(TVal& query, UVal& list, size_t n = 0)
{
    if(n == 0)
        n = list.size();
    for(size_t i = 0; i < n; i++)
    {
        if(query == list[i])
            return true;
    }

    return false;
}

/************************************************************************************************************************************************************************/
//  File and Directory helper functions.
/************************************************************************************************************************************************************************/
#define DOT_CHAR                        '.'
#define DIR_DELIM_CHAR                  '/'
#define WILDCARD_CHAR                   '*'

#define DOT_STR                         "."
#define DIR_DELIM_STR                   "/"
#define WILDCARD_STR                    "*"

inline string   file_part(const string& path) {
    return (path.substr(path.find_last_of(DIR_DELIM_CHAR)+1, path.length()));
}

inline string   file_ext(const string& path) {
    string file = file_part(path);
    return (file.substr(file.find_last_of(DOT_CHAR)+1, file.length()));
}

inline string   file_name(const string& path) {
    string file = file_part(path);
    return (file.substr(0, file.find_last_of(DOT_CHAR)));
}

inline string   path_part(const string& path) {
    return (path.substr(0,path.find_last_of(DIR_DELIM_CHAR)+1));
}

inline string   path_join(const string& path, const string& join) {
    return (path.length() == 0 ? join : (path[path.length()-1] != DIR_DELIM_CHAR && join[0] != DIR_DELIM_CHAR) ? (path+DIR_DELIM_CHAR+join) : ((path[path.length()-1] == DIR_DELIM_CHAR && join[0] != DIR_DELIM_CHAR) || (path[path.length()-1] != DIR_DELIM_CHAR && join[0] == DIR_DELIM_CHAR) ? path+join : path.substr(0, path.length()-1)+join));
}

inline string   search_part(const string& path) {
    return ((path.find(WILDCARD_CHAR) == string::npos) ? path : path.substr(0,path.find(WILDCARD_CHAR)));
}

inline bool     is_directory(const dirent* dir) {
    return (dir->d_type == DT_DIR);
}

//
inline bool file_exists (const std::string& name)
{
    if (FILE *file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    return false;
}

template <class TVal>
inline bool load_txt(vector< vector<TVal> >& data, const std::string& path, string delim = ",", const size_t& reserve = 0, const size_t& skiprows = 0 , const size_t& stoprow = 0, const size_t& axis = 1)
{
    // TODO: check for valid axis.
    std::cout << "load_txt - N path = "<< path << std::endl;
    ifstream txtfile;
    txtfile.open(path.c_str(), ios_base::in);
    if(txtfile)
    {
        string line;
        size_t ln = 0, lnread = 0;
        while (std::getline(txtfile, line, '\n'))
        {
            if(ln % 10000 == 0)
                std::cout << "read "<< ln << " lines!" << std::endl;
            ln++;
            if(ln-1 < skiprows) continue;
            else if (stoprow && ln > stoprow) break;

            lnread++;
            vector<string> split = SplitString(line, delim);

            if(data.size() == 0 && axis == 1)
            {
                data.reserve(split.size());
                for(size_t i = 0; i < split.size(); i++)
                {
                    vector<TVal> temp;
                    data.push_back(temp);
                    if( reserve > 0 )
                        data[i].reserve(reserve);
                }
            }
            if(axis == 1 &&  split.size() != data.size())
            {
                std::cout << "error @ line "<< ln <<" split size = "<< split.size() << ", data size = "<< data.size() << std::endl;
                if( split.size() < data.size())
                {
                    std::cout << "warning!! skipping line "<< ln << std::endl;
                    continue;
                }
                else
                {
                    std::cout << "warning!! trucating the extra columns!" << std::endl;
                }
            }
            vector<TVal> tempAxis0;
            for(size_t i = 0; i < split.size(); i++)
            {
                if(axis == 1 && i < data.size())
                {
                    TVal temp;
                    stringstream ss;
                    ss << split[i];
                    ss >> temp;
                    data[i].push_back( temp );
                }
                else // assume axis = 0;
                {
                    TVal temp;
                    stringstream ss;
                    ss << split[i];
                    ss >> temp;
                    tempAxis0.push_back( temp );
                }
            }
            if(axis != 1) // assume axis = 0;
            {
                data.push_back( tempAxis0 );
            }
        }
        txtfile.close();
    }
    else
    {
        std::cout << "load_txt - X (false)" << std::endl;
        return false;
    }
    std::cout << "load_txt - X (true)" << std::endl;
    return true;
}

template <class TVal>
inline vector< vector<TVal> > load_txt(const std::string& path, string delim = ",", const size_t& reserve = 0, size_t skiprows = 0 , size_t stoprow = 0)
{
    vector< vector<TVal> > data;
    load_txt<TVal>(data, path, delim, reserve, skiprows, stoprow);
    return data;
}

template <class TVal>
inline TVal peek_at_file(const std::string& path, string delim = ",", size_t skiprows = 0)
{
    vector< vector<TVal> > data;
    if(load_txt(data, path, delim, 0, skiprows, skiprows+1)){
        return (data.size() ? (data[0].size() ? data[0][0] : TVal()) : TVal());
    }
    else{
        return TVal();
    }
}
#ifdef c_plus_plus_11
template<class TVal>
inline bool dump_vector(const string& path, const vector<TVal>& v)
{
    ofstream f(path.c_str());
    if(f.is_open())
    {
        for(TVal i : v)
            f << i << endl;
        f.close();
        return true;
    }
    return false;
}

template<class TVal>
inline bool dump_txt(const string& path, const vector< vector<TVal> > & v, string delim = ",")
{
    ofstream f(path.c_str());
    if(f.is_open())
    {
        for(const vector<TVal>& i : v)
        {
            for(const TVal& j : i)
                f << j << delim;
            f << endl;
        }
        f.close();
        return true;
    }
    return false;
}
#endif

// TODO:
//  Clean up this function to make this work more expectedly.
inline vector<string> FindAllPaths(string SearchPath) /*  __uint8_t FileTypeRestriction = 0  not sure hwy this was here. */
{
    dirent* dirNav = NULL;
    DIR* dirc = NULL;
    bool bAdd = false;
    vector<string> paths;
    string parent = path_part(SearchPath);
    string wildCard = file_part(SearchPath);
    vector<string> searchPattern = SplitString(wildCard, WILDCARD_STR);

    dirc = opendir(parent.c_str());

    if(!dirc)
    {
        cerr << " Could not open : " << parent <<endl;
        return paths;
    }

    while((dirNav = readdir(dirc)) && dirc)
    {
        string name(dirNav->d_name);
//        cout <<name<<(is_directory(dirNav) ? " Dir " : " File ")<<endl; // TODO: make this scan recursively.
        bAdd = true;
        for (size_t i = 0; (i < searchPattern.size() && bAdd); i++)
        {
            if(name.find(searchPattern[i]) == string::npos)
                bAdd = false;
        }
        if (bAdd) paths.push_back((parent + name));
    }

    return paths;
}


inline long int GetSizeOfFile(FILE* file)
{
    long int szFile = 0;
    if(file)
    {
        fseek (file , 0 , SEEK_END);
        szFile = ftell (file);
        rewind (file);
    }
    return szFile;
}

inline bool CompareFilesLineByLine(string path1, string path2, bool bPrintOut = false)
{
    bool bcomp = true; char yes;
    ifstream file1; ifstream file2;
    file1.open(path1.c_str()); file2.open(path2.c_str());
    int64 ctLinesMatch = 0;
    int64 ctLinesDiff = 0;

    while(!file1.eof() || (!file2.eof() && bPrintOut))
    {
        char line1[3200];  // 3200 characters, Too small?
        char line2[3200];
        file1.getline(&line1[0], 3200, '\n');
        file2.getline(&line2[0], 3200, '\n');
        string testline1 = line1;
        string testline2 = line2;



        if(strcmp(testline1.c_str(), testline2.c_str()))
        {
            if(bPrintOut)
            {
                cout<<"Lines are different: "<<endl;
                cout<<"file 1 : "<< testline1 << endl <<"file 2 : "<<  testline2 <<endl;
                bcomp = false;

                cout << endl<< " Print next line mismatch (y/n)? " ;
                cin>>yes;

                if(strcmp(&yes, "y"))
                    bPrintOut = true;
                else
                    bPrintOut = false;
            }
            ctLinesDiff++;


        }
        else {
            ctLinesMatch ++;
        }

    }

    cout    << " lines Match : " << ctLinesMatch<<endl
            << " lines Different: " << ctLinesDiff<<endl;


    return bcomp;


}

inline bool CompareFiles(string Path1, string Path2, bool bPrintOut = false )
{
    FILE* pFile1 = NULL; FILE* pFile2 = NULL;
    unsigned char* buf1,*buf2;
    buf1 = NULL; buf2 = NULL;
    size_t szBuffer = 0, ctRead1 = 0, ctRead2 = 0, ctSmall = 0, ctBytesMatch = 0, ctBytesMisMatch = 0, ctTotalRead = 0;
    long int szFile1, szFile2, szSmall;
    szFile1 = 0; szFile2 = 0; szBuffer = 1024*1024; // 1MB buffer
    bool ret_val = true;

    pFile1 = fopen(Path1.c_str(), "r");
    pFile2 = fopen(Path2.c_str(), "r");

    szFile1 = GetSizeOfFile(pFile1);
    szFile2 = GetSizeOfFile(pFile2);
    szSmall = (szFile1 <= szFile2) ? szFile1 : szFile2;

    if(szFile1 && szFile2)
    {
        buf1 = (unsigned char*) malloc(szBuffer);
        buf2 = (unsigned char*) malloc(szBuffer);
        if(!buf1 || !buf2) return false;

        memset(buf1, 0x00, szBuffer);
        memset(buf2, 0x00, szBuffer);

        do
        {

            ctRead1 = 0; ctRead2 = 0;
            ctRead1 = fread(buf1, 1, szBuffer, pFile1);
            ctRead2 = fread(buf2, 1, szBuffer, pFile2);
            ctSmall = (ctRead1 < ctRead2) ? ctRead1 : ctRead2;

            if(memcmp(buf1, buf2, ctSmall) != 0) // there is a difference in the two buffers.
            {
                size_t i = 0; char yes;
                while (buf1[i] == buf2[i]) i++; // find the position where the files are defferent.

                ctBytesMisMatch += ctSmall - i; // Could over estimate the value of mis matched bytes.

                while(bPrintOut && buf1[i] != buf2[i])
                {
                    cout << " Mismatch occured at byte: " << ctTotalRead + i<<endl;
                    cout << " File 1: ";
                    while (buf1[i] != '\n' && i < szBuffer) // print out buffer 1.
                        cout << buf1[i++];
                    cout << buf1[i++];
                    cout << " File 2: ";
                    while (buf2[i] != '\n' && i < szBuffer) // print out buffer 1.
                        cout << buf2[i++];
                    cout << buf2[i++];

                    cout << endl<< " Print next line mismatch (y/n)? " ;
                    cin>>yes;

                    if(strcmp(&yes, "y"))
                        bPrintOut = true;
                    else
                        bPrintOut = false;
                }

            }
            else {
                ctBytesMatch += ctSmall;
            }

            ctTotalRead += ctSmall;

        }while(ctTotalRead < size_t(szSmall) && ctRead1 && ctRead2);


        free(buf1);
        free(buf2);
    }
    else
    {
        cerr<<Path1 <<" size : " <<szFile1 <<endl<<Path2<< " size : " << szFile2<<endl;
    }

    cout    << "results" << endl<< " File 1 size : " << szFile1<<endl
            << " File 2 size: " << szFile2 << endl
            << " Bytes Match : " << ctBytesMatch<<endl
            << " Bytes Different: " << ctBytesMisMatch<<endl;

    fclose(pFile1);
    fclose(pFile2);

    return ret_val;
}


/************************************************************************************************************************************************************************/
// Data analysis helper classes and functions.
/************************************************************************************************************************************************************************/

// Below are classes and functions that rely on c++11 support.
#ifdef c_plus_plus_11


class CHistogram
{
    public:
        template<class TVal>
        CHistogram(size_t nbins, vector<TVal>& data, const string& path = "", const float* const min = NULL, const float* const max = NULL) :  m_Path(path), m_nBins(nbins), m_pHist(NULL), m_pBinsEdges(NULL)
        {
            m_pHist = std::shared_ptr<size_t>(new size_t[m_nBins]);
            m_pBinsEdges = std::shared_ptr<float>(new float[m_nBins+1]);
            InitializeArrays();
            MakeHistogram(data, min, max);

        }
        ~CHistogram()
        {
            m_pHist.reset();
            m_pBinsEdges.reset();

            // cout << "hist:" << hex << m_pHist.get() << " bins: "<< m_pBinsEdges.get() << dec << endl;
        }
        void InitializeArrays()
        {
            for(size_t b = 0; b < m_nBins+1; b++)
                m_pBinsEdges.get()[b] = 0;

            for(size_t d = 0; d < m_nBins; d++)
                m_pHist.get()[d] = 0;

        }

        bool Save()
        {
            return SaveToFile(m_Path);
        }
        bool SaveToFile(const string& path)
        {
            FILE* pFile = NULL;
            pFile = fopen(path.c_str(), "w");
            bool bSuccess = (pFile);
            if(pFile)
            {
                PrintHistogram(pFile);
                fclose(pFile);
            }
            return bSuccess;
        }

        void PrintHistogram(FILE* f = stdout)
        {
            for(size_t i = 0; i < m_nBins; i++)
            {
                float width = m_pBinsEdges.get()[i+1]-m_pBinsEdges.get()[i];
                float center = m_pBinsEdges.get()[i]+(width/2.0f);
                fprintf(f, "%f %lu %f \n", center, m_pHist.get()[i], width);
            }
            fflush(f);
        }

    protected:
        template<class TVal>
        void MakeHistogram(vector<TVal>& data, const float* const min, const float* const max)
        {
            float mindata = (!min) ? int(mininum(data)) : *min;
            float maxdata = (!max) ? int(maximum(data))+1 : *max;
            assert(mindata < maxdata);
            float delta = (maxdata - mindata)/float(m_nBins);
            cout << "min = "<< mindata << " max = " << maxdata << " nbins = "<< m_nBins << endl;

            // Set the bin edges.
            for(size_t b = 0; b < m_nBins+1; b++)
                m_pBinsEdges.get()[b] = delta*b + mindata;

            // Make the histogram.
            size_t bndx = 0;
            for(size_t d = 0; d < data.size(); d++)
            {
                bndx = int((data[d]-mindata)/delta);
                m_pHist.get()[bndx]++;

                assert(BOUNDED(data[d], m_pBinsEdges.get()[bndx], m_pBinsEdges.get()[bndx+1]));
                if(!BOUNDED(data[d], m_pBinsEdges.get()[bndx], m_pBinsEdges.get()[bndx+1]))
                {
                    cout << "Error! data out of histogram bin ranges! "<< endl;
                }
            }

        }


    private:
        string m_Path;
        size_t m_nBins;
        std::shared_ptr<size_t>  m_pHist;
        std::shared_ptr<float>   m_pBinsEdges;
};





// My tuple utility functions.

// base case
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I == sizeof...(Ts), void >::type put (
                                                                        std::tuple< Ts... >& ,
                                                                        TVal& ,
                                                                        const size_t&)
{
    return;
}
// induction
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I < sizeof...(Ts), void >::type put (
                                                                        std::tuple< Ts... >& t,
                                                                        TVal& val,
                                                                        const size_t& Index)
{
    typedef typename std::tuple_element< I, std::tuple< Ts... > >::type elem_type;
    if(I == Index){
        std::get<I>(t) = elem_type(val);
        }
    else
        put< (I+1), TVal, Ts... > (t, val, Index);
}

// base case
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I == sizeof... (Ts), bool >::type is_type (
                                                                            std::tuple< Ts... >& ,
                                                                            const size_t&)
{
    return false; // Index out of range so return false.
}
// induction
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I < sizeof...(Ts), bool >::type is_type (
                                                                            std::tuple< Ts... >& t,
                                                                            const size_t& Index)
{
    if(I == Index)
        return typeid(std::get<I>(t)) == typeid(TVal);
    else
        return is_type<I+1, TVal, Ts...>(t, Index);
}

// base case
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I == sizeof... (Ts), TVal* >::type pull (
                                                                            std::tuple< Ts... >& ,
                                                                            const size_t& )
{
    return NULL; // Index out of range so return NULL.
}
// induction
template< std::size_t I = 0,  class TVal, typename... Ts >
inline typename std::enable_if< I < sizeof...(Ts), TVal* >::type pull (
                                                                            std::tuple< Ts... >& t,
                                                                            const size_t& Index)
{
    if(I == Index)
        return (TVal*)&std::get<I>(t);
    else
        return pull<I+1, TVal, Ts...>(t, Index);
}

template< typename FromVal, typename ToVal>
inline typename std::enable_if< !is_same<FromVal, ToVal>::value, void >::type assign_from (
                                                                                                FromVal&,
                                                                                                ToVal&)
{
    return;
}

template< typename FromVal, typename ToVal>
inline typename std::enable_if< is_same<FromVal, ToVal>::value, void >::type assign_from (
                                                                                                FromVal& from,
                                                                                                ToVal&  to)
{
    to = from;
}

// base case
template< std::size_t I = 0, class TVal, typename... Ts >
inline typename std::enable_if< I == sizeof... (Ts), void >::type pull2 (
                                                                            std::tuple< Ts... >&,
                                                                            TVal&,
                                                                            const size_t&)
{
    return; // Index out of range so return NULL.
}
// induction
template< std::size_t I = 0,  class TVal, typename... Ts >
inline typename std::enable_if< I < sizeof...(Ts), void >::type pull2 (
                                                                            std::tuple< Ts... >& t,
                                                                            TVal& val,
                                                                            const size_t& Index)
{
    if(I == Index)
        val = std::get<I>(t);
    else
        pull2<I+1, TVal, Ts...>(t, val, Index);
}


// base case
template< std::size_t I = 0>
inline typename std::enable_if< I == 5, void>::type count_to_five_or_less(const size_t&)
{
    return;
}
// induction
template< std::size_t I = 0>
inline typename std::enable_if< I < 5, void>::type count_to_five_or_less(const size_t& Index)
{
    if(I < Index)
        cout << I << endl;
    count_to_five_or_less<I+1>(Index);
}

#endif



}
















#endif
