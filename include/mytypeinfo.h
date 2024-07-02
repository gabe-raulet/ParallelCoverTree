#ifndef MY_TYPE_INFO_H_
#define MY_TYPE_INFO_H_

#include <type_traits>
#include <vector>
#include <tuple>
#include <array>

using namespace std;

template <class T>
struct is_tuple : false_type {};

template <class ...Ts>
struct is_tuple<tuple<Ts...>> : true_type {};

template <class T>
inline constexpr bool is_tuple_v = is_tuple<T>::value;

template <class T>
struct is_std_array : false_type {};

template <class T, size_t N>
struct is_std_array<array<T, N>> : true_type {};

template <class T>
inline constexpr bool is_array_type = is_std_array<T>::value || is_array<T>::value;

template <class T>
concept array_type = is_array_type<T>;

template <array_type T>
struct array_info;

template <class T, size_t N>
struct array_info<T[N]>
{
    using value_type = T;
    static constexpr size_t size = N;
};

template <class T, size_t N>
struct array_info<array<T, N>>
{
    using value_type = T;
    static constexpr size_t size = N;
};

template <class T>
inline constexpr size_t array_size = array_info<T>::size;

template <class T>
concept string_type = same_as<T, ::string> || same_as<T, vector<char>>;

#endif
