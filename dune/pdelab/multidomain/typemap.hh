#ifndef DUNE_MULTIDOMAIN_TYPEMAP_HH
#define DUNE_MULTIDOMAIN_TYPEMAP_HH

#include <type_traits>

namespace Dune {

namespace PDELab {

namespace MultiDomain {

// An associative map for types

namespace {

// The trick to make the map work without having to do all the type
// testing by hand is by building an inheritance chain of map_entry
// structs from the key-value pairs. Each map_entry struct overloads
// the method lookup() for its key-value combination. The parameter
// passed in is not the type itself, but a trivial struct templated
// on the type(this way, we do not require the key to be constructible).
// We return the value type.
// When looking up a key, we simply exploit the function dispatch
// mechanism of C++ in combination with the C++0x decltype operator.

template<typename T>
struct map_key {};

struct no_key {};

}

/**
 * Retrieves the mapped value for the given key from the map.
 *
 * Raises a static assertion if the key is not contained in the map.
 */
template<typename key, typename map>
struct get_map_entry
{
  typedef decltype(map().lookup(map_key<key>())) type;

  static_assert((!std::is_same<type,no_key>::value),"key type not contained in type map");
};

/**
 * Tests whether the given key is contained in the map.
 */
template<typename key, typename map>
struct map_contains
{
  static const bool value = !std::is_same<decltype(map().lookup(map_key<key>())),no_key>::value;
};

namespace {

struct empty_type_map
{

  static const bool has_duplicate_entries = false;

  template<typename T>
  no_key lookup(map_key<T>);
};

template<typename key, typename value, typename tail = empty_type_map>
struct map_entry : public tail
{
  using tail::lookup;

  static const bool has_duplicate_entries = tail::has_duplicate_entries || map_contains<key,tail>::value;

  value lookup(map_key<key>);

};

}

/**
 * Adds a key-value pair to the given map.
 */
template<typename key, typename value, typename map>
struct add_map_entry
{
  typedef map_entry<key,value,map> type;
};

namespace {

template<typename... entries>
struct embedded_key_map_helper;

template<typename entry, typename... entries>
struct embedded_key_map_helper<entry,entries...>
{
  typedef typename add_map_entry<typename entry::key,entry,typename embedded_key_map_helper<entries...>::type>::type type;
};

template<>
struct embedded_key_map_helper<>
{
  typedef empty_type_map type;
};

}

/**
 * A type map for entries which contain the key type as a typedef value_type::key.
 */
template<typename... entries>
struct embedded_key_map : public embedded_key_map_helper<entries...>::type
{
};

} // namespace MultiDomain
} // namespace PDELab
} // namespace Dune

#endif // DUNE_MULTIDOMAIN_TYPEMAP_HH
