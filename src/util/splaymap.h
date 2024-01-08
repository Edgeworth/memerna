// Copyright 2016 Eliot Courtney.
#ifndef UTIL_SPLAYMAP_H_
#define UTIL_SPLAYMAP_H_

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace mrna {

template <typename Key, typename Value>
class SplayMap {
 public:
  SplayMap() : ns_(2), root(NONE) {}
  // Returns false if already in the tree.
  template <typename ValueRef>
  bool Insert(Key key, ValueRef&& value) {
    if (Find(key)) return false;
    ++size_;

    int oldroot = root;
    root = static_cast<int>(ns_.size());
    // Case where `oldroot` is NONE will be handled magically, since none.{l, r} == none.
    if (key < ns_[oldroot].k) {
      ns_.push_back({key, ns_[oldroot].l, oldroot, std::forward<ValueRef>(value)});
      ns_[oldroot].l = NONE;
    } else {
      ns_.push_back({key, oldroot, ns_[oldroot].r, std::forward<ValueRef>(value)});
      ns_[oldroot].r = NONE;
    }
    assert(ns_[NONE].l == NONE && ns_[NONE].r == NONE);  // Keep this invariant.
    return true;
  }

  // Returns true if found.
  bool Find(Key key) {
    if (root == NONE) return false;
    int lnext = TMP;
    int rnext = TMP;
    bool found = false;
    while (!found) {
      if (key < ns_[root].k) {
        // Case: we are going left
        int l = ns_[root].l;
        if (l == NONE) break;
        if (key < ns_[l].k) {
          // Zig-Zig - Rotate right
          ns_[root].l = ns_[l].r;
          ns_[l].r = root;
          if (ns_[l].l == NONE) {  // No left child
            root = l;
            break;
          }  // Split left child
          root = ns_[l].l;
          ns_[rnext].l = l;
          rnext = l;

        } else if (ns_[l].k < key) {
          // Zig - Split left child
          ns_[rnext].l = root;
          rnext = root;
          if (ns_[l].r == NONE) {  // No right child.
            root = l;
            break;
          }  // Zag - Split right child
          root = ns_[l].r;
          ns_[lnext].r = l;
          lnext = l;

        } else {
          // Found (zig) - Split left child
          ns_[rnext].l = root;
          rnext = root;
          root = l;
          found = true;
        }
      } else if (ns_[root].k < key) {
        // Case: we are going right
        int r = ns_[root].r;
        if (r == NONE) break;
        if (ns_[r].k < key) {
          // Zig-Zig - Rotate left
          ns_[root].r = ns_[r].l;
          ns_[r].l = root;
          if (ns_[r].r == NONE) {  // No right child.
            root = r;
            break;
          }  // Split right child
          root = ns_[r].r;
          ns_[lnext].r = r;
          lnext = r;

        } else if (key < ns_[r].k) {
          // Zig - Split right child
          ns_[lnext].r = root;
          lnext = root;
          if (ns_[r].l == NONE) {  // No left child
            root = r;
            break;
          }  // Zag - Split left child.
          root = ns_[r].l;
          ns_[rnext].l = r;
          rnext = r;

        } else {
          // Found (zig) - Split right child
          ns_[lnext].r = root;
          lnext = root;
          root = r;
          found = true;
        }
      } else {
        // Found.
        found = true;
      }
    }
    // Reassemble the tree.
    ns_[lnext].r = ns_[root].l;
    ns_[rnext].l = ns_[root].r;
    ns_[root].l = ns_[TMP].r;
    ns_[root].r = ns_[TMP].l;
    assert(ns_[NONE].l == NONE && ns_[NONE].r == NONE);  // Keep this invariant.
    return found;
  }

  // Returns true if successfully deleted. This only pretend deletes the data.
  bool Delete(Key key) {
    if (Find(key)) {
      // Root now contains the key.
      if (ns_[root].r == NONE) {
        // Either, the right subtree is empty, in which case set the root to the left subtree:
        root = ns_[root].l;
      } else {
        // Or the right subtree is not empty, in which case:
        // Move the next lowest key up to the top of the right subtree with another find.
        // Since it is the next lowest, the left child of the right subtree will be NONE,
        // so we can attach the left subtree there.
        int oldroot = root;
        root = ns_[root].r;
        Find(key);
        assert(ns_[root].l == NONE);
        ns_[root].l = ns_[oldroot].l;
      }
      --size_;
      return true;
    }
    return false;
  }

  const Value& Get() {
    assert(size_ > 0);
    return ns_[root].v;
  }

  Value& operator[](Key key) {
    if (Find(key)) return Get();
    Insert(key, Value());
    return Get();
  }

  std::size_t Size() { return size_; }

  void Reserve(std::size_t s) { ns_.reserve(s); }

  void Clear() {
    ns_.resize(2);
    assert(ns_[NONE].l == NONE && ns_[NONE].r == NONE);
    root = NONE;
    size_ = 0;
  }

  // Testing / visualisation methods.
  std::string Describe() {
    std::string ans = fmt::format(
        "Tree with {} nodes. Backing node size: {}, root at index {}\n", Size(), ns_.size(), root);
    for (const auto& s : DescribeInternal(root)) ans += s + "\n";
    return ans;
  }

  std::vector<Key> Keys() { return KeysInternal(root); }

 private:
  inline static constexpr int NONE = 0, TMP = 1;

  struct Node {
    Key k;
    int l, r;
    Value v;
  };

  std::vector<Node> ns_;
  int root;
  std::size_t size_{0};

  std::vector<std::string> DescribeInternal(int node) {
    if (node == NONE) return {""};
    const auto& n = ns_[node];
    std::vector<std::string> desc;
    desc.push_back(std::to_string(n.k));
    for (const auto& s : DescribeInternal(n.l)) desc.push_back("| " + s);
    const int idx = static_cast<int>(desc.size());
    for (const auto& s : DescribeInternal(n.r)) desc.push_back("  " + s);
    desc[1][1] = '_';
    desc[idx][0] = '|';
    desc[idx][1] = '_';
    return desc;
  }

  std::vector<Key> KeysInternal(int node) {
    if (node == NONE) return {};
    auto a = KeysInternal(ns_[node].l);
    a.push_back(ns_[node].k);
    auto b = KeysInternal(ns_[node].r);
    std::vector<Key> c(a.size() + b.size());
    std::merge(a.begin(), a.end(), b.begin(), b.end(), c.begin());
    return c;
  }
};

struct Nothing {};

template <typename Key>
using SplaySet = SplayMap<Key, Nothing>;

}  // namespace mrna

#endif  // UTIL_SPLAYMAP_H_
