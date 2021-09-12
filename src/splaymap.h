// Copyright 2016 E.
#ifndef SPLAYMAP_H_
#define SPLAYMAP_H_

#include <algorithm>

#include "common.h"

namespace memerna {

template <typename Key, typename Value>
class SplayMap {
 public:
  SplayMap() : ns(2), root(NONE), size(0) {}
  // Returns false if already in the tree.
  template <typename ValueRef>
  bool Insert(Key key, ValueRef&& value) {
    if (Find(key)) return false;
    ++size;

    int oldroot = root;
    root = int(ns.size());
    // Case where |oldroot| is NONE will be handled magically, since none.{l, r} == none.
    if (key < ns[oldroot].k) {
      ns.push_back({key, ns[oldroot].l, oldroot, std::forward<ValueRef>(value)});
      ns[oldroot].l = NONE;
    } else {
      ns.push_back({key, oldroot, ns[oldroot].r, std::forward<ValueRef>(value)});
      ns[oldroot].r = NONE;
    }
    assert(ns[NONE].l == NONE && ns[NONE].r == NONE);  // Keep this invariant.
    return true;
  }

  // Returns true if found.
  bool Find(Key key) {
    if (root == NONE) return false;
    int lnext = TMP, rnext = TMP;
    bool found = false;
    while (!found) {
      if (key < ns[root].k) {
        // Case: we are going left
        int l = ns[root].l;
        if (l == NONE) break;
        if (key < ns[l].k) {
          // Zig-Zig - Rotate right
          ns[root].l = ns[l].r;
          ns[l].r = root;
          if (ns[l].l == NONE) {  // No left child
            root = l;
            break;
          } else {  // Split left child
            root = ns[l].l;
            ns[rnext].l = l;
            rnext = l;
          }
        } else if (ns[l].k < key) {
          // Zig - Split left child
          ns[rnext].l = root;
          rnext = root;
          if (ns[l].r == NONE) {  // No right child.
            root = l;
            break;
          } else {  // Zag - Split right child
            root = ns[l].r;
            ns[lnext].r = l;
            lnext = l;
          }
        } else {
          // Found (zig) - Split left child
          ns[rnext].l = root;
          rnext = root;
          root = l;
          found = true;
        }
      } else if (ns[root].k < key) {
        // Case: we are going right
        int r = ns[root].r;
        if (r == NONE) break;
        if (ns[r].k < key) {
          // Zig-Zig - Rotate left
          ns[root].r = ns[r].l;
          ns[r].l = root;
          if (ns[r].r == NONE) {  // No right child.
            root = r;
            break;
          } else {  // Split right child
            root = ns[r].r;
            ns[lnext].r = r;
            lnext = r;
          }
        } else if (key < ns[r].k) {
          // Zig - Split right child
          ns[lnext].r = root;
          lnext = root;
          if (ns[r].l == NONE) {  // No left child
            root = r;
            break;
          } else {  // Zag - Split left child.
            root = ns[r].l;
            ns[rnext].l = r;
            rnext = r;
          }
        } else {
          // Found (zig) - Split right child
          ns[lnext].r = root;
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
    ns[lnext].r = ns[root].l;
    ns[rnext].l = ns[root].r;
    ns[root].l = ns[TMP].r;
    ns[root].r = ns[TMP].l;
    assert(ns[NONE].l == NONE && ns[NONE].r == NONE);  // Keep this invariant.
    return found;
  }

  // Returns true if successfully deleted. This only pretend deletes the data.
  bool Delete(Key key) {
    if (Find(key)) {
      // Root now contains the key.
      if (ns[root].r == NONE) {
        // Either, the right subtree is empty, in which case set the root to the left subtree:
        root = ns[root].l;
      } else {
        // Or the right subtree is not empty, in which case:
        // Move the next lowest key up to the top of the right subtree with another find.
        // Since it is the next lowest, the left child of the right subtree will be NONE,
        // so we can attach the left subtree there.
        int oldroot = root;
        root = ns[root].r;
        Find(key);
        assert(ns[root].l == NONE);
        ns[root].l = ns[oldroot].l;
      }
      --size;
      return true;
    }
    return false;
  }

  const Value& Get() {
    assert(Size() > 0);
    return ns[root].v;
  }

  Value& operator[](Key key) {
    if (Find(key)) return Get();
    Insert(key, Value());
    return Get();
  }

  std::size_t Size() { return size; }

  void Reserve(std::size_t s) { ns.reserve(s); }

  void Clear() {
    ns.resize(2);
    assert(ns[NONE].l == NONE && ns[NONE].r == NONE);
    root = NONE;
    size = 0;
  }

  // Testing / visualisation methods.
  std::string Describe() {
    std::string ans = sfmt(
        "Tree with %zu nodes. Backing node size: %zu, root at index %d\n", Size(), ns.size(), root);
    for (const auto& s : DescribeInternal(root)) ans += s + "\n";
    return ans;
  }

  std::vector<Key> Keys() { return KeysInternal(root); }

 private:
  constexpr static int NONE = 0, TMP = 1;

  struct node_t {
    Key k;
    int l, r;
    Value v;
  };

  std::vector<node_t> ns;
  int root;
  std::size_t size;

  std::vector<std::string> DescribeInternal(int node) {
    if (node == NONE) return {""};
    const auto& n = ns[node];
    std::vector<std::string> desc;
    desc.push_back(std::to_string(n.k));
    for (const auto& s : DescribeInternal(n.l)) desc.push_back("| " + s);
    int idx = int(desc.size());
    for (const auto& s : DescribeInternal(n.r)) desc.push_back("  " + s);
    desc[1][1] = '_';
    desc[idx][0] = '|';
    desc[idx][1] = '_';
    return desc;
  }

  std::vector<Key> KeysInternal(int node) {
    if (node == NONE) return {};
    auto a = KeysInternal(ns[node].l);
    a.push_back(ns[node].k);
    auto b = KeysInternal(ns[node].r);
    std::vector<Key> c(a.size() + b.size());
    std::merge(a.begin(), a.end(), b.begin(), b.end(), c.begin());
    return c;
  }
};

}  // namespace memerna

#endif  // SPLAYMAP_H_
