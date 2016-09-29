#ifndef MEMERNA_SPLAYMAP_H
#define MEMERNA_SPLAYMAP_H

#include "common.h"

namespace memerna {

template<typename Key, typename Value>
class SplayMap {
public:
  SplayMap() : ns(2), vals(2), root(NONE), size(0) {}
  // Returns false if already in the tree.
  bool Insert(Key key, const Value& value) {
    if (Find(key)) return false;
    ++size;
    vals.push_back(value);

    int oldroot = root;
    root = int(ns.size());
    // Case where |oldroot| is NONE will be handled magically, since none.{l, r} == none.
    if (key < ns[oldroot].k) {
      ns.push_back({key, ns[oldroot].l, oldroot});
      ns[oldroot].l = NONE;
    } else {
      ns.push_back({key, oldroot, ns[oldroot].r});
      ns[oldroot].r = NONE;
    }
    assert(ns[NONE].l == NONE && ns[NONE].r == NONE);  // Keep this invariant.
    return true;
  }

  // Returns true if found.
  bool Find(Key key) {
    if (root == NONE) return false;
    int lnext = TMP, rnext = TMP;
    ns[TMP].l = ns[TMP].r = NONE;
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
          root = l;
          l = ns[root].l;
          if (l == NONE) break;
          // Split left child
          ns[rnext].l = root;
          rnext = root;
          root = l;
        } else if (ns[l].k < key) {
          // Zig-Zag
          // Zig - Split left child
          ns[rnext].l = root;
          rnext = root;
          root = l;
          if (ns[root].r == NONE) break;
          // Zag - Split right child
          root = ns[root].r;
          ns[lnext].r = l;
          lnext = l;
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
          root = r;
          r = ns[root].r;
          if (r == NONE) break;
          // Split right child
          ns[lnext].r = root;
          lnext = root;
          root = r;
        } else if (key < ns[r].k) {
          // Zig-Zag
          // Zig - Split right child
          ns[lnext].r = root;
          lnext = root;
          root = r;
          if (ns[root].l == NONE) break;
          // Zag - Split left child.
          root = ns[root].l;
          ns[rnext].l = r;
          rnext = r;
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
    return vals[root];
  }

  Value& operator[](Key key) {
    if (Find(key)) return Get();
    Insert(key, Value());
    return Get();
  }

  std::size_t Size() { return size; }

  std::string Describe() {
    std::string ans = sfmt(
        "Tree with %zu nodes. Backing node size: %zu, Backing vals size: %zu, root at index %d\n",
        Size(), ns.size(), vals.size(), root);
    for (const auto& s : DescribeInternal(root))
      ans += s + "\n";
    return ans;
  }
private:
  constexpr static int NONE = 0, TMP = 1;
  struct node_t {
    node_t() : k(), l(NONE), r(NONE) {}
    node_t(const Key& k_, int l_, int r_) : k(k_), l(l_), r(r_) {}
    Key k;
    int l, r;
  };
  std::vector<node_t> ns;  // TODO: add none at start, eliminate comparisons
  std::vector<Value> vals;
  int root;
  std::size_t size;

  std::vector<std::string> DescribeInternal(int node) {
    if (node == NONE) return {""};
    const auto& n = ns[node];
    std::vector<std::string> desc;
    desc.push_back(std::to_string(n.k));
    for (const auto& s : DescribeInternal(n.l))
      desc.push_back("| " + s);
    int idx = int(desc.size());
    for (const auto& s : DescribeInternal(n.r))
      desc.push_back("  " + s);
    desc[1][1] = '_';
    desc[idx][0] = '|';
    desc[idx][1] = '_';
    return desc;
  }
};

}

#endif  // MEMERNA_SPLAYMAP_H
