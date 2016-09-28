#ifndef MEMERNA_SPLAYMAP_H
#define MEMERNA_SPLAYMAP_H

#include "common.h"

namespace memerna {

template<typename Key, typename Value>
class SplayMap {
public:
  SplayMap() : ns(2), vals(2), root(NONE), size(0) {}
  // Returns false if already in the tree.
  bool Insert(const Key& key, const Value& value) {
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
  bool Find(const Key& key) {
    if (root == NONE) return false;  // TODO remove this check?
    // TODO compare performance replacing Zig-Zag with just Zig
    int lnext = TMP, rnext = TMP;
    bool found = false;
    while (1) {
      auto& n = ns[root];
      /*    n
       *   / \
       *  l   r
       * Call the new root we are going to 'x'.
       */
      if (key < n.k) {
        // Case: we are going left
        int l = n.l;
        if (l == NONE) break;
        if (key < ns[l].k) {
          // Zig-Zig
          RotateRight(root);  // Splay tree condition. l is now the parent of x
          if (ns[l].l == NONE) {
            // Node isn't in the tree (too small), so set root to the nearest one.
            root = l;
            break;
          }
          root = ns[l].l;  // Update root.
          ns[l].l = NONE;  // Split left child
          ns[rnext].l = l;  // Add to the R tree.
          rnext = l;  // Update next lowest point for R tree.

        } else if (ns[l].k < key) {
          // Zig-Zag
          // Zig:
          n.l = NONE;  // Split left child
          ns[rnext].l = root;  // Add to R tree.
          rnext = root;  // Update next lowest point for R tree.
          if (ns[l].r == NONE) {
            // Node isn't in the tree, so set root to nearest.
            root = l;
            break;
          }
          root = ns[l].r;  // Update root.
          // Zag:
          ns[l].r = NONE;  // Split right child next level down.
          ns[lnext].r = l;  // Add to L tree.
          lnext = l; // Update next lowest point.
        } else {
          // Found (zig).
          n.l = NONE;  // Split left child
          ns[rnext].l = root;  // Add to R tree.
          rnext = root;  // TODO unnecessary
          root = l;  // Update root.
          found = true;
          break;
        }
      } else if (n.k < key) {
        // Case: we are going right
        int r = n.r;
        if (r == NONE) break;
        if (ns[r].k < key) {
          // Zig-Zig
          RotateLeft(root);  // Splay tree condition. r is now the parent of x
          if (ns[r].r == NONE) {
            // Node isn't in the tree (too big), so set root to the nearest one.
            root = r;
            break;
          }
          root = ns[r].r;  // Update root.
          ns[r].r = NONE;  // Split right child
          ns[lnext].r = r;  // Add to the L tree.
          lnext = r;  // Update next lowest point for L tree.
        } else if (key < ns[r].k) {
          // Zig-Zag
          // Zig:
          n.r = NONE;  // Split right child
          ns[lnext].r = root;  // Add to L tree.
          lnext = root;  // Update next lowest point for L tree.
          if (ns[r].l == NONE) {
            // Node isn't in the tree, so set root to nearest.
            root = r;
            break;
          }
          root = ns[r].l;  // Update root.
          // Zag:
          ns[r].l = NONE;  // Split left child next level down.
          ns[rnext].l = r;  // Add to R tree.
          rnext = r; // Update next lowest point.
        } else {
          // Found (zig).
          n.r = NONE;  // Split right child
          ns[lnext].r = root;  // Add to L tree.
          lnext = root;  // TODO unnecessary
          root = r;  // Update root.
          found = true;
          break;
        }
      } else {
        // Found.
        found = true;
        break;
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
  bool Delete(const Key& key) {
    if (Find(key)) {
      // Root now contains the key.
      if (ns[root].r == NONE) {
        // Either, the right subtree is empty, in which case set the root to the left subtree:
        root = ns[root].l;
      } else {
        // Or the right subtree is not empty, in which case:
        // Move the next lowest key up to the top of the right subtree
        // with another find. Since it is the next lowest, the left child of the right subtree
        // will be NONE, so we can attach the left subtree there.
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

  Value& operator[](const Key& key) {
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

  void RotateLeft(int p) {
    int x = ns[p].r;
    assert(x != NONE);
    ns[p].r = ns[x].l;
    ns[x].l = p;
  }

  void RotateRight(int p) {
    int x = ns[p].l;
    assert(x != NONE);
    ns[p].l = ns[x].r;
    ns[x].r = p;
  }

  std::vector<std::string> DescribeInternal(int node) {
    if (node == NONE) return {""};
    const auto& n = ns[node];
    std::vector<std::string> desc;
    desc.push_back(std::to_string(n.k));
    for (const auto& s : DescribeInternal(n.l))
      desc.push_back("| " + s);
    int idx = int(desc.size());
    for (const auto& s : DescribeInternal(n.r))
      desc.push_back(s);
    desc[1][1] = '-';
    desc[idx] = "|-" + desc[idx];
    return desc;
  }
};

}

#endif  // MEMERNA_SPLAYMAP_H
