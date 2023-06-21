/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <deque>

using namespace std;

template <int Dim>
bool smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim)
{
  if (first[curDim] == second[curDim]) {
    return first < second;
  } else {
    return (first[curDim] < second[curDim]);
  }
}

template <int Dim>
bool shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential)
{
    int dist_one = 0;
    int dist_two = 0;
    for (int i = 0; i < Dim; i++) {
      dist_one += (target[i] - currentBest[i]) * (target[i] - currentBest[i]);
    }
    for (int i = 0; i < Dim; i++) {
      dist_two += (target[i] - potential[i]) * (target[i] - potential[i]);
    }
    return dist_one > dist_two;
}

template<int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::buildKDTree(vector<Point<Dim>>& points, int left, int right, int currDim) {
    if (left > right) {
        return nullptr;
    }
    int mid = left + (right - left) / 2;
    select(points.begin() + left, points.begin() + right + 1, points.begin() + mid, ComparePoints(currDim));
    KDTreeNode* node = new KDTreeNode(points[mid]);
    node->left = buildKDTree(points, left, mid - 1, (currDim + 1) % Dim);
    node->right = buildKDTree(points, mid + 1, right, (currDim + 1) % Dim);
    return node;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
  if (newPoints.empty()) {
        return;
    }
    vector<Point<Dim>> points = newPoints;
    root = buildKDTree(points, 0, points.size() - 1, 0);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
   if (this != &rhs) {
        // Create a new root node and recursively copy the tree
        root = copyNode(rhs.root);
    }
    return *this;
}

template<int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::copyNode(const KDTreeNode* node) {
    if (node == nullptr) {
        return nullptr;
    } else {
        KDTreeNode* newNode = new KDTreeNode(node->point, node->left, node->right);
        newNode->left = copyNode(node->left);
        newNode->right = copyNode(node->right);
        return newNode;
    }
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    Point<Dim> nearest;
    double distance = std::numeric_limits<double>::infinity();
    findNearestNeighborHelper(root, query, 0, nearest, distance);
    return nearest;
}

template<int Dim>
void KDTree<Dim>::findNearestNeighborHelper(const KDTreeNode* node, const Point<Dim>& query,
                                             int level, Point<Dim>& nearest, double& distance) const {
    if (node == nullptr) {
        return;
    }
    double currDist = distanceSquared(node->point, query);
    if (currDist < distance) {
        nearest = node->point;
        distance = currDist;
    } else if (currDist == distance) {
        if (node->point < nearest) {
            nearest = node->point;
        }
    }
    int nextDim = (level + 1) % Dim;
    bool leftFirst = (query[level] < node->point[level]);
    if (leftFirst) {
        findNearestNeighborHelper(node->left, query, nextDim, nearest, distance);
        if ((node->point[level] - query[level]) * (node->point[level] - query[level]) < distance) {
            findNearestNeighborHelper(node->right, query, nextDim, nearest, distance);
        } else if ((node->point[level] - query[level]) * (node->point[level] - query[level]) == distance) {
            if (node->point < nearest) {
                findNearestNeighborHelper(node->right, query, nextDim, nearest, distance);
            } else {
                findNearestNeighborHelper(node->left, query, nextDim, nearest, distance);
            }
        }
    } else {
        findNearestNeighborHelper(node->right, query, nextDim, nearest, distance);
        if ((node->point[level] - query[level]) * (node->point[level] - query[level]) < distance) {
            findNearestNeighborHelper(node->left, query, nextDim, nearest, distance);
        } else if ((node->point[level] - query[level]) * (node->point[level] - query[level]) == distance) {
            if (node->point < nearest) {
                findNearestNeighborHelper(node->left, query, nextDim, nearest, distance);
            } else {
                findNearestNeighborHelper(node->right, query, nextDim, nearest, distance);
            }
        }
    }
}

template <int Dim>
int KDTree<Dim>::distanceSquared(const Point<Dim>& first,
                                 const Point<Dim>& second) const {
    int distanceSquared = 0;
    for (int i = 0; i < Dim; i++) {
        distanceSquared += (first[i] - second[i]) * (first[i] - second[i]);
    }
    return distanceSquared;
}

template <typename RandIter, typename Comparator>
void select(RandIter start, RandIter end, RandIter k, Comparator cmp)
{
  while (true) {
        if (start == end) {
            return;
        }

        // Choose a pivot element at random from the range
        RandIter pivot = start + std::rand() % (end - start);

        // Swap the pivot element with the last element in the range
        std::iter_swap(pivot, end - 1);

        // Partition the elements around the pivot element
        RandIter mid = start;
        for (RandIter i = start; i != end - 1; ++i) {
            if (cmp(*i, *(end - 1))) {
                std::iter_swap(i, mid);
                ++mid;
            }
        }

        // Swap the pivot element back to its final position
        std::iter_swap(mid, end - 1);

        // Check the position of the pivot element
        if (k == mid) {
            return;
        } else if (k < mid) {
            end = mid;
        } else {
            start = mid + 1;
        }
    }
}

