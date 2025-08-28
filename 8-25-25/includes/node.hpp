#ifndef NODE_H
#define NODE_H

#include <memory>

template <typename T>
struct Node {
    std::unique_ptr<T> data;
    Node* next;
    Node* previous;
    Node() : data(std::make_unique<T>()), next(nullptr), previous(nullptr) {}
    Node(const T& v) : data(std::make_unique<T>(v)), next(nullptr), previous(nullptr) {}
};

#endif // NODE_H