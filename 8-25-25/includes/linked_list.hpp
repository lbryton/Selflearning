// linked_list.hpp
#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include "node.hpp"
#include <iostream>

template <typename T>
class LinkedList {

private:
  Node<T>* head;
  Node<T>* tail;

  void insert_between(Node<T>* prevNode, Node<T>* currNode, Node<T>* nextNode) {
    currNode->next = nextNode;
    currNode->previous = prevNode;
    prevNode->next = currNode;
    nextNode->previous = currNode;
  }

  void delete_between(Node<T>* prevNode, Node<T>* currNode, Node<T>* nextNode) {
    prevNode->next = nextNode;
    nextNode->previous = prevNode;

    while(currNode != nextNode) {
      Node<T>* delNode = currNode;
      currNode = currNode ->next;
      delete delNode;
    }
  }

public:
  LinkedList()
  {
    head = new Node<T>();
    tail = new Node<T>();
    head->next = tail;
    tail->previous = head;
  }

  ~LinkedList() 
  {
    Node<T>* current = head;
    while(current) {
        Node<T>* next = current->next;
        delete current;
        current = next;
    }
  }

  void push_back(const T& value) 
  {
    Node<T>* newNode = new Node<T>(value);
    Node<T>* prevNode = tail->previous;
    insert_between(prevNode, newNode, tail);
  }

  void push_front(const T& value) {
    Node<T>* newNode = new Node<T>(value);
    Node<T>* nextNode = head->next;
    insert_between(head, newNode, nextNode);
  }

  void pop_front() {
    Node<T>* prevNode = head;
    Node<T>* currNode = prevNode->next;
    if (currNode == tail) {
      return;
    }

    Node<T>* nextNode = currNode->next;

    delete_between(prevNode, currNode, nextNode);
  }
  
  void pop_back() {
    Node<T>* nextNode = tail;
    Node<T>* currNode = nextNode->previous;
    if (currNode == head) {
      return;
    }
    Node<T>* prevNode = currNode->previous;

    delete_between(prevNode, currNode, nextNode);

  }

  bool insert(unsigned position, const T& value) {
    Node<T>* newNode = new Node<T>(value);

    Node<T>* prevNode = head->next;
    if (prevNode == tail) {
      return false;
    }

    for (int i = 0; i < position; ++i) {
      prevNode = prevNode->next;
      if (prevNode == tail) {
        return false;
      }
    }

    Node<T>* nextNode = prevNode->next;
    insert_between(prevNode, newNode, nextNode);

    return true;
  }

  bool remove(unsigned position) {

    Node<T>* prevNode = head;
    Node<T>* currNode = head->next;

    if (currNode == tail) {
      return false;
    }

    for (int i = 0; i < position; ++i) {
      prevNode = currNode;
      currNode = currNode->next;
      if (currNode == tail) {
        return false;
      }
    }

    Node<T>* nextNode = currNode->next;
    delete_between(prevNode, currNode, nextNode);
    return true;
  }

  void display() const 
  {
    Node<T>* current = head->next;
    std::cout << "head -> ";
    while(current != tail) {
        std::cout << *(current->data) << " -> ";
        current = current->next;
    }
    std::cout << "nullptr\n";
  }
};

#endif // LINKED_LIST_H
