#include "linked_list.hpp"


int main() {
  LinkedList<int> list;

  list.push_back(4);
  list.push_back(5);
  list.push_back(6);

  list.push_front(2);
  list.push_front(1);

  list.insert(1,3);
  list.remove(5);

  list.pop_front();
  list.pop_back();
  list.pop_front();
  list.pop_front();
  list.pop_front();
  list.pop_front();
  list.pop_front();

  list.display();
}