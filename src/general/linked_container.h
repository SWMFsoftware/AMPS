

#include "stack.h"

#ifndef _LINKED_CONTAINER_H_
#define _LINKED_CONTAINER_H_

template <class T>
class cLinkedContainer {
public:

  class cNode {
  public:
    T data;
    cNode* next;
  };
  
  cStack<cNode> NodeStack;
  cNode *front_node;
  
  cLinkedContainer() {
    front_node=NULL;
    NodeStack.explicitConstructor();
  }
  
  ~cLinkedContainer() {
    NodeStack.clear();
  }
  
  class iterator {
  private:
    cNode *CurrentNode;
    
  public:
    iterator() {
      CurrentNode=NULL;
    }
    
    iterator(cNode* t) {
      CurrentNode=t;
    }
    
    iterator& operator=(cNode* pNode) {
      CurrentNode = pNode;
      return *this;
    }
    
    // Prefix ++ overload
    iterator& operator++() {
      if (CurrentNode!=NULL) CurrentNode = CurrentNode->next;
      return *this;
    }
    
    // Postfix ++ overload
    iterator operator++(int) {
      iterator it = *this;
      ++*this;
      return it;
    }
    
    bool operator!=(const iterator& it) {
      return CurrentNode != it.CurrentNode;
    }
    
    T operator*() {
      return CurrentNode->data;
    }
  };

  iterator end() {
    iterator it(NULL);
    
    return it;
  }
  
  iterator begin() {
    iterator it(front_node);
    
    return it;
  }
  
  void push_front(T& new_element) {
    cNode *t=NodeStack.newElement();

    t->data=new_element; 
    
    t->next=front_node;
    front_node=t;
  }
};




#endif


