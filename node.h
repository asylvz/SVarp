#ifndef __NODE
#define __NODE

using namespace std;

template<class T>
class Node
{
private:

public:
	T data;
	Node<T>* next;
};


template<class T>
class List {
  public:
    int empty();
    void add(T data);
    T remove();
    
	List() 
	{ 
      sentinel = new Node<T>(); 
      sentinel->next = sentinel; 
    }

  private:
    Node<T>* sentinel;
};

template<class T>
void List<T>::add(T data)
{
  Node<T> *p;

  p = new Node<T>();
  p->data = data;
  p->next = sentinel->next;
  sentinel->next = p;
}

template<class T>
T List<T>::remove()
{
  T data;
  Node<T> *node;

  if (sentinel->next == sentinel) {
    cerr << "ERROR: `remove' called with empty list.\n";
    exit(1);
  }

  node = sentinel->next;
  data = node->data;

  sentinel->next = node->next;
  delete node;

  return data;
}

#endif
