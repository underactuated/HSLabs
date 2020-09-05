using namespace std;

typedef float heapkey_t;

template <class Tc> struct heap;

typedef heap<double*> heapp_t;

template <class Tc>
struct heap{
  struct keyval{
    heapkey_t key;
    Tc* val;
    keyval (heapkey_t key_, Tc* val_){key=key_; val=val_;}
    bool operator < (const keyval& kv) const {
      return key < kv.key;
    }
  };
  heap();
  keyval& operator [] (int i) {
    return keyvals[i];
  }
  void push(heapkey_t key, Tc* val);
  void pop();
  heapkey_t key();
  Tc* val();
  int size();
  void sort();
  private:
    vector<keyval> keyvals;
};

template <class Tc>
heap<Tc>::heap(){
  make_heap(keyvals.begin(),keyvals.end());
}

template <class Tc>
void heap<Tc>::push(heapkey_t key, Tc* p){
  keyval kv (key, p);
  keyvals.push_back(kv);
  push_heap(keyvals.begin(),keyvals.end());
}

template <class Tc>
void heap<Tc>::pop(){
  pop_heap(keyvals.begin(),keyvals.end());
  keyvals.pop_back();
}

template <class Tc>
heapkey_t heap<Tc>::key(){
  return keyvals.front().key;
}

template <class Tc>
Tc* heap<Tc>::val(){
  return keyvals.front().val;
}

template <class Tc>
int heap<Tc>::size(){
  return keyvals.size();
}

template <class Tc> 
void heap<Tc>::sort(){
  sort_heap(keyvals.begin(),keyvals.end());
}

