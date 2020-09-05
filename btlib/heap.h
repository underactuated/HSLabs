using namespace std;

template <class T_k, class T_v>
struct heap{
  T_k key();
  T_v pull_val();
  void clear();
  void insert(T_k key, T_v val);
  int size();
  void print_keys(); // temp
  bool is_empty(){return !bool(keyvals.size());}
  struct comp{
    bool operator()(pair<T_k,T_v>& kv1, pair<T_k,T_v>& kv2){
      return (kv1.first < kv2.first);
    }
  };
  private:
    vector<pair<T_k, T_v> > keyvals;
};

template <class T_k, class T_v>
void heap<T_k,T_v>::insert(T_k key, T_v val){
  keyvals.push_back(pair<T_k,T_v> (key,val));
  push_heap(keyvals.begin(),keyvals.end(),comp());
}

template <class T_k, class T_v>
T_v heap<T_k,T_v>::pull_val(){
  pop_heap(keyvals.begin(),keyvals.end(),comp());
  T_v val = keyvals.back().second;
  keyvals.pop_back();
  return val;
}

template <class T_k, class T_v>
T_k heap<T_k,T_v>::key(){
  return keyvals.front().first;
}

template <class T_k, class T_v>
void heap<T_k,T_v>::clear(){
  keyvals.clear();
}

template <class T_k, class T_v>
int heap<T_k,T_v>::size(){
  return keyvals.size();
}

template <class T_k, class T_v>
void heap<T_k,T_v>::print_keys(){
  int n = size();
  for(int i=0;i<n;i++){
    if(i){cout << " ";}
    cout << keyvals[i].first;
  }
  cout << endl;
}

